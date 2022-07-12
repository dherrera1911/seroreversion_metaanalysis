library(dplyr)
library(lubridate)
library(stringr)
library(ggplot2)
library(ggpubr)
library(matrixStats)
library(tidybayes)
library(rstan)

####################
######## Auxiliary functions to analyize bayesian fits to seroreversion data
####################

### Functions to facilitate obtaining initial values of MCMC

# Basic logistic regression fit, for initial slope and mean vals
mean_slope_intercept <- function(nSeropositives, nSamples, timeVec) {
  sensitivityVec <- nSeropositives/nSamples
  # Get logit values of sensitivities for initial intercept and slope
  logitVals <- log((sensitivityVec/100)/(1-sensitivityVec/100))
  # mean logits at time 1
  meanIntercept <- mean(logitVals[timeVec<=1 & logitVals!=Inf])
  sdIntercept <- sd(logitVals[timeVec<=1 & logitVals!=Inf])
  # slopes
  normalizedTime <- timeVec/mean(timeVec)
  # fit glm and get slope
  logReg <- glm(cbind(nSeropositives, nSamples-nSeropositives) ~ normalizedTime, family="binomial")
  meanSlope <- logReg$coefficients[2]
  return(list(slope=meanSlope, intercept=meanIntercept, sdIntercept=sdIntercept))
}

# Sample the initial values for a STAN fit, given input instructions.
#
# nChains: number of chains in the model
# paramListName: list with the names of the parameters # in the model
# lowerBound, uperBound: vectors with lower (and upper) bounds of uniform
#   intervals for sampling each parameter in paramListName
# paramSize: size of vector that is each parameter. is 1 for scalar parameters
sample_initial_values <- function(nChains, paramListName, lowerBound, upperBound,
                                  paramSize) {
  initList <- list()
  for (ch in c(1:nChains)) {
    initList[[ch]] <- list()
    for (p in c(1:length(paramListName))) {
      initList[[ch]][[paramListName[p]]] <- runif(paramSize[p], min=lowerBound[p],
                                                  max=upperBound[p])
    }
  }
  return(initList)
}

# Take values and apply the inverse logit function to them
inverse_logit <- function(x) {
  expPart <- exp(x)
  invLogit <- expPart/(1+expPart)
  return(invLogit)
}

### Functions to extract or process the posterior samples of a STAN fit

# Generate a sensitivity across time profile from the posterior
# samples of the model parameters.
#
# posteriorDf: Dataframe of posterior samples, obtained with function
#   tidybayes::gather_draws()
# timeVec: Time points at which to evaluate sensitivity
# slopeName: Name that the slope parameter takes in posteriorDf
# interceptName: Name that the intercept parameter takes in posteriorDf
# timeNormalization: Constant used to normalize time in model fitting
seroreversion_samples <- function(posteriorDf, timeVec,
                               slopeName="timeSlope",
                               interceptName="intercept",
                               timeNormalization=1) {
  stdTimeVec <- timeVec/timeNormalization
  slopeVec <- dplyr::filter(posteriorDf, .variable==slopeName)[[".value"]]
  interceptVec <- dplyr::filter(posteriorDf, .variable==interceptName)[[".value"]]
  # fill a matrix with the sensitivities(time) samples as each row
  linMat <- slopeVec %*% t(stdTimeVec)
  linMat <- linMat + interceptVec %*% t(rep(1, length(stdTimeVec)))
  fitProp <- inverse_logit(linMat)
  # extract statistics on the samples
  meanProp <- colMeans(fitProp) 
  ciProp <- matrixStats::colQuantiles(fitProp, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
  # put into dataFrame and export
  sensitivityDf <- data.frame(time=timeVec,
                     sensitivityMean=meanProp,
                     sensitivityMedian=ciProp[,3],
                     sensitivityL=ciProp[,1],
                     sensitivityQL=ciProp[,2],
                     sensitivityQH=ciProp[,4],
                     sensitivityH=ciProp[,5])
  return(sensitivityDf)
}

# Get the slopes resulting from the combination of different
# effects fitted in a model. E.g. the slopes for tests that
# target N antigen, S antigen, or the combination of the two.
#
# posteriorDf: Dataframe of posterior samples, obtained with function
#   tidybayes::gather_draws()
# charCombs: List where each element has the names of the effects
#   that are to be combined.
#   E.g. charCombs <- list(c("S"), c("N"), c("S", "N"))
#   returns 3 slopes: assays targeting S alone, N alone, and S and N. 
get_combination_slopes <- function(posteriorDf, charCombs) {
  parDf <- dplyr::select(posteriorDf, .draw, parName, .value) %>%
    dplyr::filter(., !is.na(parName)) %>%
    tidyr::pivot_wider(., names_from="parName", values_from=".value")
  charSlopes <- data.frame(.draw=parDf$.draw)
  for (cc in charCombs) {
    slopeVec <- dplyr::select(parDf, cc) %>%
      as.matrix() %>%
      rowSums(.)
    tempDf <- data.frame(tempName=slopeVec) 
    names(tempDf) <- paste(cc, collapse="-")
    charSlopes <- cbind(charSlopes, tempDf)
  }
  return(charSlopes)
}

# Get the average sensitivity profile for different combinations
# of effects in a fit.
# E.g. the sensitivity profiles for tests that
# target N antigen, S antigen, or the combination of the two.
#
# posteriorDf: Dataframe of posterior samples, obtained with function
#   tidybayes::gather_draws()
# timeVec: Times at which to evaluate sensitivity
# interceptName: Name that the intercept parameter has in the model
# charCombs: List where each element has the names of the effects
#   that are to be combined.
#   E.g. charCombs <- list(c("S"), c("N"), c("S", "N"))
#   returns 3 slopes: assays targeting S alone, N alone, and S and N. 
# timeNormalization: constant used to normalize time in the fit
mean_param_seroreversion <- function(posteriorDf, timeVec,
                               interceptName="intercept",
                               charCombs=NA,
                               timeNormalization=1) {
  stdTimeVec <- timeVec/timeNormalization
  chars <- unique(posteriorDf$parName)
  chars <- chars[!is.na(chars)]
  if (is.na(charCombs)) {charCombs <- chars}
  sensitivityDf <- NULL
  # Get the slopes of combinations of parameters
  combinationSlopes <- get_combination_slopes(posteriorDf, charCombs)
  # extract the shared intercepts
  intercept <- dplyr::filter(posteriorDf, .variable==interceptName)$.value
  for (n in c(1:length(charCombs))) {
    combName <- paste(charCombs[[n]], collapse="-")
    # compute the slope samples of this particular combination
    charSlope <- combinationSlopes[[combName]]
    # fill a matrix with the sensitivities(time) samples as each column
    linMat <- charSlope %*% t(stdTimeVec)
    linMat <- linMat + intercept %*% t(rep(1, length(stdTimeVec)))
    expMat <- exp(linMat)
    fitProp <- expMat / (1+expMat)
    # get the mean and 95% CI sensitivities
    meanProp <- colMeans(fitProp) 
    ciProp <- matrixStats::colQuantiles(fitProp, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
    tempDf <- data.frame(chars=paste(charCombs[[n]], collapse="-"),
                         time=timeVec,
                         sensitivityMean=meanProp,
                         sensitivityMedian=ciProp[,3],
                         sensitivityL=ciProp[,1],
                         sensitivityQL=ciProp[,2],
                         sensitivityQH=ciProp[,4],
                         sensitivityH=ciProp[,5])
    sensitivityDf <- rbind(sensitivityDf, tempDf)
  }
  return(sensitivityDf)
}


# Put together the assay-specific effect on slope with the
# effects of assay parameters, to get the overall slope
# for each assay at each posterior sample. Return a dataframe
# with overall assay slopes and intercepts
#assay_char_slope <- function(posteriorDf, charCombs, assayRefDf) {
#  # Get the slopes of relevant combinations of parameters
#  charSlopes <- get_combination_slopes(posteriorDf, charCombs)
#  # Get info on what assays to return
#  assayVec <- unique(posteriorDf$assay)
#  assayVec <- assayVec[!is.na(assayVec)]
#  assayN <- length(assayVec)
#  assayCharSlopesDf <- NULL
#  for (a in c(1:assayN)) {
#    assayName <- assayVec[a]
#    # Extract assay characteristic from reference
#    assayChar <- with(assayRefDf, chars[testName==assayName])
#    assayDf <- dplyr::filter(posteriorDf, assay==assayName)
#    assayDf$.value[assayDf$.variable=="assaySlope"] <-
#      assayDf$.value[assayDf$.variable=="assaySlope"] + 
#      charSlopes[[assayChar]]
#    assayCharSlopesDf <- rbind(assayCharSlopesDf, assayDf)
#  }
#  return(assayCharSlopesDf)
#}
#
#
#assay_char_sensitivity_posterior <- function(posteriorDf, timeVec,
#                                             charCombs,
#                                             assayRefDf,
#                                             timeNormalization=4) {
#  # Put together assay slopes with characteristic effects
#  assaysPosterior <- assay_char_slope(posteriorDf, charCombs, assayRefDf)
#  # make some variables to fill/use
#  assayVec <- sort(unique(assaysPosterior$assay))
#  assayFitDf <- NULL
#  for (n in c(1:length(assayVec))) {
#    singleAssayPosterior <- dplyr::filter(assaysPosterior, assay==assayVec[n])
#    sensitivityPosterior <- seroreversion_samples(singleAssayPosterior, timeVec,
#                                                  slopeName="assaySlope",
#                                                  interceptName="assayIntercept",
#                                                  timeNormalization=4)
#    assayChar <- assayRefDf$chars[assayRefDf$testName==assayVec[n]]
#    testDf <- data.frame(time=timeVec,
#                         sensitivity=sensitivityPosterior$prop_mean,
#                         sensitivityL=sensitivityPosterior$prop_L,
#                         sensitivityH=sensitivityPosterior$prop_H,
#                         testName=assayVec[n],
#                         chars=assayChar)
#    assayFitDf <- rbind(assayFitDf, testDf)
#  }
#  return(assayFitDf)
#}


# Take the posterior samples fitted on a subset of the data, and
# apply them to predict the results for another set of validation data.
# Takes into account the binomial error in the sample to be predicted
#
# posteriorTraces: Posterior samples of the parameters in the STAN model,
#   obtained with tidybayes::gather_draws()
# validationDf: Dataframe with validation samples, to predict with the model
# binomSamples: Number of times that the binomial outcome of positive tests
#   is sampled for each posterior sample. I.e. number of binomial
#   draws of the validation datapoint for each MCMC sample.
sensitivity_prediction_validation <- function(posteriorTraces,
                                              validationDf,
                                              binomSamples=5,
                                              timeNormalization=1) {
  # Extract some useful info
  valTests <- unique(validationDf$testName)
  nDraws <- max(posteriorTraces$.draw)
  # Extract parameters general to test population (most used for assays that 
  # are in validation Df but not in fitted Df)
  # Mean intercept and slope values
  interceptVec <- dplyr::filter(posteriorTraces, .variable=="intercept")$.value
  slopeDf <- dplyr::filter(posteriorTraces, .variable=="charSlope")
  multiChar <- max(slopeDf$par)>1
  if (multiChar) {
    allChars <- unique(slopeDf$parName)
  }
  # SD of between-study variability
  studySigmaVec <- dplyr::filter(posteriorTraces, .variable=="studySigma")$.value
  # SD of assay intercept variability
  interceptSigmaVec <- dplyr::filter(posteriorTraces, .variable=="interceptSigma")$.value
  # SD of assay slope variability
  slopeSigmaVec <- dplyr::filter(posteriorTraces, .variable=="slopeSigma")$.value
  # Add rows in the validation dataframe, to fill with predictions
  validationDf$predictedPos <- NA
  validationDf$predictedPosL <- NA
  validationDf$predictedPosH <- NA
  validationDf$sensitivityMeanPred <- NA
  validationDf$sensitivityMedianPred <- NA
  validationDf$sensitivityLPred <- NA
  validationDf$sensitivityHPred <- NA
  # Go one test at a time
  for (vt in c(1:length(valTests))) {
    # Either extract the assay-specific samples from the fit, if the assay
    # is in both validation and fitted datasets, or generate the samples
    # from the across-assay parameters, if it wasn't.
    #
    # Rows in the validation dataframe to use/fill for this test
    valRows <- which(validationDf$testName==valTests[vt])
    if (valTests[vt] %in% unique(posteriorTraces$assay)) {
      # If this assay is among the fitted ones, EXTRACT assay-specific samples
      # Filter posterior of this specific assay
      testPosterior <- dplyr::filter(posteriorTraces, (assay==valTests[vt]) |
                                     (studyAssay==valTests[vt]))
      # Intercept of the assay
      assayInterceptVec <- dplyr::filter(testPosterior, .variable=="assayIntercept")$.value
      # Slope of the assay
      assaySlopeVec <- dplyr::filter(testPosterior, .variable=="assaySlope")$.value
      # What studies were fitted for this assay
      fittedStudies <- unique(testPosterior$studyCitation)[-1]
    } else {
      # If this assay is NOT among the fitted ones, GENERATE assay-specific samples
      # Intercept of the assay
      assayInterceptVec <- rnorm(nDraws) * interceptSigmaVec + interceptVec
      fittedStudies <- NULL
      # Slope of the assay. Have to put together characteristic effects on the slope
      if (!multiChar) {
        assaySlopeVec <- slopeDf$.value + rnorm(nDraws) * slopeSigmaVec 
      } else {
        # extract slope fixed effect
        charLogicVec <- which(as.logical(validationDf[valRows[1],][allChars]))
        assayChars <- allChars[charLogicVec]
        assaySlope <- get_combination_slopes(posteriorDf=ungroup(slopeDf),
                                             charCombs=list(assayChars))
        # add fixed effect and random effect
        assaySlopeVec <- as.numeric(assaySlope[,2]) + rnorm(nDraws) * slopeSigmaVec
      }
    }
    # Go row by row, generating the predictions
    for (vr in valRows) {
      # Add the effect of this row's study (i.e. between-study variability)
      # to the intercept
      rowStudy <- strsplit(validationDf$citationID[vr], "-")[[1]][2]
      if (rowStudy %in% fittedStudies) {
        studyInterceptVec <- dplyr::filter(testPosterior, studyCitation==rowStudy)$.value
        totalIntercept <- assayInterceptVec + studyInterceptVec
      } else {
        studyInterceptVec <- rnorm(nDraws) * studySigmaVec
        totalIntercept <- assayInterceptVec + studyInterceptVec
      }
      # Now that we have the final vector of intercept samples and
      # slope samples for this specific row, estimate the sensitivities
      # at this time. Extract mean, 95% CI
      linVec <- totalIntercept + assaySlopeVec * validationDf$testTime[vr] / timeNormalization
      sensitivityVec <- exp(linVec)/(1+exp(linVec))
      validationDf$sensitivityMeanPred[vr] <- mean(sensitivityVec)*100
      validationDf$sensitivityMedianPred[vr] <- median(sensitivityVec)*100
      sensitivityQuants <- quantile(sensitivityVec, probs=c(0.025, 0.975))
      validationDf$sensitivityLPred[vr] <- sensitivityQuants[1]*100
      validationDf$sensitivityHPred[vr] <- sensitivityQuants[2]*100
      # Use the sensitivity to sample the expected number of positive samples
      # and 95% CI among the tested samples
      samplePositives <- NULL
      for (ss in c(1:binomSamples)) {
        samplePositives <- c(samplePositives,
                             rbinom(n=nDraws, size=validationDf$nSamples[vr],
                                    prob=sensitivityVec))
      }
      ciRow <- quantile(samplePositives, probs=c(0.025, 0.5, 0.975))
      validationDf$predictedPosL[vr] <- ciRow[1]
      validationDf$predictedPos[vr] <- ciRow[2]
      validationDf$predictedPosH[vr] <- ciRow[3]
    }
  }
  return(validationDf)
}

