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

######## Initial values for the fit

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

## Define function to make list of initial values, for STAN
sample_initial_values <- function(nChains, paramListName, lowerUni, upperUni, paramSize) {
  initList <- list()
  for (ch in c(1:nChains)) {
    initList[[ch]] <- list()
    for (p in c(1:length(paramListName))) {
      initList[[ch]][[paramListName[p]]] <- runif(paramSize[p], min=lowerUni[p],
                                                  max=upperUni[p])
    }
  }
  return(initList)
}


# Get posterior samples of sensitivity across time
seroreversion_samples <- function(posteriorDf, timeVec,
                               slopeName="timeSlope",
                               interceptName="intercept",
                               timeNormalization=4) {
  stdTimeVec <- timeVec/timeNormalization
  slopeVec <- dplyr::filter(posteriorDf, .variable==slopeName)[[".value"]]
  interceptVec <- dplyr::filter(posteriorDf, .variable==interceptName)[[".value"]]
  # fill a matrix with the sensitivities(time) samples as each row
  linMat <- slopeVec %*% t(stdTimeVec)
  linMat <- linMat + interceptVec %*% t(rep(1, length(stdTimeVec)))
  expMat <- exp(linMat)
  fitProp <- expMat / (1+expMat)
  # extract statistics on the samples
  meanProp <- colMeans(fitProp) 
  ciProp <- matrixStats::colQuantiles(fitProp, probs=c(0.025, 0.975))
  # put everything into the output list
  sampleVec <- sort(rep(c(1:nrow(fitProp)), ncol(fitProp)))
  timeVecLong <- rep(timeVec, nrow(fitProp))
  timeIndVec <- rep(c(1:length(timeVec)), nrow(fitProp))
  samplesDf <- data.frame(sample=sampleVec,
                          proportion=as.vector(fitProp),
                          time=timeVecLong, timeInd=timeIndVec)
  sampleList <- list(samples=samplesDf, prop_mean=meanProp,
                     prop_L=ciProp[,1], prop_H=ciProp[,2])
  return(sampleList)
}
  

# Make a data.frame containing the draws of the slopes for different
# combinations of assay parameters
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

# Get sensitivity across time for different levels of
# a parameter. posteriorDf has the samples of the parameters,
# timeVec gives the times at which to evaluate sensitivity.
# charCombs is a list containing character vectors that
# indicate what parameters to add up to return their sensitivity
# (e.g. an element c("S", "N") will return the sensitivities
# for tests with the slope_S + slope_N)
mean_param_seroreversion <- function(posteriorDf, timeVec,
                               slopeName="charSlope",
                               interceptName="intercept",
                               charCombs=NA,
                               timeNormalization=4) {
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
    ciProp <- matrixStats::colQuantiles(fitProp, probs=c(0.025, 0.975))
    tempDf <- data.frame(chars=paste(charCombs[[n]], collapse="-"),
                         time=timeVec,
                         sensitivity=meanProp,
                         sensitivityL=ciProp[,1],
                         sensitivityH=ciProp[,2])
    sensitivityDf <- rbind(sensitivityDf, tempDf)
    # This commented code makes a Df with all the sensitivity samples
    #sampleVec <- sort(rep(c(1:ncol(fitSampleMat)), nrow(fitSampleMat)))
    #timeVecLong <- rep(timeVec, ncol(fitSampleMat))
    #timeIndVec <- rep(c(1:length(timeVec)), ncol(fitSampleMat))
    #samplesDf <- data.frame(sample=sampleVec,
    #                        proportion=as.vector(fitSampleMat),
    #                        time=timeVecLong, timeInd=timeIndVec)
    #sampleList <- list(prop_mean=meanProp,
    #                   prop_L=ciProp[,1], prop_H=ciProp[,2])
  }
  return(sensitivityDf)
}


# Put together the assay-specific effect on slope with the
# effects of assay parameters, to get the overall slope
# for each assay at each posterior sample. Return a dataframe
# with overall assay slopes and intercepts
assay_char_slope <- function(posteriorDf, charCombs, assayRefDf) {
  # Get the slopes of relevant combinations of parameters
  charSlopes <- get_combination_slopes(posteriorDf, charCombs)
  # Get info on what assays to return
  assayVec <- unique(posteriorDf$assay)
  assayVec <- assayVec[!is.na(assayVec)]
  assayN <- length(assayVec)
  assayCharSlopesDf <- NULL
  for (a in c(1:assayN)) {
    assayName <- assayVec[a]
    # Extract assay characteristic from reference
    assayChar <- with(assayRefDf, chars[testName==assayName])
    assayDf <- dplyr::filter(posteriorDf, assay==assayName)
    assayDf$.value[assayDf$.variable=="assaySlope"] <-
      assayDf$.value[assayDf$.variable=="assaySlope"] + 
      charSlopes[[assayChar]]
    assayCharSlopesDf <- rbind(assayCharSlopesDf, assayDf)
  }
  return(assayCharSlopesDf)
}

assay_char_sensitivity_posterior <- function(posteriorDf, timeVec,
                                             charCombs,
                                             assayRefDf,
                                             timeNormalization=4) {
  # Put together assay slopes with characteristic effects
  assaysPosterior <- assay_char_slope(posteriorDf, charCombs, assayRefDf)
  # make some variables to fill/use
  assayVec <- sort(unique(assaysPosterior$assay))
  assayFitDf <- NULL
  for (n in c(1:length(assayVec))) {
    singleAssayPosterior <- dplyr::filter(assaysPosterior, assay==assayVec[n])
    sensitivityPosterior <- seroreversion_samples(singleAssayPosterior, timeVec,
                                                  slopeName="assaySlope",
                                                  interceptName="assayIntercept",
                                                  timeNormalization=4)
    assayChar <- assayRefDf$chars[assayRefDf$testName==assayVec[n]]
    testDf <- data.frame(time=timeVec,
                         sensitivity=sensitivityPosterior$prop_mean,
                         sensitivityL=sensitivityPosterior$prop_L,
                         sensitivityH=sensitivityPosterior$prop_H,
                         testName=assayVec[n],
                         chars=assayChar)
    assayFitDf <- rbind(assayFitDf, testDf)
  }
  return(assayFitDf)
}

# Take the posterior samples fitted on a subset of the data, and
# apply them to predict the results for another set of validation data
sensitivity_prediction_validation <- function(posteriorTraces,
                                              validationDf,
                                              binomSamples=5) {
  # Do one test at a time for efficiency
  valTests <- unique(validationDf$testName)
  nDraws <- max(posteriorTraces$.draw)
  validationDf$predictedPosL <- NA
  validationDf$predictedPosH <- NA
  validationDf$sensitivityMeanPred <- NA
  validationDf$sensitivityLPred <- NA
  validationDf$sensitivityHPred <- NA
  for (vt in c(1:length(valTests))) {
    # Filter posterior of this specific test
    testPosterior <- dplyr::filter(posteriorTraces, (assay==valTests[vt]) |
                                   (studyAssay==valTests[vt]))
    # extract test parameters
    assayInterceptVec <- dplyr::filter(testPosterior, .variable=="assayIntercept")$.value
    assaySlopeVec <- dplyr::filter(testPosterior, .variable=="assaySlope")$.value
    studySigmaVec <- dplyr::filter(posteriorTraces, .variable=="studySigma")$.value
    fittedStudies <- unique(testPosterior$studyCitation)[-1]
    valRows <- which(validationDf$testName==valTests[vt])
    for (vr in valRows) {
      # add a study effect to the assay intercept
      rowStudy <- strsplit(valSampleDf$citationID[vr], "-")[[1]][2]
      if (rowStudy %in% fittedStudies) {
        studyInterceptVec <- dplyr::filter(testPosterior, studyCitation==rowStudy)$.value
        totalIntercept <- assayInterceptVec + studyInterceptVec
      } else {
        studyInterceptVec <- rnorm(nDraws) * studySigmaVec
        totalIntercept <- assayInterceptVec + studyInterceptVec
      }
      # Get the estimated vector of sensitivities. Store the mean and 95% ci
      linVec <- totalIntercept + assaySlopeVec * validationDf$testTime[vr]
      sensitivityVec <- exp(linVec)/(1+exp(linVec))
      validationDf$sensitivityMeanPred[vr] <- mean(sensitivityVec)*100
      sensitivityQuants <- quantile(sensitivityVec, probs=c(0.025, 0.975))
      validationDf$sensitivityLPred[vr] <- sensitivityQuants[1]*100
      validationDf$sensitivityHPred[vr] <- sensitivityQuants[2]*100
      # Sample the original data (nSeropositives) from the reported sensitivities, and store
      samplePositives <- NULL
      for (ss in c(1:binomSamples)) {
        samplePositives <- c(samplePositives,
                             rbinom(n=nDraws, size=valSampleDf$nSamples[vr],
                                    prob=sensitivityVec))
      }
      ciRow <- quantile(samplePositives, probs=c(0.025, 0.975))
      validationDf$predictedPosL[vr] <- ciRow[1]
      validationDf$predictedPosH[vr] <- ciRow[2]
    }
  }
  return(validationDf)
}


