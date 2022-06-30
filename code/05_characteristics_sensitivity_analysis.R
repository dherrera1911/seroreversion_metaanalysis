##################################
##################################
# 
# This script fits a model to data on time-varying
# sensitivity of different serology assays. Only data
# of serology testing on previously diagnosed individuals
# is used.
# 
# For the studies without known time from diagnosis to
# sero-testing, the diagnosis time was estimated (in
# script 05).
# 
# The script also does some fitting leaving out some of
# the data.
# 
##################################
##################################

library(dplyr)
library(tidyr)
library(lubridate)
library(rstan)
library(bayesplot)
library(tidybayes)
library(stringr)

source("./functions_auxiliary.R")
source("./functions_seroreversion_fit_analysis.R")

# Fitting parameters
nChains <- 4
nCores <- 2
nIter <- 4000

############
# 1) Load data generated in script 06
############
seroFitted <- read.csv("../data/analysis_results/PCR_to_serotest_all.csv",
                       stringsAsFactors=FALSE)

# Add binary columns indicating presence or absence of characteristics
seroFitted$spikeAntigen <- stringr::str_detect(seroFitted$antigenTarget, "S") |
  stringr::str_detect(seroFitted$antigenTarget, "RBD")
seroFitted$rbdAntigen <- stringr::str_detect(seroFitted$antigenTarget, "RBD")
seroFitted$nuceloAntigen <- stringr::str_detect(seroFitted$antigenTarget, "N")
seroFitted$LFA <- stringr::str_detect(seroFitted$technique, "LFIA")
seroFitted$ELISA <- stringr::str_detect(seroFitted$technique, "ELISA")
seroFitted$chemlum <- (stringr::str_detect(seroFitted$technique, "CMIA") |
                       stringr::str_detect(seroFitted$technique, "CLIA"))
seroFitted$IgM <- stringr::str_detect(seroFitted$antibodyTarget, "IgM")
seroFitted$IgA <- stringr::str_detect(seroFitted$antibodyTarget, "IgA")


############
# 2) Make matrix with assay characters to fit
############

### Filter rows without value for the characteristic of interest
nChars <- 4
seroFitted <- dplyr::filter(seroFitted, !is.na(seroFitted$antigenTarget))
assayVec <- as.factor(seroFitted$testName)
studies <- as.factor(paste(seroFitted$testName, seroFitted$citationID))
charMatrix <- matrix(as.integer(c(seroFitted$spikeAntigen, seroFitted$nuceloAntigen,
                                  seroFitted$rbdAntigen,
                                  seroFitted$spikeAntigen*seroFitted$nuceloAntigen,
                                  as.integer(assayVec))), ncol=nChars+1)
charsName <- c("S", "N", "RBD", "SN")
fileIdentifier <- "antigenTarget"

#nChars <- 2
#seroFitted <- dplyr::filter(seroFitted, !is.na(technique))
#assayVec <- as.factor(seroFitted$testName)
#studies <- as.factor(paste(seroFitted$testName, seroFitted$citationID))
#charMatrix <- matrix(as.integer(c(seroFitted$LFA, !seroFitted$LFA)), ncol=nChars)
#charsName <- c("LFA", "Rest")
#fileIdentifier <- "technique"
#

#nChars <- 4
#seroFitted <- dplyr::filter(seroFitted, !is.na(antibodyTarget))
#assayVec <- as.factor(seroFitted$testName)
#studies <- as.factor(paste(seroFitted$testName, seroFitted$citationID))
#charMatrix <- matrix(as.integer(c(rep(1, nrow(seroFitted)),
#                                  seroFitted$IgM, seroFitted$IgA,
#                                  seroFitted$IgM & seroFitted$IgA)), ncol=nChars)
#charsName <- c("IgG", "IgM", "IgA", "Total")
#fileIdentifier <- "antibody"
#

#nChars <- 2
#seroFitted <- dplyr::filter(seroFitted, !is.na(antibodyTarget))
#assayVec <- as.factor(seroFitted$testName)
#studies <- as.factor(paste(seroFitted$testName, seroFitted$citationID))
#charMatrix <- matrix(as.integer(c(rep(1, nrow(seroFitted)),
#                                  (seroFitted$IgM | seroFitted$IgA))), ncol=nChars)
#charsName <- c("IgG", "Others")
#fileIdentifier <- "antibody2"

#nChars <- 5
#seroFitted <- dplyr::filter(seroFitted, !(is.na(antigenTarget) & !is.na(technique)))
#assayVec <- as.factor(seroFitted$testName)
#studies <- as.factor(paste(seroFitted$testName, seroFitted$citationID))
#charMatrix <- matrix(as.integer(c(seroFitted$spikeAntigen, seroFitted$nuceloAntigen,
#                                  seroFitted$rbdAntigen,
#                                  seroFitted$spikeAntigen*seroFitted$nuceloAntigen,
#                                  seroFitted$LFA,
#                                  as.integer(assayVec))), ncol=nChars+1)
#charsName <- c("S", "N", "RBD", "SN", "LFA")
#fileIdentifier <- "fullModel"

# Shorten and sort character matrix to match index of assays
charMatrix <- unique(charMatrix)
charMatrix <- charMatrix[order(charMatrix[,nChars+1]),]
charMatrix <- charMatrix[,c(1:nChars)]


###################
# 3) Get initial values for the fitting
###################
firstGuess <- mean_slope_intercept(nSeropositive=seroFitted$nSeropositives,
                                   nSamples=seroFitted$nSamples,
                                   timeVec=seroFitted$testTime)

# Set the ranges for the initial values of the parameters
paramListName <- c("charSlope", "intercept", "slopeSigma", "interceptSigma",
  "studySigma", "assaySlope", "assayIntercept", "studyIntercept")
lowerUni <- c(firstGuess$slope, firstGuess$intercept, 0.2,
              firstGuess$sdIntercept*0.9*0.5,
              firstGuess$sdIntercept*0.9*0.5,
              firstGuess$slope,
              firstGuess$intercept, -0.1)
upperUni <- c(firstGuess$slope+0.05, firstGuess$intercept+1, 0.3,
              firstGuess$sdIntercept*1.1*0.5,
              firstGuess$sdIntercept*1.1*0.5,
              firstGuess$slope+0.05,
              firstGuess$intercept+1, 0.1)

nTests <- length(unique(seroFitted$testName))
nStudies <- length(unique(paste(seroFitted$testName, seroFitted$citationID)))
paramSize <- c(nChars, 1, 1, 1, 1, nTests, nTests, nStudies) # size of param vectors

# Sample the initial values to use
initList <- sample_initial_values(nChains=nChains, paramListName=paramListName,
                           lowerUni=lowerUni, upperUni=upperUni,
                           paramSize=paramSize)

###################
# 4) Fit the model to the data
###################

# compile model
outcome_reg <- rstan::stan_model("./sensitivity_change_assay_chars2.stan",
                                 model_name="time_change_sensitivity",
                                 warn_pedantic=TRUE)

# Make a list with the input we pass to STAN
assayDataList <- list(N=nrow(seroFitted),
                  K=length(unique(seroFitted$testName)),
                  M=length(unique(studies)),
                  C=nChars,
                  characteristics=charMatrix,
                  assay=as.integer(assayVec),
                  study=as.integer(studies),
                  timeVec=seroFitted$testTime/4,
                  nPositive=seroFitted$nSeropositives,
                  nTested=seroFitted$nSamples)

# Fit model
model <- rstan::sampling(outcome_reg, data=assayDataList,
                           chains=nChains, iter=nIter, refresh=0,
                           verbose=TRUE, cores=nCores, init=initList)

###################
# 5) Extract and save the posterior samples of the fit
###################
posteriorTraces <- tidybayes::gather_draws(model,
                                           intercept, slopeSigma,
                                           interceptSigma,
                                           studySigma,
                                           assayIntercept[loc],
                                           assaySlope[loc],
                                           studyIntercept[stu],
                                           charSlope[par])

posteriorTraces$parName <- charsName[posteriorTraces$par]
posteriorTraces$assay <- levels(assayVec)[posteriorTraces$loc]
studyList <- strsplit(as.character(levels(studies)), " n-")
assayVec <- NULL
studyVec <- NULL
for (i in c(1:length(studyList))) {
  assayVec <- c(assayVec, studyList[[i]][1])
  studyVec <- c(studyVec, studyList[[i]][2])
}
posteriorTraces$studyAssay <- assayVec[posteriorTraces$stu]
posteriorTraces$studyCitation <- studyVec[posteriorTraces$stu]

write.csv(posteriorTraces, paste("../data/analysis_results/7_assay_characteristics_",
                              fileIdentifier, ".csv", sep=""), row.names=FALSE)


####################
## 6) Fit the model to subsets of the data, to test performance on left out data
####################

# Tests with multiple time points, which we will use for testing
multiSamples <- group_by(seroFitted, testName) %>%
  summarize(., nTimes=length(unique(testTime)) > 1) %>%
  dplyr::filter(., nTimes)

# Set the testing groups
nLeaveOut <- 2
nGroups <- ceiling(nrow(multiSamples)/nLeaveOut)
valGroup <- rep(c(1:nGroups), nLeaveOut)[1:nrow(multiSamples)]
valGroup <- sample(valGroup)

validationDf <- NULL

for (ng in c(1:nGroups)) {
  leaveOutTests <- multiSamples$testName[valGroup==ng]
  meanTimes <- dplyr::filter(seroFitted, testName %in% leaveOutTests) %>%
    group_by(., testName) %>%
    summarize(., medianTime= mean(testTime))

  # Test on the first half of the assay data, and the second half separately
  for (testHalf in c(0, 1)) {
    # Create fitting and testing data frames
    valSampleDf <- NULL
    fittingDf <- seroFitted
    for (lo in c(1:nrow(meanTimes))) {
      tName <- meanTimes$testName[lo]
      timesVec <- fittingDf[fittingDf$testName==tName, "testTime"]
      if (testHalf == 0) {
        validationRows <- which((fittingDf$testName == tName &
                                  fittingDf$testTime <= meanTimes$medianTime[lo]))
      } else {
        validationRows <- which((fittingDf$testName == tName &
                                  fittingDf$testTime > meanTimes$medianTime[lo]))
      } 
      valSampleDf <- rbind(valSampleDf, fittingDf[validationRows,])
      fittingDf <- fittingDf[-validationRows,]
    }
    valSampleDf$validationTimeGroup <- testHalf

    # Prepare data and initializations for fitting data
    nStudies <- length(unique(paste(fittingDf$testName, fittingDf$citationID)))
    paramSize <- c(nChars, 1, 1, 1, 1, nTests, nTests, nStudies) # size of param vectors

    # Sample the initial values to use
    initList <- sample_initial_values(nChains=nChains,
                                      paramListName=paramListName,
                                      lowerUni=lowerUni, upperUni=upperUni,
                                      paramSize=paramSize)

    # Make a list with the input we pass to STAN
    assayVec <- as.factor(fittingDf$testName)
    studies <- as.factor(paste(fittingDf$testName, fittingDf$citationID))
    assayDataList <- list(N=nrow(fittingDf),
                      K=length(unique(fittingDf$testName)),
                      M=length(unique(studies)),
                      C=nChars,
                      characteristics=charMatrix,
                      assay=as.integer(assayVec),
                      study=as.integer(studies),
                      timeVec=fittingDf$testTime/4,
                      nPositive=fittingDf$nSeropositives,
                      nTested=fittingDf$nSamples)

    # Fit model
    modelVal <- rstan::sampling(outcome_reg, data=assayDataList,
                               chains=nChains, iter=nIter, refresh=0,
                               verbose=TRUE, cores=nCores, init=initList)

    posteriorTraces <- tidybayes::gather_draws(modelVal,
                                               intercept, slopeSigma,
                                               interceptSigma,
                                               studySigma,
                                               assayIntercept[loc],
                                               assaySlope[loc],
                                               studyIntercept[stu],
                                               charSlope[par])

    # Put assay labels to the samples of the assay (not any specific study)
    posteriorTraces$assay <- levels(assayVec)[posteriorTraces$loc]
    studyList <- strsplit(as.character(levels(studies)), " n-")
    # Put assay and study labels to the samples of specific studies
    studyAssayVec <- NULL
    studyVec <- NULL
    for (i in c(1:length(studyList))) {
      studyAssayVec <- c(studyAssayVec, studyList[[i]][1])
      studyVec <- c(studyVec, studyList[[i]][2])
    }
    posteriorTraces$studyAssay <- studyAssayVec[posteriorTraces$stu]
    posteriorTraces$studyCitation <- studyVec[posteriorTraces$stu]

    tempPred <- sensitivity_prediction_validation(posteriorTraces=posteriorTraces,
                                             validationDf=valSampleDf,
                                             binomSamples=5)
    validationDf <- rbind(validationDf, tempPred)
  }
}

write.csv(validationDf, paste("../data/analysis_results/7_assay_characteristics_",
                              fileIdentifier, "_validation.csv", sep=""), row.names=FALSE)

