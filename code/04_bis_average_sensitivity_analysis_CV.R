###################################
# 
# This script fits a model to data on time-varying
# sensitivity of different serology assays.
#
# The model is fitted without accounting for test characteristics
#
# This script does the same analysis as
# script 04_average_sensitivity_analysis.R, but separating the
# dataset into groups, and performing cross-validation.
# The results of this script are the reported percentages
# of the cross-validation performance in the associated paper,
# for section 1 of the Results. They are also shown in Fig S1.
#
# The output is a csv file
# "../data/analysis_results/04_predicted_sensitivities_grouped_CV.csv"
# with a similar structure to the # input data file "PCR_to_serotest_all.csv",
# but with additional columns that contain the predicted sensitivities,
# as well as the predicted experimental outcomes (nPositives, which take
# into account the binomial sampling error too).
#
# Script authored by Daniel Herrera-Esposito.
# For questions, contact me at dherrera1911[at]gmail.com
# 
# Final version revised 10/03/2023
# 
##################################

library(dplyr)
library(tidyr)
library(lubridate)
library(rstan)
library(bayesplot)
library(tidybayes)

source("./functions_auxiliary.R")
source("./functions_seroreversion_fit_analysis.R")

### Fitting parameters
nChains <- 4
nCores <- 2
nIter <- 4000
warmup <- 1000
nChar <- 1
# "grouped" is the cross-validation procedure described in the associated paper
crossValidationType <- "grouped" 
testsPerGroup <- 1 # Don't mind this parameter. When it's larger than 1 it
# saves time in the cross validation, but gives less precise results. It was
# used during debugging.

#crossValidationType <- "random"
#pointsPerGroup <- 7

seroFitted <- read.csv("../data/processed_data/PCR_to_serotest_all.csv",
                       stringsAsFactors=FALSE)
### Compile model
outcome_reg <- rstan::stan_model("./sensitivity_change_chars.stan",
                                 model_name="time_change_sensitivity",
                                 warn_pedantic=TRUE)
# load model if available
#outcome_reg <- readRDS("./sensitivity_change_chars.RDS")


### Make the cross validation batch groups
if (crossValidationType == "grouped") {
  # GROUPED
  #### Make index indicating data in each validation group
  # Validation data are grouped, so that either the first or
  # second half of the data for a given test are removed at
  # once, so CV involves extrapolation.
  validationGroup <- grouped_val_batches(seroFitted, testsPerGroup)
  # attach CV indexing to main data frame
  seroFitted$validationGroup <- validationGroup
} else if (crossValidationType=="random") {
  # RANDOM
  # Just divide the rows randomly into batches with the
  # number of points indicated above
  nGroups <- ceiling(nrow(seroFitted)/pointsPerGroup)
  validationGroup <- rep(c(1:nGroups), pointsPerGroup)[1:nrow(seroFitted)]
  validationGroup <- sample(validationGroup)
}


### Set some initialization parameters for the fits
firstGuess <- mean_slope_intercept(nSeropositive=seroFitted$nSeropositives,
                                   nSamples=seroFitted$nSamples,
                                   timeVec=seroFitted$testTime)
paramListName <- c("charSlope", "intercept", "slopeSigma", "interceptSigma",
  "studySigma", "assaySlope", "assayIntercept", "studyIntercept")
lowerBound <- c(firstGuess$slope, firstGuess$intercept, 0.2,
              firstGuess$sdIntercept*0.9*0.5,
              firstGuess$sdIntercept*0.9*0.5,
              firstGuess$slope,
              firstGuess$intercept, -0.1)
upperBound <- c(firstGuess$slope+0.05, firstGuess$intercept+1, 0.3,
              firstGuess$sdIntercept*1.1*0.5,
              firstGuess$sdIntercept*1.1*0.5,
              firstGuess$slope+0.05,
              firstGuess$intercept+1, 0.1)

### Fit and validate
nValGroups <- length(unique(validationGroup))
allValidationDf <- NULL
for (vg in c(1:nValGroups)) {
  fitDf <- dplyr::filter(seroFitted, validationGroup!=vg)
  valDf <- dplyr::filter(seroFitted, validationGroup==vg)

  # Make initialization values
  nTests <- length(unique(fitDf$testName))
  nStudies <- length(unique(paste(fitDf$testName, fitDf$citationID)))
  paramSize <- c(1, 1, 1, 1, 1, nTests, nTests, nStudies) # size of param vectors
  initList <- sample_initial_values(nChains=nChains, paramListName=paramListName,
                             lowerBound=lowerBound, upperBound=upperBound,
                             paramSize=paramSize)

  # Make data list to input to STAN
  charMatrix <- as.matrix(rep(1, length(unique(fitDf$testName))), ncol=1)
  assayVec <- as.factor(fitDf$testName)
  studies <- as.factor(paste(fitDf$testName, fitDf$citationID))
  assayDataList <- list(N=nrow(fitDf),
                    K=length(unique(fitDf$testName)),
                    M=length(unique(studies)),
                    C=nChar,
                    characteristics=charMatrix,
                    assay=as.integer(assayVec),
                    study=as.integer(studies),
                    timeVec=fitDf$testTime,
                    nPositive=fitDf$nSeropositives,
                    nTested=fitDf$nSamples)

  # Fit model
  model <- rstan::sampling(outcome_reg, data=assayDataList,
                             chains=nChains, iter=nIter, refresh=0,
                             warmup=warmup,
                             verbose=TRUE, cores=nCores, init=initList)

  # Extract and tidy posterior
  posteriorTraces <- tidybayes::gather_draws(model, 
                                             intercept, slopeSigma,
                                             interceptSigma,
                                             studySigma,
                                             assayIntercept[loc],
                                             assaySlope[loc],
                                             studyIntercept[stu],
                                             charSlope[par])

  posteriorTraces$assay <- levels(assayVec)[posteriorTraces$loc]
  studyList <- strsplit(as.character(levels(studies)), " n-")
  studyAssayVec <- NULL
  studyVec <- NULL
  for (i in c(1:length(studyList))) {
    studyAssayVec <- c(studyAssayVec, studyList[[i]][1])
    studyVec <- c(studyVec, studyList[[i]][2])
  }
  posteriorTraces$studyAssay <- studyAssayVec[posteriorTraces$stu]
  posteriorTraces$studyCitation <- studyVec[posteriorTraces$stu]

  # Predict the outcome of the left out samples, and append to validation Df
  valOutput <- sensitivity_prediction_validation(posteriorTraces, valDf,
    binomSamples=5, timeNormalization=1)
  allValidationDf <- rbind(valOutput, allValidationDf)
}

valName <- paste("../data/analysis_results/04_predicted_sensitivities_",
                 crossValidationType, "_CV.csv", sep="")
write.csv(allValidationDf, valName, row.names=FALSE)

