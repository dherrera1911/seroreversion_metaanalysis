##################################
##################################
# 
# This script fits a model to data on time-varying
# sensitivity of different serology assays. Only data
# of serology testing on previously diagnosed individuals
# is used.
#
# The model is fitted without accounting for the assays technical
# characteristics.
# 
##################################
##################################

library(dplyr)
library(tidyr)
library(lubridate)
library(rstan)
library(bayesplot)
library(tidybayes)

source("./functions_auxiliary.R")
source("./functions_seroreversion_fit_analysis.R")

# Fitting parameters
nChains <- 4
nCores <- 2
nIter <- 4000
warmup <- 1000

############
# Load data generated in script 06
############
seroFitted <- read.csv("../data/processed_data/PCR_to_serotest_all.csv",
                       stringsAsFactors=FALSE)

# Same slope shared by all assays
nChar <- 1
charMatrix <- as.matrix(rep(1, length(unique(seroFitted$testName))), ncol=1)

###################
# Get initial values for the fitting
###################
firstGuess <- mean_slope_intercept(nSeropositive=seroFitted$nSeropositives,
                                   nSamples=seroFitted$nSamples,
                                   timeVec=seroFitted$testTime)

# Set the ranges for the initial values of the parameters
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

nTests <- length(unique(seroFitted$testName))
nStudies <- length(unique(paste(seroFitted$testName, seroFitted$citationID)))
paramSize <- c(1, 1, 1, 1, 1, nTests, nTests, nStudies) # size of param vectors

# Sample the initial values to use
initList <- sample_initial_values(nChains=nChains, paramListName=paramListName,
                           lowerBound=lowerBound, upperBound=upperBound,
                           paramSize=paramSize)

###################
# Fit the model to the data
###################

# compile model
#regression_model <- rstan::stan_model("./sensitivity_change_chars.stan",
#                                 model_name="time_change_sensitivity",
#                                 warn_pedantic=TRUE)
# load model if available
regression_model <- readRDS("./sensitivity_change_chars.RDS")

# Make a list with the input we pass to STAN
assayVec <- as.factor(seroFitted$testName)
studies <- as.factor(paste(seroFitted$testName, seroFitted$citationID))
assayDataList <- list(N=nrow(seroFitted),
                  K=length(unique(seroFitted$testName)),
                  M=length(unique(studies)),
                  C=nChar,
                  characteristics=charMatrix,
                  assay=as.integer(assayVec),
                  study=as.integer(studies),
                  timeVec=seroFitted$testTime,
                  nPositive=seroFitted$nSeropositives,
                  nTested=seroFitted$nSamples)

# Fit model
fittedModel <- rstan::sampling(regression_model, data=assayDataList,
                           chains=nChains, iter=nIter, refresh=0,
                           warmup=warmup,
                           verbose=TRUE, cores=nCores, init=initList)

###################
# Extract and save the posterior samples of the fit
###################
posteriorTraces <- tidybayes::gather_draws(fittedModel, 
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

write.csv(posteriorTraces, "../data/analysis_results/04_model_posteriors_samples.csv",
          row.names=FALSE)


###################
# Save the predicted sensitivity profile for each test,
# and the test specific parameters
###################

# Get the samples of the sensitivity across time posterior for each test
testNames <- levels(assayVec)
timeVec <- seq(0.5, 14, 0.1)
assayFitDf <- NULL
summaryDf <- NULL
for (tN in testNames) {
  # get posterior 
  testTrace <- dplyr::filter(posteriorTraces, assay==tN)
  testProfileDf <- seroreversion_samples(testTrace, timeVec,
                                  slopeName="assaySlope",
                                  interceptName="assayIntercept",
                                  timeNormalization=1)
  testProfileDf$testName <- tN
  assayFitDf <- rbind(assayFitDf, testProfileDf)
}
# Get sensitivity across time for the average test
averagePosterior <- dplyr::filter(posteriorTraces, is.na(assay))
averageDf <- seroreversion_samples(averagePosterior, timeVec,
                                   slopeName="charSlope",
                                   interceptName="intercept",
                                   timeNormalization=1)
averageDf$testName <- NA
assayFitDf <- rbind(averageDf, assayFitDf)
# Export sensitivity profiles
write.csv(assayFitDf, "../data/analysis_results/04_assay_sensitivity_curve.csv",
          row.names=FALSE)

# Get summary statistics of the fitted parameters
parameterSummary <- ungroup(posteriorTraces) %>%
  group_by(., .variable, assay) %>%
  summarize(., paramMean=mean(.value), paramMedian=median(.value),
            paramSD=sd(.value),
            paramL=quantile(.value, probs=0.025),
            paramQuartileL=quantile(.value, probs=0.25),
            paramQuartileH=quantile(.value, probs=0.75),
            paramH=quantile(.value, probs=0.975))

write.csv(parameterSummary, "../data/analysis_results/04_parameter_summary.csv",
          row.names=FALSE)

