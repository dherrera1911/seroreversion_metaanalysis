##################################
# 
# This script fits a model to the subset of data
# that are RBD-targeting and not LFA assays, and fits
# two slopes: a slope for the first 3 months after diagnosis
# (or as determined in slopeCutoff), and a slope after the
# first 3 months.
# 
# The structure of the script and the outputs is analogous
# to scripts 04 and 05. The results obtained in this script are
# the ones presented in Section Supplementary D and Figure S2.
#
# The main reported results are the parameter summaries that are
# saved into "../data/analysis_results/06_parameter_summary.csv"
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
library(stringr)

source("./functions_auxiliary.R")
source("./functions_seroreversion_fit_analysis.R")

# Fitting parameters
nChains <- 4
nCores <- 2
nIter <- 4000
warmup <- 1000

slopeCutoff <- 3

############
# Load data generated in script 06
############
seroFitted <- read.csv("../data/processed_data/PCR_to_serotest_all.csv",
                       stringsAsFactors=FALSE)

filterCols <- c("RBD", "Rest")

############
# Add binary columns to dataframe, to work as indicator variables
############
seroFitted$RBD <- stringr::str_detect(seroFitted$antigenTarget, "RBD")
seroFitted$LFA <- stringr::str_detect(seroFitted$technique, "LFIA")
seroFitted$Rest <- !seroFitted$LFA

for (fc in filterCols) {
  seroFitted <- seroFitted[which(seroFitted[[fc]]),]
}

remainingTests <- unique(seroFitted$testName)

# Filter those that don't have times at both sides of the cutoff
for (t in c(1:length(remainingTests))) {
  testTimes <- seroFitted$testTime[seroFitted$testName == remainingTests[t]]
  if (!(any(testTimes<=3) & any(testTimes>3))) {
    seroFitted <- dplyr::filter(seroFitted, !(testName==remainingTests[t]))
  }
}

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
paramListName <- c("charSlope", "meanSlopeLate", "intercept",
                   "slopeSigmaEarly", "slopeSigmaLate",
                   "interceptSigma", "studySigma",
                   "assaySlopeEarly", "assaySlopeLate", "assayIntercept",
                   "studyIntercept")

interceptFloor <- -2
interceptCeil <- -1
sdFloor <- 0.4
sdCeil <- 0.5
sdFloorSlope <- 0.2
sdCeilSlope <- 0.3
lowerBound <- c(firstGuess$slope, firstGuess$slope, interceptFloor,
                sdFloorSlope, sdFloorSlope,
                sdFloor, sdFloor,
                firstGuess$slope, firstGuess$slope, interceptFloor, -0.1)
upperBound <- c(firstGuess$slope+0.05, firstGuess$slope+0.05, interceptCeil,
                sdCeilSlope, sdCeilSlope,
                sdCeil, sdCeil,
                firstGuess$slope+0.05, firstGuess$slope+0.05, interceptCeil, 0.1)

nTests <- length(unique(seroFitted$testName))
nStudies <- length(unique(paste(seroFitted$testName, seroFitted$citationID)))
paramSize <- c(nChar, 1, 1, 1, 1, 1, 1, nTests, nTests, nTests, nStudies) # size of param vectors

# Sample the initial values to use
initList <- sample_initial_values(nChains=nChains, paramListName=paramListName,
                           lowerBound=lowerBound, upperBound=upperBound,
                           paramSize=paramSize)


###################
# Fit the model to the data
###################

# compile model
regression_model <- rstan::stan_model("./stan_models/sensitivity_change_lateSlope.stan",
                                 model_name="time_change_sensitivity",
                                 warn_pedantic=TRUE)
# load model if available
#regression_model <- readRDS("./stan_models/sensitivity_change_lateSlope.RDS")

# Make a list with the input we pass to STAN
# Separate times into early time vector and late time vector
earlyTime <- seroFitted$testTime
earlyTime[earlyTime>slopeCutoff] <- slopeCutoff
lateTime <- seroFitted$testTime
lateTime[earlyTime<slopeCutoff] <- slopeCutoff
lateTime <- lateTime - slopeCutoff

assayVec <- as.factor(seroFitted$testName)
studies <- as.factor(paste(seroFitted$testName, seroFitted$citationID))
assayDataList <- list(N=nrow(seroFitted),
                  K=length(unique(seroFitted$testName)),
                  M=length(unique(studies)),
                  C=nChar,
                  characteristics=charMatrix,
                  assay=as.integer(assayVec),
                  study=as.integer(studies),
                  timeVecEarly=earlyTime,
                  timeVecLate=lateTime,
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
                                           intercept,
                                           interceptSigma,
                                           slopeSigmaEarly,
                                           slopeSigmaLate,
                                           studySigma,
                                           meanSlopeLate,
                                           assayIntercept[loc],
                                           assaySlopeEarly[loc],
                                           assaySlopeLate[loc],
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

write.csv(posteriorTraces,
          "../data/analysis_results/06_late_slope_posterior_samples.csv",
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

write.csv(parameterSummary, "../data/analysis_results/06_parameter_summary.csv",
          row.names=FALSE)

