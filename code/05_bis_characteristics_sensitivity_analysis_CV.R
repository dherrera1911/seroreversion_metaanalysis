##################################
##################################
# 
# This script fits a model to data on time-varying
# sensitivity of different serology assays. Only data
# of serology testing on previously diagnosed individuals
# is used.
# 
# The fitting is done taking into account different test
# characteristics
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
warmup <- 1000

characteristics <- "antigen"

crossValidationType <- "grouped"
testsPerGroup <- 2 # for grouped CV only
#
#crossValidationType <- "random"
#pointsPerGroup <- 7


############
# Load data generated in script 06
############
seroFitted <- read.csv("../data/processed_data/PCR_to_serotest_all.csv",
                       stringsAsFactors=FALSE)

############
# Add binary columns to dataframe, to work as indicator variables
############
seroFitted$S <- stringr::str_detect(seroFitted$antigenTarget, "S") |
  stringr::str_detect(seroFitted$antigenTarget, "RBD")
seroFitted$RBD <- stringr::str_detect(seroFitted$antigenTarget, "RBD")
seroFitted$N <- stringr::str_detect(seroFitted$antigenTarget, "N")
seroFitted$SN <- seroFitted$S & seroFitted$N
seroFitted$LFA <- stringr::str_detect(seroFitted$technique, "LFIA")
seroFitted$Rest <- !seroFitted$LFA
seroFitted$ELISA <- stringr::str_detect(seroFitted$technique, "ELISA")
seroFitted$chemlum <- (stringr::str_detect(seroFitted$technique, "CMIA") |
                       stringr::str_detect(seroFitted$technique, "CLIA"))
seroFitted$IgM <- stringr::str_detect(seroFitted$antibodyTarget, "IgM")
seroFitted$IgA <- stringr::str_detect(seroFitted$antibodyTarget, "IgA")
seroFitted$IgG <- TRUE
seroFitted$Total <- seroFitted$IgM & seroFitted$IgA


############
# Make indicator matrix with assay characters to fit
# Uncomment chung of code corresponding to the characteristic of
# interest
############

if (characteristics=="antigen") {
  ### Fit ANTIGEN characteristics
  charsName <- c("S", "N", "RBD", "SN")
  fileIdentifier <- "antigen"
} else if (characteristics=="technique") {
  #### Fit ANALYTIC TECHNIQUE characteristics
  charsName <- c("LFA", "Rest")
  fileIdentifier <- "technique"
} else if (characteristics=="antibody") {
  #### Fit ANTIBODY SEROTYPE characteristics
  charsName <- c("IgG", "IgM", "IgA", "Total")
  fileIdentifier <- "antibody"
} else if (characteristics=="fullModel") {
  #### Fit full model characteristics
  charsName <- c("S", "N", "RBD", "LFA")
  fileIdentifier <- "fullModel"
}

# Make characteristics matrix from selected columns
charMatrix <- make_characteristics_matrix(seroFitted, charsName)
naRows <- which(is.na(rowSums(charMatrix)))

# Find which assays have NA values, and filter out
assayVecTemp <- as.factor(seroFitted$testName)
naAssays <- levels(assayVecTemp)[naRows]
seroFitted <- dplyr::filter(seroFitted, !(testName %in% naAssays))
charMatrix <- charMatrix[-naRows,]

# make useful vectors
assayVec <- as.factor(seroFitted$testName)
studies <- as.factor(paste(seroFitted$testName, seroFitted$citationID))
nChars <- length(charsName)

###############
# compile model
###############
#outcome_reg <- rstan::stan_model("./sensitivity_change_chars.stan",
#                                 model_name="time_change_sensitivity",
#                                 warn_pedantic=TRUE)
# load model if available
outcome_reg <- readRDS("./sensitivity_change_chars.RDS")


###############
### Make the cross validation batch groups
###############
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


###############
### Set some initialization parameters for the fits
###############
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

###########################
### Fit and validate
###########################
nValGroups <- length(unique(validationGroup))
allValidationDf <- NULL

for (vg in c(1:nValGroups)) {
  fitDf <- dplyr::filter(seroFitted, validationGroup!=vg)
  valDf <- dplyr::filter(seroFitted, validationGroup==vg)

  #############
  # Make initialization values
  #############
  nTests <- length(unique(fitDf$testName))
  nStudies <- length(unique(paste(fitDf$testName, fitDf$citationID)))
  paramSize <- c(1, 1, 1, 1, 1, nTests, nTests, nStudies) # size of param vectors
  initList <- sample_initial_values(nChains=nChains, paramListName=paramListName,
                             lowerBound=lowerBound, upperBound=upperBound,
                             paramSize=paramSize)

  #############
  # Make data list to input to STAN
  #############
  charMatrix <- make_characteristics_matrix(fitDf, charsName)
  assayVec <- as.factor(fitDf$testName)
  studies <- as.factor(paste(fitDf$testName, fitDf$citationID))
  assayDataList <- list(N=nrow(fitDf),
                    K=length(unique(fitDf$testName)),
                    M=length(unique(studies)),
                    C=nChars,
                    characteristics=charMatrix,
                    assay=as.integer(assayVec),
                    study=as.integer(studies),
                    timeVec=fitDf$testTime,
                    nPositive=fitDf$nSeropositives,
                    nTested=fitDf$nSamples)

  #############
  # Fit model
  #############
  model <- rstan::sampling(outcome_reg, data=assayDataList,
                             chains=nChains, iter=nIter,
                             warmup=warmup, refresh=0,
                             verbose=TRUE, cores=nCores, init=initList)

  #############
  # Extract and tidy posterior
  #############
  posteriorTraces <- tidybayes::gather_draws(model, 
                                             intercept, slopeSigma,
                                             interceptSigma,
                                             studySigma,
                                             assayIntercept[loc],
                                             assaySlope[loc],
                                             studyIntercept[stu],
                                             charSlope[par])

  posteriorTraces$assay <- levels(assayVec)[posteriorTraces$loc]
  posteriorTraces$parName <- charsName[posteriorTraces$par]
  studyList <- strsplit(as.character(levels(studies)), " n-")
  studyAssayVec <- NULL
  studyVec <- NULL
  for (i in c(1:length(studyList))) {
    studyAssayVec <- c(studyAssayVec, studyList[[i]][1])
    studyVec <- c(studyVec, studyList[[i]][2])
  }
  posteriorTraces$studyAssay <- studyAssayVec[posteriorTraces$stu]
  posteriorTraces$studyCitation <- studyVec[posteriorTraces$stu]

  #############
  # Predict the outcome of the left out samples, and append to validation Df
  #############
  valOutput <- sensitivity_prediction_validation(posteriorTraces, valDf,
    binomSamples=5)
  allValidationDf <- rbind(valOutput, allValidationDf)
}

valName <- paste("../data/analysis_results/05_characteristics_",
                 fileIdentifier, "_", crossValidationType, "_CV.csv", sep="")
write.csv(allValidationDf, valName, row.names=FALSE)


