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
# options of characteristics to fit:
# antigen, antibody, technique, design, fullModel
characteristics <- "fullModel"

############
# Load data generated in script 03
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
seroFitted$sandwich <- stringr::str_detect(seroFitted$design, "sandwich")
seroFitted$indirect <- stringr::str_detect(seroFitted$design, "indirect")
seroFitted$competitive <- stringr::str_detect(seroFitted$design, "competitive")
seroFitted$notSandwich <- seroFitted$indirect | seroFitted$competitive


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
} else if (characteristics=="design") {
  charsName <- c("sandwich", "notSandwich", "LFA")
  fileIdentifier <- "design"
} else if (characteristics=="fullModel") {
  #### Fit full model characteristics
  charsName <- c("S", "N", "RBD", "LFA", "sandwich")
  fileIdentifier <- "fullModel"
}

# Make characteristics matrix from selected columns
charMatrix <- make_characteristics_matrix(seroFitted, charsName)
naRows <- which(is.na(rowSums(charMatrix)))

# Find which assays have NA values, and filter out
if (length(naRows)!=0) {
  assayVecTemp <- as.factor(seroFitted$testName)
  naAssays <- levels(assayVecTemp)[naRows]
  seroFitted <- dplyr::filter(seroFitted, !(testName %in% naAssays))
  charMatrix <- charMatrix[-naRows,]
}

# make useful vectors
assayVec <- as.factor(seroFitted$testName)
studies <- as.factor(paste(seroFitted$testName, seroFitted$citationID))
nChars <- length(charsName)

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
paramSize <- c(nChars, 1, 1, 1, 1, nTests, nTests, nStudies) # size of param vectors

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
assayDataList <- list(N=nrow(seroFitted),
                  K=length(unique(seroFitted$testName)),
                  M=length(unique(studies)),
                  C=nChars,
                  characteristics=charMatrix,
                  assay=as.integer(assayVec),
                  study=as.integer(studies),
                  timeVec=seroFitted$testTime,
                  nPositive=seroFitted$nSeropositives,
                  nTested=seroFitted$nSamples)

# Fit model
model <- rstan::sampling(regression_model, data=assayDataList,
                           chains=nChains, iter=nIter, refresh=0,
                           warmup=warmup,
                           verbose=TRUE, cores=nCores, init=initList)

###################
# Extract and save the posterior samples of the fit
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
studyAssayVec <- NULL
studyVec <- NULL
for (i in c(1:length(studyList))) {
  studyAssayVec <- c(studyAssayVec, studyList[[i]][1])
  studyVec <- c(studyVec, studyList[[i]][2])
}
posteriorTraces$studyAssay <- studyAssayVec[posteriorTraces$stu]
posteriorTraces$studyCitation <- studyVec[posteriorTraces$stu]
posteriorTraces <- ungroup(posteriorTraces)

write.csv(posteriorTraces, paste("../data/analysis_results/05_characteristics_",
                fileIdentifier, "_posterior_samples.csv", sep=""), row.names=FALSE)


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
  testRow <- which(seroFitted$testName==tN)[1]
  testCharacteristics <- seroFitted[charsName][testRow,]
  testProfileDf <- merge(testProfileDf, testCharacteristics)
  assayFitDf <- rbind(assayFitDf, testProfileDf)
}

# Get the sensitivity profiles for the average test of
# each combination of characteristics
uniqueCharRows <- unique(charMatrix)
meanCharProfileDf <- NULL
for (cc in c(1:nrow(uniqueCharRows))) {
  charRow <- as.logical(uniqueCharRows[cc,])
  charCombs <- charsName[charRow]
  charProfileDf <- mean_param_seroreversion(posteriorTraces, timeVec,
                           interceptName="intercept",
                           charCombs=list(charCombs),
                           timeNormalization=1) %>%
    dplyr::select(., -chars)
  for (cn in c(1:length(charsName))) {
    charProfileDf[[charsName[cn]]] <- charRow[cn]
  }
  meanCharProfileDf <- rbind(meanCharProfileDf, charProfileDf)
}
meanCharProfileDf$testName <- NA
assayFitDf <- rbind(meanCharProfileDf, assayFitDf)

# Export sensitivity profiles
write.csv(assayFitDf, paste("../data/analysis_results/05_characteristics_",
            fileIdentifier, "_assay_sensitivity_curve.csv", sep=""), row.names=FALSE)

# Get summary statistics of the fitted parameters
parameterSummary <- ungroup(posteriorTraces) %>%
  group_by(., .variable, assay, parName) %>%
  summarize(., paramMean=mean(.value), paramMedian=median(.value),
            paramSD=sd(.value),
            paramL=quantile(.value, probs=0.025),
            paramQuartileL=quantile(.value, probs=0.25),
            paramQuartileH=quantile(.value, probs=0.75),
            paramH=quantile(.value, probs=0.975))

# Get the summary of parameter combinations and add them
# to the parmeter summary
multiCharDf <- NULL
for (cc in c(1:nrow(uniqueCharRows))) {
  charRow <- as.logical(uniqueCharRows[cc,])
  if (sum(charRow)>1) {
    charCombs <- charsName[charRow]
    combSlopes <- get_combination_slopes(posteriorTraces, list(charCombs))
    combName <- names(combSlopes)[2]
    names(combSlopes)[2] <- ".value"
    combSlopes <- dplyr::mutate(combSlopes, .variable="charSlope",
                                assay=NA, parName=combName)
    combSlopesDf <- combSlopes %>%
      group_by(., .variable, assay, parName) %>%
      summarize(., paramMean=mean(.value), paramMedian=median(.value),
            paramSD=sd(.value),
            paramL=quantile(.value, probs=0.025),
            paramQuartileL=quantile(.value, probs=0.25),
            paramQuartileH=quantile(.value, probs=0.75),
            paramH=quantile(.value, probs=0.975))
    parameterSummary <- rbind(parameterSummary, combSlopesDf)
  }
}

write.csv(parameterSummary, paste("../data/analysis_results/05_characteristics_",
            fileIdentifier, "_parameter_summary.csv", sep=""), row.names=FALSE)


