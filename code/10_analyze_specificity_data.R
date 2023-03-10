##################################
# 
# This script analyzes fits a model to assay specificity,
# using assay characteristics. The outputs are similar to
# scripts 04 and 05. The results are reported on
# Supplementary Section G of the associated paper
# https://www.medrxiv.org/content/10.1101/2022.09.08.22279731v3
# now in press at Eurosurveillance
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
# options of characteristics to fit:
# antigen, antibody, technique, design, fullModel, antigen_technique
characteristics <- "fullModel"

############
# Load data generated in script 03
############
specificityDf <- read.csv("../data/processed_data/assay_specificity.csv",
                       stringsAsFactors=FALSE)

############
# Add binary columns to dataframe, to work as indicator variables
############
specificityDf$RBD <- stringr::str_detect(specificityDf$antigenTarget, "RBD")
specificityDf$S <- stringr::str_detect(specificityDf$antigenTarget, "S") &
  !specificityDf$RBD
specificityDf$N <- stringr::str_detect(specificityDf$antigenTarget, "N") &
  !(specificityDf$RBD | specificityDf$S)

specificityDf$LFA <- stringr::str_detect(specificityDf$technique, "LFIA")
specificityDf$Indirect <- ! specificityDf$LFA &
  stringr::str_detect(specificityDf$design, "indirect")
specificityDf$Direct <- ! specificityDf$LFA &
  stringr::str_detect(specificityDf$design, "sandwich")
specificityDf$Neutralization <- ! specificityDf$LFA &
  stringr::str_detect(specificityDf$design, "competitive")

############
# Make indicator matrix with assay characters to fit
# Uncomment chung of code corresponding to the characteristic of
# interest
############

if (characteristics=="antigen") {
  ### Fit ANTIGEN characteristics
  charsName <- c("N", "S", "RBD")
  fileIdentifier <- "antigen"
} else if (characteristics=="technique") {
  #### Fit ANALYTIC TECHNIQUE characteristics
  charsName <- c("LFA", "Indirect", "Direct", "Neutralization")
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
  charsName <- c("N", "S", "RBD", "LFA", "Direct", "Neutralization")
  fileIdentifier <- "fullModel"
} else if (characteristics=="antigen_technique") {
  #### Fit model with all but design
  charsName <- c("N", "S", "RBD", "LFA")
  fileIdentifier <- "antigen_technique"
}

# Find rows without data and filter out
specificityDf <- dplyr::filter(specificityDf, !is.na(trueNegatives))

# Make characteristics matrix from selected columns
charMatrix <- make_characteristics_matrix(specificityDf, charsName)
naRows <- which(is.na(rowSums(charMatrix)))

# Find which assays have NA values, and filter out
if (length(naRows)!=0) {
  assayVecTemp <- as.factor(specificityDf$testName)
  naAssays <- levels(assayVecTemp)[naRows]
  specificityDf <- dplyr::filter(specificityDf, !(testName %in% naAssays))
  charMatrix <- charMatrix[-naRows,]
}

# make useful vectors
assayVec <- as.factor(specificityDf$testName)
studies <- as.factor(paste(specificityDf$testName, specificityDf$citation))
nChars <- length(charsName)

###################
# Fit the model to the data
###################

# compile model
#regression_model <- rstan::stan_model("./stan_models/specificity_chars.stan",
#                                 model_name="specificity",
#                                 warn_pedantic=TRUE)
# load model if available
regression_model <- readRDS("./stan_models/specificity_chars.RDS")

# Make a list with the input we pass to STAN
assayDataList <- list(N=nrow(specificityDf),
                  K=length(unique(specificityDf$testName)),
                  M=length(unique(studies)),
                  C=nChars,
                  characteristics=charMatrix,
                  assay=as.integer(assayVec),
                  study=as.integer(studies),
                  nNegative=specificityDf$trueNegatives,
                  nTested=specificityDf$nSamplesSpec)

# Fit model
model <- rstan::sampling(regression_model, data=assayDataList,
                           chains=nChains, iter=nIter, refresh=0,
                           warmup=warmup,
                           verbose=TRUE, cores=nCores)

###################
# Extract and save the posterior samples of the fit
###################
posteriorTraces <- tidybayes::gather_draws(model,
                                           interceptSigma,
                                           studySigma,
                                           assayIntercept[loc],
                                           studyIntercept[stu],
                                           charIntercept[par])

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

write.csv(posteriorTraces, paste("../data/analysis_results/10_specificity_",
                fileIdentifier, "_posterior_samples.csv", sep=""), row.names=FALSE)


# Get the sensitivity profiles for the average test of
# each combination of characteristics
uniqueCharRows <- unique(charMatrix)
meanCharProfileDf <- NULL
assayFitDf <- NULL
for (cc in c(1:nrow(uniqueCharRows))) {
  charRow <- as.logical(uniqueCharRows[cc,])
  charCombs <- charsName[charRow]
  charProfileDf <- mean_param_specificity(posteriorTraces,
                           charCombs=list(charCombs)) %>%
    dplyr::select(., -chars)
  rownames(charProfileDf) <- NULL
  for (cn in c(1:length(charsName))) {
    charProfileDf[[charsName[cn]]] <- charRow[cn]
  }
  meanCharProfileDf <- rbind(meanCharProfileDf, charProfileDf)
}
meanCharProfileDf$testName <- NA

write.csv(meanCharProfileDf, paste("../data/analysis_results/10_specificity_",
                fileIdentifier, "_specificity_averages.csv", sep=""), row.names=FALSE)

###################
# Save the predicted sensitivity profile for each test,
# and the test specific parameters
###################

# Get summary statistics of the fitted parameters
parameterSummary <- ungroup(posteriorTraces) %>%
  group_by(., .variable, assay, parName) %>%
  summarize(., paramMean=mean(.value), paramMedian=median(.value),
            paramSD=sd(.value),
            paramL=quantile(.value, probs=0.025),
            paramQuartileL=quantile(.value, probs=0.25),
            paramQuartileH=quantile(.value, probs=0.75),
            paramH=quantile(.value, probs=0.975))

write.csv(parameterSummary, paste("../data/analysis_results/10_specificity_",
            fileIdentifier, "_parameter_summary.csv", sep=""), row.names=FALSE)


