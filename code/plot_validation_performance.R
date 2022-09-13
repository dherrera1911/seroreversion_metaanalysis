####################################################
#
# This script doesn't make plots, but computes the
# performance at validation for the different models
#
####################################################


library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(bayesplot)
library(tidybayes)
library(viridis)
library(RColorBrewer)
library(ggthemes)
library(lemon)
library(cowplot)
library(ggpubr)
source("./functions_auxiliary.R")
source("./functions_seroreversion_fit_analysis.R")

regLineSize=1
ribbonAlpha=0.2
locLineSize=0.4
locLineAlpha=0.4
locPointSize=1.6
sampleLineSize=0.1
sampleLineAlpha=0.1
dataAlpha <- 0.5

# Load the CV of the basic model, done with grouped validation
basicModelGroupedCV <- read.csv("../data/analysis_results/04_predicted_sensitivities_grouped_CV.csv",
                         stringsAsFactors=FALSE) %>%
dplyr::mutate(., inInterval=(nSeropositives>=predictedPosL) & (nSeropositives<=predictedPosH))

# Compute the proportion of data points within the CrI
basicModelGroupedCVPerformance <- mean(basicModelGroupedCV$inInterval)

testSamples <- arrange(testSamples, nSamples)
fewSampleTests <- testSamples$testName[testSamples$nSamples<=9]
nFewSamples <- sum(testSamples$nSamples[testSamples$nSamples<=9])
nManySamples <- sum(testSamples$nSamples[testSamples$nSamples>9])

fewSampleVal_basic <- dplyr::filter(basicModelGroupedCV, testName %in% fewSampleTests)

basicFewSamplePerformance <- mean(fewSampleVal_basic$inInterval)


## Load the CV of the basic model, done with non-grouped validation
#basicModelRandomCV <- read.csv("../data/analysis_results/04_predicted_sensitivities_random_CV.csv",
#                         stringsAsFactors=FALSE) %>%
#dplyr::mutate(., inInterval=(nSeropositives>=predictedPosL) & (nSeropositives<=predictedPosH))
#
## Compute the proportion of data points within the CrI
#basicModelRandomCVPerformance <- mean(basicModelGroupedCV$inInterval)
#
## Check the validation performance in assays that have few datapoints
#testSamples <- group_by(basicModelGroupedCV, testName) %>%
#  summarize(., multiTime=length(unique(testTime)) > 1,
#            nTimes=length(unique(testTime)),
#            meanTime=mean(testTime),
#            nSamples=n())
#


# Load full model, grouped CV results
fullModelGroupedCV <- read.csv("../data/analysis_results/05_characteristics_fullModel_grouped_CV.csv",
                         stringsAsFactors=FALSE) %>%
dplyr::mutate(., inInterval=(nSeropositives>=predictedPosL) & (nSeropositives<=predictedPosH))

# Compute the proportion of data points within the CrI
fullModelGroupedCVPerformance <- mean(fullModelGroupedCV$inInterval)

# Load full model, random CV results
fullModelRandomCV <- read.csv("../data/analysis_results/04_predicted_sensitivities_random_CV.csv",
                         stringsAsFactors=FALSE) %>%
dplyr::mutate(., inInterval=(nSeropositives>=predictedPosL) & (nSeropositives<=predictedPosH))

# Compute the proportion of data points within the CrI
#fullModelRandomCVPerformance <- mean(fullModelGroupedCV$inInterval)

# Full model on tests with few data points
fewSampleVal_full <- dplyr::filter(fullModelGroupedCV, testName %in% fewSampleTests)
fullFewSamplePerformance <- mean(fewSampleVal_full$inInterval)


# Antigen
antigenValGr <- read.csv("../data/analysis_results/05_characteristics_antigen_grouped_CV.csv",
                         stringsAsFactors=FALSE) %>%
dplyr::mutate(., inInterval=(nSeropositives>=predictedPosL) &
              (nSeropositives<=predictedPosH))

antigenValGrPerformance <- mean(antigenValGr$inInterval)

antigenValRn <- read.csv("../data/analysis_results/04_predicted_sensitivities_random_CV.csv",
                         stringsAsFactors=FALSE) %>%
dplyr::mutate(., inInterval=(nSeropositives>=predictedPosL) &
              (nSeropositives<=predictedPosH))

antigenValRnPerformance <- mean(antigenValGr$inInterval)

# Full model on tests with few data points
fewSampleVal_antigen <- dplyr::filter(antigenValGr, testName %in% fewSampleTests)
antigenFewSamplePerformance <- mean(fewSampleVal_antigen$inInterval)


# Technique
#techniqueValGr <- read.csv("../data/analysis_results/05_characteristics_technique_grouped_CV.csv",
#                         stringsAsFactors=FALSE) %>%
#dplyr::mutate(., inInterval=(nSeropositives>=predictedPosL) &
#              (nSeropositives<=predictedPosH))
#
#techniqueValGrPerformance <- mean(techniqueValGr$inInterval)
#
#techniqueValRn <- read.csv("../data/analysis_results/04_predicted_sensitivities_random_CV.csv",
#                         stringsAsFactors=FALSE) %>%
#dplyr::mutate(., inInterval=(nSeropositives>=predictedPosL) &
#              (nSeropositives<=predictedPosH))
#
#techniqueValRnPerformance <- mean(techniqueValGr$inInterval)
#
## Full model on tests with few data points
#fewSampleVal_technique <- dplyr::filter(techniqueValGr, testName %in% fewSampleTests)
#techniqueFewSamplePerformance <- mean(fewSampleVal_technique$inInterval)
#


# Compare the sizes of the CrI for the basic vs full model
basic2 <- filter(basicModelGroupedCV, testName %in% fullModelGroupedCV$testName) %>%
  dplyr::mutate(., intervalWidthBasic=(sensitivityHPred-sensitivityLPred)) %>%
  dplyr::select(., testName, citationID, testTime, intervalWidthBasic)

full2 <- fullModelGroupedCV %>%
  dplyr::mutate(., intervalWidthFull=(sensitivityHPred-sensitivityLPred)) %>%
  select(., testName, citationID, testTime, intervalWidthFull)

comparisonDf <- merge(basic2, full2)

narrowerFull <- with(comparisonDf, intervalWidthFull<intervalWidthBasic)
mean(narrowerFull)

