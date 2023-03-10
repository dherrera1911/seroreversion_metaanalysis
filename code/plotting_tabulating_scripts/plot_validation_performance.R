####################################################
#
# This script doesn't make plots, but computes the
# performance at cross validation for the different models
# that are reported in the paper, as well as the
# comparison of CI width between the model that
# includes assay characteristics and the model that doesn't
# 
# Script authored by Daniel Herrera-Esposito.
# For questions, contact me at dherrera1911[at]gmail.com
# 
# Final version revised 10/03/2023
# 
####################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
source("../functions_auxiliary.R")
source("../functions_seroreversion_fit_analysis.R")

# Load the CV of the basic model, done with grouped validation
basicModelCV <- read.csv("../../data/analysis_results/04_predicted_sensitivities_grouped_CV.csv",
                         stringsAsFactors=FALSE) %>%
dplyr::mutate(., inInterval=(nSeropositives>=predictedPosL) & (nSeropositives<=predictedPosH))

# Compute the proportion of data points within the CrI
basicModelCVPerformance <- mean(basicModelCV$inInterval)

testSamples <- group_by(basicModelCV, testName) %>%
  summarize(., multiTime=length(unique(testTime)) > 1,
            nTimes=length(unique(testTime)),
            meanTime=mean(testTime),
            nSamples=n())
testSamples <- arrange(testSamples, nSamples)
fewSampleTests <- testSamples$testName[testSamples$nSamples<=9]
nFewSamples <- sum(testSamples$nSamples[testSamples$nSamples<=9])
nManySamples <- sum(testSamples$nSamples[testSamples$nSamples>9])

fewSampleVal_basic <- dplyr::filter(basicModelCV, testName %in% fewSampleTests)

basicFewSamplePerformance <- mean(fewSampleVal_basic$inInterval)

# Load full model, grouped CV results
fullModelCV <- read.csv("../../data/analysis_results/05_characteristics_fullModel_grouped_CV.csv",
                         stringsAsFactors=FALSE) %>%
dplyr::mutate(., inInterval=(nSeropositives>=predictedPosL) & (nSeropositives<=predictedPosH))

# Compute the proportion of data points within the CrI
fullModelCVPerformance <- mean(fullModelCV$inInterval)

# Full model on tests with few data points
fewSampleVal_full <- dplyr::filter(fullModelCV, testName %in% fewSampleTests)
fullFewSamplePerformance <- mean(fewSampleVal_full$inInterval)


# Compare the sizes of the CrI for the basic vs full model
basic2 <- filter(basicModelCV, testName %in% fullModelCV$testName) %>%
  dplyr::mutate(., intervalWidthBasic=(sensitivityHPred-sensitivityLPred)) %>%
  dplyr::select(., testName, citationID, testTime, intervalWidthBasic)

full2 <- fullModelCV %>%
  dplyr::mutate(., intervalWidthFull=(sensitivityHPred-sensitivityLPred)) %>%
  select(., testName, citationID, testTime, intervalWidthFull)

comparisonDf <- merge(basic2, full2)

narrowerFull <- with(comparisonDf, intervalWidthFull<intervalWidthBasic)
mean(narrowerFull)

