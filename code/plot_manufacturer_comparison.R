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


fullModelProfile <- read.csv("../data/analysis_results/05_characteristics_fullModel_assay_sensitivity_curve.csv",
                           stringsAsFactors=FALSE)


testNames <- unique(fullModelProfile$testName)
testNames <- testNames[!is.na(testNames)]

assayChars <- read.csv("../data/raw_data/assay_characteristics.csv",
                       stringsAsFactors=FALSE, na.strings=c("", "--")) %>%
  dplyr::rename(., testName=test_name_long, antigen_target=antigen.target,
  assay_type=assay.type, isotype_target=isotype.target) %>%
  dplyr::filter(., testName %in% testNames)

decayTime70 <- NA

sensitivityFraction <- 0.5
for (tN in c(1:length(testNames))) {
  rowNumber<- which(assayChars$testName == testNames[tN])
  testRow <- assayChars[rowNumber,]
  rowName <- testRow$testName
  rowSensitivityReported <- testRow$reported_sensitivity
  testProfile <- dplyr::filter(fullModelProfile, testName==rowName) %>%
    dplyr::mutate(., sensitivityMean=sensitivityMean*100) %>%
    dplyr::filter(., sensitivityMean >= rowSensitivityReported*sensitivityFraction)
  if (!nrow(testProfile)==0) {
    assayChars$decayTime70[rowNumber] <- max(testProfile$time)
  }
  de
}

decayTimes <- sort(assayChars$decayTime70[!is.na(assayChars$decayTime70)])

frac6months <- mean(decayTimes<6)




