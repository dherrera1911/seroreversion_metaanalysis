##################################
# 
# This script compares the results of time-varying assay sensitivity
# to the reference sensitivities reported by manufacturers (or other sources).
# 
# The script generates as output the file
# "../data/analysis_results/07_manufacturer_comparison.csv"
# that shows for every assay, how many months it takes for
# sensitivity to drop to 75% and 50% of the reference
# sensitivity. The source of the reference sensitivity is
# indicated in baselineSource (M=manufacturer, F=FDA, O=other).
#
# The results from this script are reported in the Results section 3
# in the associated manuscript
# https://www.medrxiv.org/content/10.1101/2022.09.08.22279731v3
# now in press in Eurosurveillance
#
# 
# Script authored by Daniel Herrera-Esposito.
# For questions, contact me at dherrera1911[at]gmail.com
# 
# Final version revised 10/03/2023
# 
##################################


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

############
# Get the names of all the tests included in the analysis
############
seroFitted <- read.csv("../data/processed_data/PCR_to_serotest_all.csv",
                       stringsAsFactors=FALSE)
testNames <- unique(seroFitted$testName)

### Load sensitivity profiles obtained in model fitting
fullModelProfile <- read.csv("../data/analysis_results/05_characteristics_fullModel_assay_sensitivity_curve.csv",
                           stringsAsFactors=FALSE)
techniqueProfile <- read.csv("../data/analysis_results/05_characteristics_technique_assay_sensitivity_curve.csv",
                           stringsAsFactors=FALSE)
antigenProfile <- read.csv("../data/analysis_results/05_characteristics_technique_assay_sensitivity_curve.csv",
                           stringsAsFactors=FALSE)

# Load the manufacturer reported data
reportedSens <- read.csv("../data/raw_data/assay_performance_manufacturer.csv",
                       stringsAsFactors=FALSE, na.strings=c("", "--")) %>%
  dplyr::rename(., testName=test_name, sensitivityMean=sensitivity_estimate,
                sensitivityL=sen_lower95CI, sensitivityH=sen_upper95CI,
                specificityMean=specificity_estimate, specificityL=spec_lower95CI,
                specificityH=spec_upper95CI)

# Tidy specificity numbers of manufacturers
reportedSens$specificityMean[reportedSens$specificityMean=="see note to right"] <- NA
reportedSens$specificityMean <- as.numeric(reportedSens$specificityMean)
reportedSens$specificityL <- as.numeric(reportedSens$specificityL)
reportedSens$specificityH <- as.numeric(reportedSens$specificityH)

# Tidy the names of the data frame with manufacturer data
originalNames <- names(reportedSens)
names(reportedSens)[stringr::str_detect(originalNames, "source")] <- "baselineSource"
names(reportedSens)[stringr::str_detect(originalNames, "testing.seropositive")] <- "nSeropositives"
names(reportedSens)[stringr::str_detect(originalNames, "total.number.*sensitivity")] <- "nSamples"
names(reportedSens)[stringr::str_detect(originalNames, "total.number.*specificity")] <- "nSamplesSpec"
names(reportedSens)[stringr::str_detect(originalNames, "number.of.COVID.*seronegative")] <- "trueNegatives"
reportedSens <- dplyr::filter(reportedSens, !is.na(baselineSource))

for (r in c(1:nrow(reportedSens))) {
  # Calculate sensitivity for rows with raw numbers but no estimate
  if (with(reportedSens, is.na(sensitivityL[r]) & !is.na(nSamples[r]))) {
    seroprevConfint <- binomial_confint(reportedSens[[r,"nSamples"]],
                                        reportedSens[[r, "nSeropositives"]])
    reportedSens[[r,"sensitivityL"]] <- signif(seroprevConfint$lower*100, digits=4)
    reportedSens[[r,"sensitivityH"]] <- signif(seroprevConfint$upper*100, digits=4)
    reportedSens[[r,"sensitivityMean"]] <- signif(reportedSens[[r,"nSeropositives"]]/
                                              reportedSens[[r,"nSamples"]]*100, digits=4)
  }
  # Calculate specificity for rows with raw numbers but no estimate
  if (with(reportedSens, is.na(specificityL[r]) & !is.na(trueNegatives[r]))) {
    specConfint <- binomial_confint(reportedSens[[r,"nSamplesSpec"]],
                                        reportedSens[[r, "trueNegatives"]])
    reportedSens[[r,"specificityL"]] <- signif(specConfint$lower*100, digits=4)
    reportedSens[[r,"specificityH"]] <- signif(specConfint$upper*100, digits=4)
    reportedSens[[r,"specificityMean"]] <- signif(reportedSens[[r,"trueNegatives"]]/
                                              reportedSens[[r,"nSamplesSpec"]]*100, digits=4)
  }
}


### COMPARE BASELINE TO ESTIMATES OBTAINED FROM SEROREVERSION FITTING
sensitivityFractions <- c(0.75, 0.5)
baselineDf <- NULL

for (t in c(1:length(testNames))) {
  tN <- testNames[t]
  # Get the baseline for this test
  testRows <- which(reportedSens$testName==tN)
  manufacturerRow <- with(reportedSens, which(testName==tN & baselineSource=="M"))
  fdaRow <- with(reportedSens, which(testName==tN & baselineSource=="F"))
  otherRow <- with(reportedSens, which(testName==tN & baselineSource=="O"))
  if (length(manufacturerRow>0)) {
    # Extract sensitivity
    testBaselineMean <- reportedSens$sensitivityMean[manufacturerRow[1]]
    testBaselineL <- reportedSens$sensitivityL[manufacturerRow[1]]
    testBaselineH <- reportedSens$sensitivityH[manufacturerRow[1]]
    testSource <- "M"
    # Extract specificity
    testSpecificityMean <- reportedSens$specificityMean[manufacturerRow[1]]
    testSpecificityL <- reportedSens$specificityL[manufacturerRow[1]]
    testSpecificityH <- reportedSens$specificityH[manufacturerRow[1]]
  } else if (length(fdaRow>0)) {
    # Extract sensitivity
    testBaselineMean <- reportedSens$sensitivityMean[fdaRow[1]]
    testBaselineL <- reportedSens$sensitivityL[fdaRow[1]]
    testBaselineH <- reportedSens$sensitivityH[fdaRow[1]]
    testSource <- "F"
    # Extract specificity
    testSpecificityMean <- reportedSens$specificityMean[fdaRow[1]]
    testSpecificityL <- reportedSens$specificityL[fdaRow[1]]
    testSpecificityH <- reportedSens$specificityH[fdaRow[1]]
  } else if (length(otherRow>0)) {
    # Extract sensitivity
    testBaselineMean <- reportedSens$sensitivityMean[otherRow[1]]
    testBaselineL <- reportedSens$sensitivityL[otherRow[1]]
    testBaselineH <- reportedSens$sensitivityH[otherRow[1]]
    testSource <- "O"
    # Extract specificity
    testSpecificityMean <- reportedSens$specificityMean[otherRow[1]]
    testSpecificityL <- reportedSens$specificityL[otherRow[1]]
    testSpecificityH <- reportedSens$specificityH[otherRow[1]]
  } else {
    testBaselineMean <- NA
    testBaselineL <- NA
    testBaselineH <- NA
    testSource <- NA
  }
  # Get the sensitivity profile for this test
  if (tN %in% fullModelProfile$testName) {
    testProfile <- dplyr::filter(fullModelProfile, testName==tN) %>%
      dplyr::mutate(., sensitivityMean=sensitivityMean*100)
  } else if (tN %in% antigenProfile) {
    testProfile <- dplyr::filter(antigenProfile, testName==tN) %>%
      dplyr::mutate(., sensitivityMean=sensitivityMean*100)
  } else if (tN %in% techniqueProfile) {
    testProfile <- dplyr::filter(techniqueProfile, testName==tN) %>%
      dplyr::mutate(., sensitivityMean=sensitivityMean*100)
  }
  # Get the cutoff points
  if (!is.na(testBaselineMean)) {
    temp <- dplyr::filter(testProfile, sensitivityMean >=
                             testBaselineMean*sensitivityFractions[1])
    cutoff1 <- max(temp$time)
    temp <- dplyr::filter(testProfile, sensitivityMean >=
                             testBaselineMean*sensitivityFractions[2])
    cutoff2 <- max(temp$time)
    month2Frac <- with(testProfile, sensitivityMean[time==2]/testBaselineMean)
  } else {
    cutoff1 <- NA
    cutoff2 <- NA
  }
  dfRow <- data.frame(testName=tN, baselineMean=testBaselineMean,
                      baselineSource=testSource,
                      baselineL=testBaselineL, baselineH=testBaselineH,
                      time75=cutoff1, time50=cutoff2, month2Frac=month2Frac,
                      specificityMean=testSpecificityMean,
                      specificityL=testSpecificityL, specificityH=testSpecificityH)
  baselineDf <- rbind(baselineDf, dfRow)
}

## Add assay characteristics to dataframe
assayChars <- read.csv("../data/raw_data/assay_characteristics.csv",
                       stringsAsFactors=FALSE, na.strings=c("", "--")) %>%
  dplyr::rename(., testName=test_name_long, antigenTarget=antigen.target,
  isotypeTarget=isotype.target, technique=assay.type, 
  design=test.design) %>%
  dplyr::filter(., testName %in% testNames) %>%
  dplyr::select(., testName, test_name_simplified, antigenTarget,
                design, technique)

baselineDf <- merge(baselineDf, assayChars)

write.csv(baselineDf, "../data/analysis_results/07_manufacturer_comparison.csv",
          row.names=FALSE)
