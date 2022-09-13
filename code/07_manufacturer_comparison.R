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
# Load data generated in script 03
############
seroFitted <- read.csv("../data/processed_data/PCR_to_serotest_all.csv",
                       stringsAsFactors=FALSE)
testNames <- unique(seroFitted$testName)


### Load sensitivity profiles
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
                sensitivityL=sen_lower95CI, sensitivityH=sen_upper95CI)

originalNames <- names(reportedSens)
names(reportedSens)[stringr::str_detect(originalNames, "source")] <- "baselineSource"
names(reportedSens)[stringr::str_detect(originalNames, "testing.seropositive")] <- "nSeropositives"
names(reportedSens)[stringr::str_detect(originalNames, "total.number.*sensitivity")] <- "nSamples"

reportedSens <- dplyr::filter(reportedSens, !is.na(baselineSource))

for (r in c(1:nrow(reportedSens))) {
  if (with(reportedSens, is.na(sensitivityL[r]) & !is.na(nSamples[r]))) {
    seroprevConfint <- binomial_confint(reportedSens[[r,"nSamples"]],
                                        reportedSens[[r, "nSeropositives"]])
    reportedSens[[r,"sensitivityL"]] <- signif(seroprevConfint$lower*100, digits=4)
    reportedSens[[r,"sensitivityH"]] <- signif(seroprevConfint$upper*100, digits=4)
    reportedSens[[r,"sensitivityMean"]] <- signif(reportedSens[[r,"nSeropositives"]]/
                                              reportedSens[[r,"nSamples"]]*100, digits=4)
  }
}


# compare baseline to seroreversion
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
    testBaselineMean <- reportedSens$sensitivityMean[manufacturerRow[1]]
    testBaselineL <- reportedSens$sensitivityL[manufacturerRow[1]]
    testBaselineH <- reportedSens$sensitivityH[manufacturerRow[1]]
    testSource <- "M"
  } else if (length(fdaRow>0)) {
    testBaselineMean <- reportedSens$sensitivityMean[fdaRow[1]]
    testBaselineL <- reportedSens$sensitivityL[fdaRow[1]]
    testBaselineH <- reportedSens$sensitivityH[fdaRow[1]]
    testSource <- "F"
  } else if (length(otherRow>0)) {
    testBaselineMean <- reportedSens$sensitivityMean[otherRow[1]]
    testBaselineL <- reportedSens$sensitivityL[otherRow[1]]
    testBaselineH <- reportedSens$sensitivityH[otherRow[1]]
    testSource <- "O"
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
                      time75=cutoff1, time50=cutoff2, month2Frac=month2Frac)
  baselineDf <- rbind(baselineDf, dfRow)
}

write.csv(baselineDf, "../data/analysis_results/07_manufacturer_comparison.csv",
          row.names=FALSE)


