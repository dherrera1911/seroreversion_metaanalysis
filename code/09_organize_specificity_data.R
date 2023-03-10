##################################
# 
# This script organizes and puts together the specificity
# data reported by different sources.
#
# It generates the file "../data/processed_data/assay_specificity.csv"
# which is then used to make the specificity analysis.
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
source("./functions_auxiliary.R")

############
# Load data generated in script 03
############
seroFitted <- read.csv("../data/processed_data/PCR_to_serotest_all.csv",
                       stringsAsFactors=FALSE)
testNames <- unique(seroFitted$testName)

# Load the manufacturer reported data
manufacturerDf <- read.csv("../data/raw_data/assay_performance_manufacturer.csv",
                       stringsAsFactors=FALSE, na.strings=c("", "--")) %>%
  dplyr::rename(., testName=test_name, sensitivityMean=sensitivity_estimate,
                sensitivityL=sen_lower95CI, sensitivityH=sen_upper95CI,
                specificityMean=specificity_estimate, specificityL=spec_lower95CI,
                specificityH=spec_upper95CI)

# Tidy specificity numbers of manufacturers
manufacturerDf$specificityMean[manufacturerDf$specificityMean=="see note to right"] <- NA
manufacturerDf$specificityMean <- as.numeric(manufacturerDf$specificityMean)
manufacturerDf$specificityL <- as.numeric(manufacturerDf$specificityL)
manufacturerDf$specificityH <- as.numeric(manufacturerDf$specificityH)

# Tidy the names of the data frame with manufacturer data
originalNames <- names(manufacturerDf)
names(manufacturerDf)[stringr::str_detect(originalNames, "source")] <- "baselineSource"
names(manufacturerDf)[stringr::str_detect(originalNames, "testing.seropositive")] <- "nSeropositives"
names(manufacturerDf)[stringr::str_detect(originalNames, "total.number.*sensitivity")] <- "nSamples"
names(manufacturerDf)[stringr::str_detect(originalNames, "total.number.*specificity")] <- "nSamplesSpec"
names(manufacturerDf)[stringr::str_detect(originalNames, "number.of.COVID.*seronegative")] <- "trueNegatives"
manufacturerDf <- dplyr::filter(manufacturerDf, !is.na(baselineSource))

for (r in c(1:nrow(manufacturerDf))) {
  # Calculate sensitivity for rows with raw numbers but no estimate
  if (with(manufacturerDf, is.na(sensitivityL[r]) & !is.na(nSamples[r]))) {
    seroprevConfint <- binomial_confint(manufacturerDf[[r,"nSamples"]],
                                        manufacturerDf[[r, "nSeropositives"]])
    manufacturerDf[[r,"sensitivityL"]] <- signif(seroprevConfint$lower*100, digits=4)
    manufacturerDf[[r,"sensitivityH"]] <- signif(seroprevConfint$upper*100, digits=4)
    manufacturerDf[[r,"sensitivityMean"]] <- signif(manufacturerDf[[r,"nSeropositives"]]/
                                              manufacturerDf[[r,"nSamples"]]*100, digits=4)
  }
  # Calculate specificity for rows with raw numbers but no estimate
  if (with(manufacturerDf, is.na(specificityL[r]) & !is.na(trueNegatives[r]))) {
    specConfint <- binomial_confint(manufacturerDf[[r,"nSamplesSpec"]],
                                        manufacturerDf[[r, "trueNegatives"]])
    manufacturerDf[[r,"specificityL"]] <- signif(specConfint$lower*100, digits=4)
    manufacturerDf[[r,"specificityH"]] <- signif(specConfint$upper*100, digits=4)
    manufacturerDf[[r,"specificityMean"]] <- signif(manufacturerDf[[r,"trueNegatives"]]/
                                              manufacturerDf[[r,"nSamplesSpec"]]*100, digits=4)
  }
}

specificityDf <- dplyr::filter(manufacturerDf, testName %in% testNames)

## Add assay characteristics to dataframe
assayChars <- read.csv("../data/raw_data/assay_characteristics.csv",
                       stringsAsFactors=FALSE, na.strings=c("", "--")) %>%
  dplyr::rename(., testName=test_name_long, antigenTarget=antigen.target,
  isotypeTarget=isotype.target, technique=assay.type, 
  design=test.design) %>%
  dplyr::filter(., testName %in% testNames) %>%
  dplyr::select(., testName, test_name_simplified, antigenTarget,
                design, technique)

specificityDf <- merge(specificityDf, assayChars)

# Fit beta distribution to studies without raw data, to estimate raw data
nonRaw <- which(is.na(specificityDf$nSamplesSpec) &
                !is.na(specificityDf$specificityL))

fittedBetas <- fit_beta_ci(meanEstimate=specificityDf$sensitivityMean[nonRaw]/100,
                            lower=specificityDf$sensitivityL[nonRaw]/100,
                            upper=specificityDf$sensitivityH[nonRaw]/100)

totalCount <- with(fittedBetas, shape1+shape2)
meanSensitivity <- with(fittedBetas, shape1/(shape1+shape2)*100)
confintSensitivity <- binomial_confint(round(totalCount), round(fittedBetas$shape1))
# Add estimated raw data to the dataframe
specificityDf$nSamplesSpec[nonRaw] <- with(fittedBetas, round(shape1+shape2))
specificityDf$trueNegatives[nonRaw] <- round(fittedBetas$shape1)

write.csv(specificityDf, "../data/processed_data/assay_specificity.csv",
          row.names=FALSE)


