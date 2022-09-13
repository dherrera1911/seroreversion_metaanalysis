####################################################
#
# Plot sensitivity across time, both for the 'averages'
# of the different types of assays, and the assay-specific
# sensitivity profiles.
# Generates plots for figures 2B, 3B, 4B, 5B, 6B
# and supplementary Figure S1
# 
# The inputs for plotting and statistical analysis are
# the posterior means and intervals of sensitivity
# across time, that are computed when fitting the
# models, the original data points (to overlay with the fits),
# and a dataset with original data and the results of
# cross-validation.
#
# Files with sensitivity profiles are in ../data/processed_data/
# and have names 'XXX_sensitivity_curve.csv', where XXX
# indicates the model fitted.
# 
# The original data points are in file
# ../data/processed_data/PCR_to_serotest_all.csv
#
# The results of crossvalidation are in
# '../data/analysis_results/XXX_grouped_CV', where
# XXX indicates the model that was used for crossvalidation
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
library(gridExtra)
library(ggthemes)
library(lemon)
library(cowplot)
library(ggpubr)
source("./functions_auxiliary.R")
source("./functions_seroreversion_fit_analysis.R")

#################################
#################################
# Tabulate the mean sensitivity profiles of assay types
#################################
#################################

# Load sensitivity profile of basic model
basicModelProfile <- read.csv("../data/analysis_results/04_assay_sensitivity_curve.csv",
                               stringsAsFactors=FALSE)
fullModelProfile <- read.csv("../data/analysis_results/05_characteristics_fullModel_assay_sensitivity_curve.csv",
                           stringsAsFactors=FALSE)

tabulatedTimes <- seq(1, 14, 1)

meanProfiles <- filter(fullModelProfile, is.na(testName) & time %in% tabulatedTimes)
meanProfiles$sensitivityMean <- as.character(round(meanProfiles$sensitivityMean*100))
meanProfiles$sensitivityL <- as.character(round(meanProfiles$sensitivityL*100))
meanProfiles$sensitivityQL <- as.character(round(meanProfiles$sensitivityQL*100))
meanProfiles$sensitivityH <- as.character(round(meanProfiles$sensitivityH*100))
meanProfiles$sensitivityQH <- as.character(round(meanProfiles$sensitivityQH*100))


meanProfiles$sensitivityStr <- with(meanProfiles,
          paste(sensitivityMean, ' (', sensitivityL, '-', sensitivityH, ')', sep=""))

# additional column to organize plots
meanProfiles$antigen[with(meanProfiles, N)] <- "N"
meanProfiles$antigen[with(meanProfiles, S)] <- "S"
meanProfiles$antigen[with(meanProfiles, RBD)] <- "RBD"


tidySens <- select(meanProfiles, time, antigen, LFA, sandwich, sensitivityStr) %>%
  pivot_wider(., names_from='time', id_cols=c(antigen, LFA, sandwich),
                    values_from=sensitivityStr, names_prefix="")

write.csv(tidySens, "../data/analysis_results/sensitivity_profile_table.csv",
row.names=FALSE)



#################################
#################################
# Tabulate the assay-specific sensitivity profiles
#################################
#################################

assayChars <- read.csv("../data/raw_data/assay_characteristics.csv", stringsAsFactors=FALSE)

assayProfiles <- filter(fullModelProfile, !is.na(testName) & time %in% tabulatedTimes)
assayProfiles$sensitivityMean <- as.character(round(assayProfiles$sensitivityMean*100))
assayProfiles$sensitivityL <- as.character(round(assayProfiles$sensitivityL*100))
assayProfiles$sensitivityQL <- as.character(round(assayProfiles$sensitivityQL*100))
assayProfiles$sensitivityH <- as.character(round(assayProfiles$sensitivityH*100))
assayProfiles$sensitivityQH <- as.character(round(assayProfiles$sensitivityQH*100))


assayProfiles$sensitivityStr <- with(assayProfiles,
          paste(sensitivityMean, ' (', sensitivityL, '-', sensitivityH, ')', sep=""))

assayTidySens <- select(assayProfiles, time, testName, sensitivityStr) %>%
  pivot_wider(., names_from='time', id_cols=testName,
                    values_from=sensitivityStr, names_prefix="")

write.csv(assayTidySens, "../data/analysis_results/sensitivity_profile_table_AssaysFull.csv",
row.names=FALSE)


# Get assays missing from the full model
missingAssaysInd <- which(!basicModelProfile$testName %in% assayProfiles$testName)
missingAssays <- unique(basicModelProfile$testName[missingAssaysInd])

missingAssayChars <- dplyr::filter(assayChars, test_name_long %in% missingAssays) %>%
  dplyr::select(., test_name_long, antigen.target, test.design, assay.type)

# Load technique model profile and save those assays separately
designModelProfile <- read.csv("../data/analysis_results/05_characteristics_design_assay_sensitivity_curve.csv",
                               stringsAsFactors=FALSE)

designModelAssays <- filter(designModelProfile, testName %in% missingAssays & !is.na(testName))

designModelAssays <- filter(designModelAssays, !is.na(testName) & time %in% tabulatedTimes)
designModelAssays$sensitivityMean <- as.character(round(designModelAssays$sensitivityMean*100))
designModelAssays$sensitivityL <- as.character(round(designModelAssays$sensitivityL*100))
designModelAssays$sensitivityQL <- as.character(round(designModelAssays$sensitivityQL*100))
designModelAssays$sensitivityH <- as.character(round(designModelAssays$sensitivityH*100))
designModelAssays$sensitivityQH <- as.character(round(designModelAssays$sensitivityQH*100))

designModelAssays$sensitivityStr <- with(designModelAssays,
          paste(sensitivityMean, ' (', sensitivityL, '-', sensitivityH, ')', sep=""))

assayTidySens_design <- select(designModelAssays, time, testName, sensitivityStr) %>%
  pivot_wider(., names_from='time', id_cols=testName,
                    values_from=sensitivityStr, names_prefix="")

write.csv(assayTidySens_design, "../data/analysis_results/sensitivity_profile_table_AssaysDesign.csv",
row.names=FALSE)


# Load antigen model profile and save those assays separately
antigenModelProfile <- read.csv("../data/analysis_results/05_characteristics_antigen_technique_assay_sensitivity_curve.csv",
                               stringsAsFactors=FALSE)

antigenModelAssays <- filter(antigenModelProfile, testName %in% missingAssays & !is.na(testName))

antigenModelAssays <- filter(antigenModelAssays, !is.na(testName) & time %in% tabulatedTimes)
antigenModelAssays$sensitivityMean <- as.character(round(antigenModelAssays$sensitivityMean*100))
antigenModelAssays$sensitivityL <- as.character(round(antigenModelAssays$sensitivityL*100))
antigenModelAssays$sensitivityQL <- as.character(round(antigenModelAssays$sensitivityQL*100))
antigenModelAssays$sensitivityH <- as.character(round(antigenModelAssays$sensitivityH*100))
antigenModelAssays$sensitivityQH <- as.character(round(antigenModelAssays$sensitivityQH*100))

antigenModelAssays$sensitivityStr <- with(antigenModelAssays,
          paste(sensitivityMean, ' (', sensitivityL, '-', sensitivityH, ')', sep=""))

assayTidySens_antigen <- select(antigenModelAssays, time, testName, sensitivityStr) %>%
  pivot_wider(., names_from='time', id_cols=testName,
                    values_from=sensitivityStr, names_prefix="")

write.csv(assayTidySens_antigen, "../data/analysis_results/sensitivity_profile_table_AssaysAntigen.csv",
row.names=FALSE)


# Get basic model fit to the last assay
basicModelAssays <- filter(basicModelProfile, testName %in% missingAssays[6] & !is.na(testName))

basicModelAssays <- filter(basicModelAssays, !is.na(testName) & time %in% tabulatedTimes)
basicModelAssays$sensitivityMean <- as.character(round(basicModelAssays$sensitivityMean*100))
basicModelAssays$sensitivityL <- as.character(round(basicModelAssays$sensitivityL*100))
basicModelAssays$sensitivityQL <- as.character(round(basicModelAssays$sensitivityQL*100))
basicModelAssays$sensitivityH <- as.character(round(basicModelAssays$sensitivityH*100))
basicModelAssays$sensitivityQH <- as.character(round(basicModelAssays$sensitivityQH*100))

basicModelAssays$sensitivityStr <- with(basicModelAssays,
          paste(sensitivityMean, ' (', sensitivityL, '-', sensitivityH, ')', sep=""))

assayTidySens_basic <- select(basicModelAssays, time, testName, sensitivityStr) %>%
  pivot_wider(., names_from='time', id_cols=testName,
                    values_from=sensitivityStr, names_prefix="")

write.csv(assayTidySens_basic, "../data/analysis_results/sensitivity_profile_table_AssaysBasic.csv",
row.names=FALSE)

