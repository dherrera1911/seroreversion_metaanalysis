####################################################
#
# Tabulate the sensitivity across time, both for the 'averages'
# of the different types of assays, and the assay-specific
# sensitivity profiles, rounding up the percentages and in a
# tidy form.
# Generates tables 1 and S1.
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
library(bayesplot)
library(tidybayes)
library(viridis)
library(RColorBrewer)
library(gridExtra)
library(ggthemes)
library(lemon)
library(cowplot)
library(ggpubr)
source("../functions_auxiliary.R")
source("../functions_seroreversion_fit_analysis.R")

#################################
#################################
# Tabulate the mean sensitivity profiles of assay types
#################################
#################################

# Load sensitivity profile of basic model
basicModelProfile <- read.csv("../../data/analysis_results/04_assay_sensitivity_curve.csv",
                               stringsAsFactors=FALSE)
fullModelProfile <- read.csv("../../data/analysis_results/05_characteristics_fullModel_assay_sensitivity_curve.csv",
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


tidySens <- select(meanProfiles, time, antigen, LFA, Direct, Neutralization, sensitivityStr) %>%
  pivot_wider(., names_from='time', id_cols=c(antigen, LFA, Direct, Neutralization),
                    values_from=sensitivityStr, names_prefix="")

write.csv(tidySens, "../data/analysis_results/sensitivity_profile_table.csv",
row.names=FALSE)

#################################
#################################
# Tabulate the assay-specific sensitivity profiles
#################################
#################################

assayChars <- read.csv("../../data/raw_data/assay_characteristics.csv", stringsAsFactors=FALSE)

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

assayNamesTable <- dplyr::filter(assayChars, test_name_long %in% assayTidySens$testName) %>%
  dplyr::select(., test_name_long, test_name_simplified) %>%
  dplyr::rename(., testName=test_name_long)

assayTidySens <- merge(assayNamesTable, assayTidySens)

write.csv(assayTidySens, "../../data/analysis_results/sensitivity_profile_table_AssaysFull.csv",
  row.names=FALSE)


# Get assays missing from the full model
missingAssaysInd <- which(!basicModelProfile$testName %in% assayProfiles$testName)
missingAssays <- unique(basicModelProfile$testName[missingAssaysInd])

# Get the characteristics of missing assays, to see what model
# to use for their sensitivity
missingAssayChars <- dplyr::filter(assayChars, test_name_long %in% missingAssays) %>%
  dplyr::select(., test_name_long, antigen.target, test.design, assay.type)

# Load technique model profile and save those assays separately
designModelProfile <- read.csv("../../data/analysis_results/05_characteristics_design_assay_sensitivity_curve.csv",
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

assayNamesTable <- dplyr::filter(assayChars, test_name_long %in% assayTidySens_design$testName) %>%
  dplyr::select(., test_name_long, test_name_simplified) %>%
  dplyr::rename(., testName=test_name_long)

assayTidySens_design <- merge(assayNamesTable, assayTidySens_design)

write.csv(assayTidySens_design, "../../data/analysis_results/sensitivity_profile_table_AssaysDesign.csv",
row.names=FALSE)


# Load antigen model profile and save those assays separately
antigenModelProfile <- read.csv("../../data/analysis_results/05_characteristics_antigen_technique_assay_sensitivity_curve.csv",
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

assayNamesTable <- dplyr::filter(assayChars, test_name_long %in% assayTidySens_antigen$testName) %>%
  dplyr::select(., test_name_long, test_name_simplified) %>%
  dplyr::rename(., testName=test_name_long)

assayTidySens_antigen <- merge(assayNamesTable, assayTidySens_antigen)


write.csv(assayTidySens_antigen, "../../data/analysis_results/sensitivity_profile_table_AssaysAntigen.csv",
row.names=FALSE)

# Make a table with tidy assay characteristics
assayProfiles$Antigen[with(assayProfiles, N)] <- "Nucleocapsid"
assayProfiles$Antigen[with(assayProfiles, S)] <- "Spike"
assayProfiles$Antigen[with(assayProfiles, RBD)] <- "Receptor-binding domain"

assayProfiles$Design[assayProfiles$LFA] <- "LFA"
assayProfiles$Design[with(assayProfiles, !LFA & !Direct=="No" &
                             !Neutralization)] <- "Quantitative-Indirect"
assayProfiles$Design[assayProfiles$Direct] <- "Quantitative-Direct"
assayProfiles$Design[assayProfiles$Neutralization] <- "Quantitative-Competitive"

assayCharsTidy <- dplyr::select(assayProfiles, testName, Antigen, Design) %>%
  unique(.)

write.csv(assayCharsTidy, "../../data/analysis_results/tidy_assay_characteristics.csv",
row.names=FALSE)


