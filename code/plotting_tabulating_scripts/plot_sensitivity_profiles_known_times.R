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
source("../functions_auxiliary.R")
source("../functions_seroreversion_fit_analysis.R")

regLineSize=1
ribbonAlpha=0.2
locLineSize=0.4
locLineAlpha=0.4
locPointSize=2
sampleLineSize=0.2
sampleLineAlpha=0.1
dataAlpha <- 0.5

# Load raw data
seroFitted <- read.csv("../../data/processed_data/PCR_to_serotest_all.csv",
                    stringsAsFactors=FALSE)

# Put characteristics columns to raw data
seroFitted$RBD <- stringr::str_detect(seroFitted$antigenTarget, "RBD")
seroFitted$S <- stringr::str_detect(seroFitted$antigenTarget, "S") &
  !seroFitted$RBD
seroFitted$N <- stringr::str_detect(seroFitted$antigenTarget, "N") &
  !(seroFitted$RBD | seroFitted$S)

seroFitted$LFA <- stringr::str_detect(seroFitted$technique, "LFIA")
seroFitted$Indirect <- ! seroFitted$LFA &
  stringr::str_detect(seroFitted$design, "indirect")
seroFitted$Direct <- ! seroFitted$LFA &
  stringr::str_detect(seroFitted$design, "sandwich")
seroFitted$Neutralization <- ! seroFitted$LFA &
  stringr::str_detect(seroFitted$design, "competitive")

seroFitted <- dplyr::filter(seroFitted, timeKnown)

binaryString <- c("No", "Yes")

# Load validation results
validationDf <- read.csv("../../data/analysis_results/05_characteristics_fullModel_grouped_CV.csv",
                         stringsAsFactors=FALSE) %>%
dplyr::mutate(., inInterval=(nSeropositives>=predictedPosL) &
              (nSeropositives<=predictedPosH)) %>%
dplyr::filter(., timeKnown)

#################################
#################################
# Individual assays sensitivity profiles (Figure S1)
#################################
#################################

nPlots <- 3 # Number of sheets into which to divide the assay-plots
nTests <- length(unique(seroFitted$testName))
testLists <- split(c(1:nTests), cut(seq_along(c(1:nTests)), nPlots,
                                    labels = FALSE)) 

#################################
#################################
# Average sensitivity profiles of assay types
# Makes Figures 2, S3B, S4B, S5B, S6B
#################################
#################################


############
# Full Model (Fig 2)
############

fullModelProfile <- read.csv("../../data/analysis_results/08_characteristics_fullModel_assay_sensitivity_curve_known_times.csv",
                           stringsAsFactors=FALSE) %>%
  dplyr::mutate(., averageSens=is.na(testName))

# additional column to organize plots
fullModelProfile$Antigen[with(fullModelProfile, N)] <- "Nucleocapsid"
fullModelProfile$Antigen[with(fullModelProfile, S)] <- "Spike"
fullModelProfile$Antigen[with(fullModelProfile, RBD)] <- "Receptor-binding domain"

levelOrder <- c("Nucleocapsid", "Spike", "Receptor-binding domain")
fullModelProfile$Antigen <- factor(fullModelProfile$Antigen, levels=levelOrder)

binaryString <- c("No", "Yes")
fullModelProfile$LFA <- binaryString[as.integer(fullModelProfile$LFA)+1]
fullModelProfile$Direct <- binaryString[as.integer(fullModelProfile$Direct)+1]
fullModelProfile$Neutralization <- binaryString[as.integer(fullModelProfile$Neutralization)+1]
fullModelProfile$Indirect <- with(fullModelProfile, LFA=="No" & Direct=="No" & Neutralization=="No")
fullModelProfile$Indirect <- binaryString[as.integer(fullModelProfile$Indirect)+1]

fullModelProfile$Design[with(fullModelProfile, LFA=="Yes")] <- "LFA"
fullModelProfile$Design[with(fullModelProfile, LFA=="No" & Direct=="No" &
                             Neutralization=="No")] <- "Quantitative-Indirect"
fullModelProfile$Design[with(fullModelProfile, Direct=="Yes")] <- "Quantitative-Direct"
fullModelProfile$Design[with(fullModelProfile, Neutralization=="Yes")] <- "Quantitative-Competitive"

levelOrder2 <- c("LFA", "Quantitative-Indirect", "Quantitative-Direct", "Quantitative-Competitive")
fullModelProfile$Design <- factor(fullModelProfile$Design, levels=levelOrder2)

fullModelAverages <- filter(fullModelProfile, averageSens)
fullModelAssays <- filter(fullModelProfile, !averageSens)

testfullModel <- dplyr::select(fullModelAssays, testName, Antigen, LFA, Indirect,
                               Direct, Neutralization, Design) %>%
  unique()


fullModelPoints <- seroFitted %>%
  filter(., testName %in% fullModelProfile$testName) %>%
  select(., -LFA, -Indirect, -Direct, -Neutralization) %>%
  merge(., testfullModel) %>%
  dplyr::mutate(., time=testTime)

fullModelProfilePlot <- fullModelAverages %>%
  ggplot(., aes(x=time, y=sensitivityMean*100)) +
  geom_line(size=regLineSize, color="red") +
  geom_ribbon(aes(ymin=sensitivityL*100, ymax=sensitivityH*100),
              alpha=ribbonAlpha, color=NA, fill="red") +
  geom_point(data=fullModelPoints,
             aes(x=time, y=sensitivityMean, size=sqrt(nSamples)/10),
             color="black", alpha=0.2) +
  geom_line(data=fullModelAssays, aes(x=time, y=sensitivityMean*100, group=testName),
            color="black", alpha=0.7) +
  #facet_wrap(~Antigen+Design, labeller=label_both) +
  facet_grid(Antigen~Design, switch="y") +
  theme_bw() +
  #theme(legend.position="top", strip. background = element_blank()) +
  theme(legend.position="top", strip.placement = "outside") +
  guides(color=FALSE, size=FALSE, fill=FALSE) +
  xlim(0, 15) +
  xlab("Diagnosis to test (months)") +
  ylab("Sensitivity (%)")


ggsave("../../data/figures/characteristics_profiles_fullModel_known_times.png", fullModelProfilePlot,
       units="cm", width=18, height=15)

