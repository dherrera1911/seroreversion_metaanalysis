####################################################
# Extract some statistics, and make some plots of how
# the specificity of assays depends on assay characteristics
#
# Generates plot S7, and the accompanying statistics of
# Supplementary Section G.
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
library(rstan)
library(bayesplot)
library(tidybayes)
library(viridis)
library(RColorBrewer)
library(ggthemes)
library(lemon)
library(cowplot)
library(ggpubr)
source("../functions_auxiliary.R")
source("../functions_seroreversion_fit_analysis.R")


# read specificity data csv
specificityDf <- read.csv("../../data/processed_data/assay_specificity.csv",
                          stringsAsFactors=FALSE)

# Put characteristics columns to raw data
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

specificityDf$Antigen[with(specificityDf, N)] <- "Nucleocapsid"
specificityDf$Antigen[with(specificityDf, S)] <- "Spike"
specificityDf$Antigen[with(specificityDf, RBD)] <- "Receptor-binding domain"

levelOrder <- c("Nucleocapsid", "Spike", "Receptor-binding domain")
specificityDf$Antigen <- factor(specificityDf$Antigen, levels=levelOrder)

binaryString <- c("No", "Yes")
specificityDf$LFA <- binaryString[as.integer(specificityDf$LFA)+1]
specificityDf$Direct <- binaryString[as.integer(specificityDf$Direct)+1]
specificityDf$Neutralization <- binaryString[as.integer(specificityDf$Neutralization)+1]
specificityDf$Indirect <- with(specificityDf, LFA=="No" & Direct=="No" & Neutralization=="No")
specificityDf$Indirect <- binaryString[as.integer(specificityDf$Indirect)+1]

specificityDf$Design[with(specificityDf, LFA=="Yes")] <- "LFA"
specificityDf$Design[with(specificityDf, LFA=="No" & Direct=="No" &
                             Neutralization=="No")] <- "Quantitative-Indirect"
specificityDf$Design[with(specificityDf, Direct=="Yes")] <- "Quantitative-Direct"
specificityDf$Design[with(specificityDf, Neutralization=="Yes")] <- "Quantitative-Competitive"

levelOrder2 <- c("LFA", "Quantitative-Indirect", "Quantitative-Direct", "Quantitative-Competitive")
specificityDf$Design <- factor(specificityDf$Design, levels=levelOrder2)


# Get estimated mean specificities for each kind of assay
specProfiles <- read.csv("../../data/analysis_results/10_specificity_fullModel_specificity_averages.csv",
                        stringsAsFactors=FALSE)
# Add antigen string 
specProfiles$Antigen[with(specProfiles, N)] <- "Nucleocapsid"
specProfiles$Antigen[with(specProfiles, S)] <- "Spike"
specProfiles$Antigen[with(specProfiles, RBD)] <- "Receptor-binding domain"
specProfiles$Antigen <- factor(specProfiles$Antigen, levels=levelOrder)
# Add design string
specProfiles$Design[with(specProfiles, LFA)] <- "LFA"
specProfiles$Design[with(specProfiles, !LFA & !Direct &
                             !Neutralization)] <- "Quantitative-Indirect"
specProfiles$Design[with(specProfiles, Direct)] <- "Quantitative-Direct"
specProfiles$Design[with(specProfiles, Neutralization)] <- "Quantitative-Competitive"
specProfiles$Design <- factor(specProfiles$Design, levels=levelOrder2)


# Plot the specificity of assays, indicating assay characteristics
specificityPlot <- dplyr::filter(specificityDf, !is.na(Antigen) &
                                 !is.na(Design)) %>%
  ggplot(., aes(Antigen, specificityMean)) +
  facet_grid(~Design, switch="y") +
  geom_point() +
  geom_jitter(width=0.2) +
  geom_pointrange(data=specProfiles,
                  aes(y=specificityMean*100, ymin=specificityL*100,
                      ymax=specificityH*100), color="red", fatten=2,
                  size=1) +
  theme_bw() +
  ylab("Specificity (%)")

ggsave("../../data/figures/specificity_fit.png", specificityPlot,
       units="cm", width=20, height=12)

# Do some statistics on the fitted characteristics effects
specPosteriors <- read.csv("../../data/analysis_results/10_specificity_fullModel_posterior_samples.csv",
                           stringsAsFactors=FALSE)

### Full model
fullModelWide <- specPosteriors %>%
  dplyr::filter(., !is.na(par)) %>%
  pivot_wider(., id_cols=.draw, names_from=parName, values_from=.value)

prob_LFAlarger0 <- mean(fullModelWide$LFA>0)
prob_Directlarger0 <- mean(fullModelWide$Direct>0)
prob_Neutralizationlarger0 <- mean(fullModelWide$Neutralization>0)
prob_SlargerN <- with(fullModelWide, mean(S>N))
prob_RBDlargerS <- with(fullModelWide, mean(RBD>S))
prob_RBDlargerN <- with(fullModelWide, mean(RBD>N))


### See specificity of Orient Gene Biotech
#assayName <- "Orient Gene Biotech COVID-19 IgG/IgM Rapid Test anti-spike"
#testTrace <- dplyr::filter(specPosteriors, assay==assayName)
#lala <- inverse_logit(testTrace$.value)

