####################################################
#
# FIGURES 2, S5 of associated paper
# "Dynamics of SARS-CoV-2 seroassay sensitivity: a systematic
# review and modeling study"
#
# Plot the whisker plot that shows the slope and
# intercept for each assay, obtained from the basic
# model (without test characteristics). Can do with either
# all data (Fig 2), or data with only known diagnosis-to-test
# times (Fig S5).
# 
# The inputs for plotting and statistical analysis are
# the posterior samples of the fitted models, which
# are in .csv files with the names
# '04_parameter_summary.csv'.
#
# Script authored by Daniel Herrera-Esposito.
# For questions, contact me at dherrera1911[at]gmail.com
# 
# Final version revised 10/03/2023
# 
###############################

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
source("../functions_auxiliary.R")
source("../functions_seroreversion_fit_analysis.R")

regLineSize=1
ribbonAlpha=0.2
locLineSize=0.4
locLineAlpha=0.4
locPointSize=1.6
sampleLineSize=0.1
sampleLineAlpha=0.1
dataAlpha <- 0.5

# Load the pre-computed summary of parameters
knownTimesOnly <- TRUE
if (!knownTimesOnly) {
  assaySummaryDf <- read.csv("../../data/analysis_results/04_parameter_summary.csv",
                             stringsAsFactors=FALSE) %>%
    dplyr::filter(., !is.na(assay))
  knownTimesStr <- ""
} else {
  assaySummaryDf <- read.csv("../../data/analysis_results/04_parameter_summary_known_times.csv",
                             stringsAsFactors=FALSE) %>%
    dplyr::filter(., !is.na(assay))
  knownTimesStr <- "_known_times"
}

assayChars <- read.csv("../../data/raw_data/assay_characteristics.csv",
                       stringsAsFactors=FALSE)

# Add columns to assay summary, with characteristics from the assay chars table
assayInd <- match(assaySummaryDf$assay, assayChars$test_name_long)
assaySummaryDf$assay2 <- assayChars$test_name_simplified[assayInd]
assaySummaryDf$Antigen <- assayChars$antigen.target[assayInd]
assaySummaryDf$Antigen[assaySummaryDf$Antigen == "--"] <- NA
assaySummaryDf$Binding <- assayChars$test.design[assayInd]
assaySummaryDf$Binding[assaySummaryDf$Design == ""] <- NA
assaySummaryDf$Technique <- assayChars$assay.type[assayInd]

############
# Add binary columns to dataframe, to work as indicator variables
############
assaySummaryDf$RBD <- stringr::str_detect(assaySummaryDf$Antigen, "RBD")
assaySummaryDf$S <- stringr::str_detect(assaySummaryDf$Antigen, "S") &
  !assaySummaryDf$RBD
assaySummaryDf$N <- stringr::str_detect(assaySummaryDf$Antigen, "N") &
  !(assaySummaryDf$RBD | assaySummaryDf$S)
assaySummaryDf$Direct <- assaySummaryDf$Binding == "sandwich"
assaySummaryDf$Indirect <- assaySummaryDf$Binding == "indirect"
assaySummaryDf$Neutralization <- assaySummaryDf$Binding == "competitive"
assaySummaryDf$LFA <- stringr::str_detect(assaySummaryDf$Technique, "LFIA")

# Turn binary columns into names for the plot
assaySummaryDf$Antigen[assaySummaryDf$RBD] <- "Receptor-binding domain"
assaySummaryDf$Antigen[assaySummaryDf$S] <- "Spike"
assaySummaryDf$Antigen[assaySummaryDf$N] <- "Nucleocapsid"
assaySummaryDf$Antigen[is.na(assaySummaryDf$Antigen)] <- "Unknown"
assaySummaryDf$Antigen <- factor(assaySummaryDf$Antigen, levels=c("Nucleocapsid", "Spike",
                                            "Receptor-binding domain", "Unknown")) 

assaySummaryDf$Design[with(assaySummaryDf, LFA)] <- "LFA"
assaySummaryDf$Design[with(assaySummaryDf, !LFA & Indirect)] <- "Quantitative-Indirect"
assaySummaryDf$Design[with(assaySummaryDf, !LFA & Direct)] <- "Quantitative-Direct"
assaySummaryDf$Design[with(assaySummaryDf, !LFA & Neutralization)] <- "Quantitative-Competitive"
assaySummaryDf$Design[is.na(assaySummaryDf$Design)] <- "Unknown"
levelOrder2 <- c("LFA", "Quantitative-Indirect", "Quantitative-Direct", "Quantitative-Competitive",
  "Unknown")
assaySummaryDf$Design <- factor(assaySummaryDf$Design, levels=levelOrder2)


assaySummaryDf$shorterNames <- assaySummaryDf$assay2 %>%
  stringr::str_replace(., " \\s*\\([^\\)]+\\)", "") %>%
  stringr::str_replace(., "Anti-SARS-CoV-2", "") %>%
  stringr::str_replace(., "SARS-CoV-2", "") %>%
  stringr::str_trunc(., width=40, ellipsis="")

###########################
# Whisker plots of the assay parameters in the basic fit
###########################

# Whisker plot of the slopes
slopesDf <- dplyr::filter(assaySummaryDf, !is.na(assay) &
                          .variable=="assaySlope")
sortedAssays <- slopesDf$shorterNames[order(slopesDf$paramMean, decreasing=T)]
slopesDf$shorterNames <- factor(slopesDf$shorterNames, levels=sortedAssays)
slopeHist <- slopesDf %>%
  ggplot(., aes(x=paramMean, y=shorterNames, color=Antigen, shape=Design)) +
  geom_pointrange(aes(xmin=paramL, xmax=paramH)) +
  #facet_wrap(~Technique, ncol=1, scales="free_y") +
  theme_bw() +
  theme(axis.title.y=element_blank()) +
  scale_shape_manual(values=c(5, 16, 15, 17, 4))+
  scale_color_manual(values=c('brown1', 'darkgoldenrod2', 'deepskyblue', 'azure4'))+
  xlab("Assay slope") +
  geom_vline(xintercept=0) +
  labs(linetype="Lateral Flow Assay") +
  annotate("rect", xmin=-0.9, xmax=0.7, ymin=0.5, ymax=3.5,
           alpha=.1,fill = "blue")

fileName <- paste("../../data/figures/whisker_plots_slopes_basic_model", knownTimesStr, ".png", sep="")
ggsave(fileName, slopeHist, units="cm", width=20, height=19)

# Whisker plot of the intercepts
interceptsDf <- dplyr::filter(assaySummaryDf, !is.na(assay) &
                          .variable=="assayIntercept")
sortedAssays <- interceptsDf$shorterNames[order(interceptsDf$paramMean, decreasing=T)]
interceptsDf$shorterNames <- factor(interceptsDf$shorterNames, levels=sortedAssays)
interceptHist <- interceptsDf %>%
  ggplot(., aes(x=paramMean, y=shorterNames)) +
  geom_pointrange(aes(xmin=paramL, xmax=paramH)) +
  theme_bw() +
  theme(axis.title.y=element_blank()) +
  xlab("Assay intercept")

# Put together both whisker plots in a figure
paramPlot <- ggarrange(plotlist=list(slopeHist, interceptHist),
                       ncol=2)
fileName <- paste("../../data/figures/whisker_plots_basic_model", knownTimesStr, ".png", sep="")
ggsave(fileName, slopeHist, units="cm", width=20, height=19)
ggsave(fileName, paramPlot, units="cm", width=22, height=14)


