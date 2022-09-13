####################################################
#
# Plot the whisker plot that shows the slope and
# intercept for each assay, obtained from the basic
# model (without test characteristics).
# Generates the plots for figures 1A, 1B
# 
# The inputs for plotting and statistical analysis are
# the posterior samples of the fitted models, which
# are in .csv files with the names
# '05_characteristics_XXX_posterior_samples', where
# XXX = the model identifier. These are generated in
# script 05_characteristics_sensitivity_analysis.R
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

# Load the pre-computed summary of parameters
assaySummaryDf <- read.csv("../data/analysis_results/04_parameter_summary.csv",
                           stringsAsFactors=FALSE) %>%
  dplyr::filter(., !is.na(assay))

assayChars <- read.csv("../data/raw_data/assay_characteristics.csv",
                       stringsAsFactors=FALSE)

assayInd <- match(assaySummaryDf$assay, assayChars$test_name_long)
assaySummaryDf$assay2 <- assayChars$test_name_simplified[assayInd]
assaySummaryDf$Antigen <- assayChars$antigen.target[assayInd]
assaySummaryDf$Antigen[assaySummaryDf$Antigen == "--"] <- NA
assaySummaryDf$Design <- assayChars$test.design[assayInd]
assaySummaryDf$Design[assaySummaryDf$Design == ""] <- NA
assaySummaryDf$Technique <- assayChars$assay.type[assayInd]

############
# Add binary columns to dataframe, to work as indicator variables
############
assaySummaryDf$RBD <- stringr::str_detect(assaySummaryDf$Antigen, "RBD")
assaySummaryDf$S <- stringr::str_detect(assaySummaryDf$Antigen, "S") &
  !assaySummaryDf$RBD
assaySummaryDf$N <- stringr::str_detect(assaySummaryDf$Antigen, "N") &
  !(assaySummaryDf$RBD | assaySummaryDf$S)
assaySummaryDf$Antigen[assaySummaryDf$RBD] <- "Receptor-binding domain"
assaySummaryDf$Antigen[assaySummaryDf$S] <- "Spike"
assaySummaryDf$Antigen[assaySummaryDf$N] <- "Nucleocapsid"
assaySummaryDf$Antigen[is.na(assaySummaryDf$Antigen)] <- "Unknown"
assaySummaryDf$Antigen <- factor(assaySummaryDf$Antigen, levels=c("Nucleocapsid", "Spike",
                                            "Receptor-binding domain", "Unknown")) 
assaySummaryDf$LFA <- stringr::str_detect(assaySummaryDf$Technique, "LFIA")
binaryString <- c("No", "Yes")
assaySummaryDf$LFA <- binaryString[as.integer(assaySummaryDf$LFA)+1]
assaySummaryDf$Sandwich <- stringr::str_detect(assaySummaryDf$Design, "sandwich")
assaySummaryDf$Sandwich <- binaryString[as.integer(assaySummaryDf$Sandwich)+1]
assaySummaryDf$Sandwich[is.na(assaySummaryDf$Sandwich)] <- "Unknown"

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
  ggplot(., aes(x=paramMean, y=shorterNames, color=Antigen, shape=Sandwich,
                linetype=LFA)) +
  geom_pointrange(aes(xmin=paramL, xmax=paramH)) +
  #facet_wrap(~Technique, ncol=1, scales="free_y") +
  theme_bw() +
  theme(axis.title.y=element_blank()) +
  xlab("Assay slope") +
  geom_vline(xintercept=0) +
  labs(linetype="Lateral Flow Assay") +
  annotate("rect", xmin=-0.9, xmax=0.7, ymin=0.5, ymax=3.5,
           alpha=.1,fill = "blue")

ggsave("../data/figures/whisker_plots_slopes_basic_model2.png", slopeHist,
  units="cm", width=20, height=19)

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
ggsave("../data/figures/whisker_plots_basic_model.png", paramPlot,
  units="cm", width=22, height=14)


