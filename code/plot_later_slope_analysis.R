####################################################
#
# Plot the whisker plots of the early and late slopes.
# Generates the plot of Figure 7. 
#
# Uses the posterior samples of the parameters of
# the two-slope model, fitted in script
# 06_positive_slopes_analysis.R
#
# Also computes some significance numbers
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
source("./functions_auxiliary.R")

# Load the posterior samples data
twoSlopesPosterior <- read.csv("../data/analysis_results/06_late_slope_posterior_samples.csv",
                 stringsAsFactors=FALSE)

# Compute the mean of the early and late slopes
twoSlopesMeans <- dplyr::filter(twoSlopesPosterior, is.na(assay)) %>%
  dplyr::filter(., .variable %in% c("charSlope", "meanSlopeLate")) %>%
  pivot_wider(., id_cols=.draw, names_from=.variable, values_from=.value) %>%
  dplyr::mutate(., earlyLarger=charSlope>meanSlopeLate) %>%
  summarize(., meanEarlyLarger=mean(earlyLarger))

# Get the proportion of samples that are >0 for the late slope of
# each test.
assayLateSlopeSig <- dplyr::filter(twoSlopesPosterior, !is.na(assay)) %>%
  dplyr::filter(., .variable=="assaySlopeLate") %>%
  group_by(., assay) %>%
  summarize(., positiveFraction=mean(.value>0))

# Do some editing of assay names to better fit the plot
assaySummaryDf <- read.csv("../data/analysis_results/06_parameter_summary.csv",
                           stringsAsFactors=FALSE)

assaySummaryDf$shorterNames <- assaySummaryDf$assay %>%
  stringr::str_replace(., " \\s*\\([^\\)]+\\)", "") %>%
  stringr::str_replace(., "Anti-SARS-CoV-2", "") %>%
  stringr::str_replace(., "SARS-CoV-2", "") %>%
  stringr::str_trunc(., width=25, ellipsis="")

# Whisker plot of the EARLY slopes (Figure 7A)
slopesEarlyDf <- dplyr::filter(assaySummaryDf, !is.na(assay) &
                          .variable=="assaySlopeEarly")
# sort assays by value, for nice increasing plot
sortedAssays <- slopesEarlyDf$shorterNames[order(slopesEarlyDf$paramMean, decreasing=T)]
slopesEarlyDf$shorterNames <- factor(slopesEarlyDf$shorterNames, levels=sortedAssays)

slopeEarlyHist <- slopesEarlyDf %>%
  ggplot(., aes(x=paramMean, y=shorterNames)) +
  geom_pointrange(aes(xmin=paramL, xmax=paramH)) +
  theme_bw() +
  theme(axis.title.y=element_blank()) +
  xlab("Early slope") +
  geom_vline(xintercept=0)


# Whisker plot of the LATE slopes (Figure 7B)
slopesLateDf <- dplyr::filter(assaySummaryDf, !is.na(assay) &
                          .variable=="assaySlopeLate")
# sort assays by value, for nice increasing plot
sortedAssays <- slopesLateDf$shorterNames[order(slopesLateDf$paramMean, decreasing=T)]
slopesLateDf$shorterNames <- factor(slopesLateDf$shorterNames, levels=sortedAssays)

slopeLateHist <- slopesLateDf %>%
  ggplot(., aes(x=paramMean, y=shorterNames)) +
  geom_pointrange(aes(xmin=paramL, xmax=paramH)) +
  theme_bw() +
  theme(axis.title.y=element_blank()) +
  xlab("Late slope") +
  geom_vline(xintercept=0)


# Put together both whisker plots in a figure and save
paramPlot <- ggarrange(plotlist=list(slopeEarlyHist, slopeLateHist),
                       ncol=2)
ggsave("../data/figures/whisker_plots_two_slopes.png", paramPlot,
  units="cm", width=22, height=9)

