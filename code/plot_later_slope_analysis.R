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

twoSlopesPosterior <- read.csv("../data/analysis_results/06_late_slope_posterior_samples.csv",
                 stringsAsFactors=FALSE)

twoSlopesMeans <- dplyr::filter(twoSlopesPosterior, is.na(assay)) %>%
  dplyr::filter(., .variable %in% c("charSlope", "meanSlopeLate")) %>%
  pivot_wider(., id_cols=.draw, names_from=.variable, values_from=.value) %>%
  dplyr::mutate(., earlyLarger=charSlope>meanSlopeLate) %>%
  summarize(., meanEarlyLarger=mean(earlyLarger))

assayList <- c("Elecsys Anti-SARS-CoV-2 Roche spike",
               "Vitros Ortho total Ig anti-spike",
               "Wantai ELISA SARS-CoV-2 Total Antibody, IgG IgM anti-RBD")

assayLateSlopeSig <- dplyr::filter(twoSlopesPosterior, !is.na(assay)) %>%
  dplyr::filter(., .variable=="assaySlopeLate") %>%
  group_by(., assay) %>%
  summarize(., positiveFraction=mean(.value>0))

## Plot whisker plots

assaySummaryDf <- read.csv("../data/analysis_results/06_parameter_summary.csv",
                           stringsAsFactors=FALSE)

assaySummaryDf$shorterNames <- assaySummaryDf$assay %>%
  stringr::str_replace(., " \\s*\\([^\\)]+\\)", "") %>%
  stringr::str_replace(., "Anti-SARS-CoV-2", "") %>%
  stringr::str_replace(., "SARS-CoV-2", "") %>%
  stringr::str_trunc(., width=25, ellipsis="")

# Whisker plot of the EARLY slopes
slopesEarlyDf <- dplyr::filter(assaySummaryDf, !is.na(assay) &
                          .variable=="assaySlopeEarly")
sortedAssays <- slopesEarlyDf$shorterNames[order(slopesEarlyDf$paramMean, decreasing=T)]
slopesEarlyDf$shorterNames <- factor(slopesEarlyDf$shorterNames, levels=sortedAssays)

slopeEarlyHist <- slopesEarlyDf %>%
  ggplot(., aes(x=paramMean, y=shorterNames)) +
  geom_pointrange(aes(xmin=paramL, xmax=paramH)) +
  theme_bw() +
  theme(axis.title.y=element_blank()) +
  xlab("Early slope") +
  xlim(-0.6, 1.7)


# Whisker plot of the EARLY slopes
slopesLateDf <- dplyr::filter(assaySummaryDf, !is.na(assay) &
                          .variable=="assaySlopeLate")
sortedAssays <- slopesLateDf$shorterNames[order(slopesLateDf$paramMean, decreasing=T)]
slopesLateDf$shorterNames <- factor(slopesLateDf$shorterNames, levels=sortedAssays)

slopeLateHist <- slopesLateDf %>%
  ggplot(., aes(x=paramMean, y=shorterNames)) +
  geom_pointrange(aes(xmin=paramL, xmax=paramH)) +
  theme_bw() +
  theme(axis.title.y=element_blank()) +
  xlab("Late slope") +
  xlim(-0.6, 1.7)


# Put together both whisker plots in a figure
paramPlot <- ggarrange(plotlist=list(slopeEarlyHist, slopeLateHist),
                       ncol=2)
ggsave("../data/figures/whisker_plots_two_slopes.png", paramPlot,
  units="cm", width=22, height=9)



