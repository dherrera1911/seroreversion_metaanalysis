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


assaySummaryDf <- read.csv("../data/analysis_results/04_parameter_summary.csv",
                           stringsAsFactors=FALSE)

assaySummaryDf$shorterNames <- assaySummaryDf$assay %>%
  stringr::str_replace(., " \\s*\\([^\\)]+\\)", "") %>%
  stringr::str_replace(., "Anti-SARS-CoV-2", "") %>%
  stringr::str_replace(., "SARS-CoV-2", "") %>%
  stringr::str_trunc(., width=25, ellipsis="")

###########################
# Whisker plots of the assay parameters in the basic fit
###########################

# Whisker plot of the slopes
slopesDf <- dplyr::filter(assaySummaryDf, !is.na(assay) &
                          .variable=="assaySlope")
sortedAssays <- slopesDf$shorterNames[order(slopesDf$paramMean, decreasing=T)]
slopesDf$shorterNames <- factor(slopesDf$shorterNames, levels=sortedAssays)
slopeHist <- slopesDf %>%
  ggplot(., aes(x=paramMean, y=shorterNames)) +
  geom_pointrange(aes(xmin=paramL, xmax=paramH)) +
  theme_bw() +
  theme(axis.title.y=element_blank()) +
  xlab("Assay slope")

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



