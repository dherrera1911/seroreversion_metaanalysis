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
source("./functions_seroreversion_fit_analysis.R")

regLineSize=1
ribbonAlpha=0.2
locLineSize=0.4
locLineAlpha=0.4
locPointSize=1.6
sampleLineSize=0.1
sampleLineAlpha=0.1
dataAlpha <- 0.5


# Load raw data
seroFitted <- read.csv("../data/processed_data/PCR_to_serotest_all.csv",
                    stringsAsFactors=FALSE)

# Load validation results
validationDf <- read.csv("../data/analysis_results/04_predicted_sensitivities_grouped_CV.csv",
                         stringsAsFactors=FALSE) %>%
dplyr::mutate(., inInterval=(nSeropositives>=predictedPosL) &
              (nSeropositives<=predictedPosH))

#####################
# Sensitivity profiles fitted by the basic model
#####################
# Load sensitivity profile of basic model
sensProfileDf <- read.csv("../data/analysis_results/04_assay_sensitivity_curve.csv",
                               stringsAsFactors=FALSE)

# Plot the fitting results with the raw data
nPlots <- 3
nTests <- length(unique(seroFitted$testName))
testLists <- split(c(1:nTests), cut(seq_along(c(1:nTests)), nPlots,
                                    labels = FALSE)) 
sensitivityPlot <- list()
for (np in c(1:nPlots)) {
  testsPlot <- unique(seroFitted$testName)[testLists[[np]]]
  sensitivityPlot[[np]] <- dplyr::filter(validationDf, testName %in% testsPlot) %>%
    ggplot(aes(x=testTime, y=sensitivityMean, color=citationID,
               shape=inInterval)) +
    geom_pointrange(aes(ymin=sensitivityL, ymax=sensitivityH),
                    position=position_jitter(width=0.15, height=0)) +
    geom_line(data=dplyr::filter(sensProfileDf,
                                 !is.na(testName) & (testName %in% testsPlot)),
              aes(x=time, y=sensitivityMean*100), color="black", linetype="solid",
              size=regLineSize, inherit.aes=FALSE) +
    geom_ribbon(data=dplyr::filter(sensProfileDf,
                                   !is.na(testName) & testName %in% testsPlot),
                aes(x=time, ymin=sensitivityL*100, ymax=sensitivityH*100),
                alpha=ribbonAlpha, colour=NA, show.legend=FALSE,
                inherit.aes=FALSE) +
    facet_wrap(.~testName, ncol=3, labeller = label_wrap_gen(width=35)) +
    theme_bw() +
    theme(legend.position="top") +
    guides(color=FALSE, shape=guide_legend("In validation interval")) +
    xlim(0, 15) +
    xlab("Diagnosis to test (months)") +
    ylab("Sensitivity (%)")
ggsave(paste("../data/figures/seroreversion_fit", np, ".png", sep=""),
       sensitivityPlot[[np]], units="cm", width=24, height=30)
}


####
antigenProfile <- read.csv("../data/analysis_results/05_characteristics_fullModel_assay_sensitivity_curve.csv",
                           stringsAsFactors=FALSE)

sensitivityCompPlot <- list()
for (np in c(1:nPlots)) {
  testsPlot <- unique(seroFitted$testName)[testLists[[np]]]
  sensitivityCompPlot[[np]] <- sensitivityPlot[[np]] +
    geom_line(data=dplyr::filter(antigenProfile,
                                 !is.na(testName) & (testName %in% testsPlot)),
              aes(x=time, y=sensitivityMean*100), color="red", linetype="solid",
              size=regLineSize, inherit.aes=FALSE) +
    geom_ribbon(data=dplyr::filter(antigenProfile,
                                   !is.na(testName) & testName %in% testsPlot),
                aes(x=time, ymin=sensitivityL*100, ymax=sensitivityH*100),
                alpha=ribbonAlpha, colour=NA, show.legend=FALSE,
                inherit.aes=FALSE, fill="red")
ggsave(paste("../data/figures/seroreversion_fit_characteristics", np, ".png", sep=""),
       sensitivityCompPlot[[np]], units="cm", width=24, height=30)
}


validationDf2 <- read.csv("../data/analysis_results/05_characteristics_antigen_grouped_CV.csv",
                         stringsAsFactors=FALSE) %>%
  dplyr::mutate(., inInterval=(nSeropositives>=predictedPosL) &
                (nSeropositives<=predictedPosH))


