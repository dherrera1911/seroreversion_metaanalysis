####################################################
#
# Plot the density plots of the model parameters related to
# assay characteristics. Generates plots for Fig S3 and S4.
# 
# Also, tris script computes the statistical comparisons
# using the model posteriors. These are reported throughout
# the text E.g. the proportion of N < RBD
# 
# The inputs for plotting and statistical analysis are
# the posterior samples of the fitted models, which
# are in .csv files with the names
# '05_characteristics_XXX_posterior_samples', where
# XXX = the model identifier. These are generated in
# script 05_characteristics_sensitivity_analysis.R
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

# Some plotting parameters
regLineSize=1
ribbonAlpha=0.2
locLineSize=0.4
locLineAlpha=0.4
locPointSize=1.6
sampleLineSize=0.1
sampleLineAlpha=0.1
dataAlpha <- 0.5

knownTimesOnly <- FALSE

# Vector indicating what models to analyze (i.e. what characteristics they included)
if (!knownTimesOnly) {
  modelNames <- c("antigen", "antibody", "technique", "fullModel")
} else {
  modelNames <- c("fullModel")
}

########
### Loop below loads each models posteriors, and makes/saves the density plots.
########

charPlot <- list()
slopePosteriors <- list()
for (m in c(1:length(modelNames))) {
  mN <- modelNames[m]
  # COMMENT DEPENDING ON WHAT FIT TO USE, ALL DATA OR KNOWN TIMES
  if (!knownTimesOnly) {
    posteriorsName <- paste("../../data/analysis_results/05_characteristics_",
                            mN, "_posterior_samples.csv", sep="")
  } else {
    posteriorsName <- paste("../data/analysis_results/08_characteristics_",
                            mN, "_posterior_samples_known_times.csv", sep="")
  }
  slopePosteriors[[mN]] <- read.csv(posteriorsName, stringsAsFactors=FALSE) %>%
    dplyr::filter(., !is.na(parName))
  if (mN == "antigen") {
    slopePosteriors[[mN]]$parName <- factor(slopePosteriors[[mN]]$parName,
                                            levels=c("N", "S", "RBD"))
  }
  if (mN == "fullModel") {
    slopePosteriors[[mN]]$isAntigen <- slopePosteriors[[mN]]$parName %in%
      c("S", "N", "RBD")
    levelOrder <- c("N", "S", "RBD", "LFA", "Direct", "Neutralization")
    slopePosteriors[[mN]]$parName <- factor(slopePosteriors[[mN]]$parName,
                                            levels=levelOrder)
  }
  charPlot[[m]] <- slopePosteriors[[mN]] %>%
    ggplot(data=., aes(x=.value)) +
    geom_density(alpha=0.3, fill="blue") +
    facet_wrap(~parName, ncol=1, scales="free_y", shrink=TRUE) +
    theme_bw() +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    scale_y_continuous(breaks = NULL) +
    #guides(fill=guide_legend("Parameter")) +
    xlab("Value") +
    ylab("Density")
  # COMMENT DEPENDING ON WHAT FIT TO USE, ALL DATA OR KNOWN TIMES
  if (!knownTimesOnly) {
  plotName <- paste("../../data/figures/characteristics_effects_plot_",
                    mN, ".png", sep="")
  } else {
    plotName <- paste("../data/figures/characteristics_effects_plot_",
                      mN, "_known_times.png", sep="")
  }
  ggsave(plotName, charPlot[[m]], units="cm", width=7, height=12)
}


#######
# Code below computes the frequencies at which different
# comparisons occur in the posterior samples in the model.
# E.g. what proportion of posterior samples have S > N.
#
# Before each block of code, a comment indicates the model
# for which comparisons will be computed
#######

comparisonsDf <- NULL

##### These are the main significance statistics reported in
##### the paper, i.e. those corresponding to the full model
#####with all characteristics
## Full model
fullModelWide <- slopePosteriors[["fullModel"]] %>%
  pivot_wider(., id_cols=.draw, names_from=parName, values_from=.value)

prob_LFAlarger0 <- mean(fullModelWide$LFA>0)
prob_Directlarger0 <- mean(fullModelWide$Direct>0)
prob_Neutralizationlarger0 <- mean(fullModelWide$Neutralization>0)
prob_Sandwichlarger0 <- mean(fullModelWide$sandwich>0)
prob_SlargerN <- with(fullModelWide, mean(S>N))
prob_RBDlarger0 <- mean(fullModelWide$RBD>0)
prob_RBDlargerS <- with(fullModelWide, mean(RBD>S))
prob_RBDlargerN <- with(fullModelWide, mean(RBD>N))

prob_RBDDirect_positive <- with(fullModelWide, mean((RBD+Direct)>0)) 

comparisonsFullModel <- c("LFA > 0", "Sandwich > 0",
                          "S > N", "RBD > 0",
                          "RBD > S", "RBD > N")

fullModelComps <- data.frame(Model="fullModel",
                           Comparison=comparisonsFullModel,
                           Frequency=c(prob_LFAlarger0,
                                       prob_Sandwichlarger0,
                                       prob_SlargerN,
                                       prob_RBDlarger0,
                                       prob_RBDlargerS,
                                       prob_RBDlargerN))

if (! knownTimesOnly) {
  ######## Below, the comparisons using models that only take
  ######## into account subsets of the data. These are presented
  ######## in Supplementary Section E in the main paper

  ### Antigen
  antigenWide <- slopePosteriors[["antigen"]] %>%
    pivot_wider(., id_cols=.draw, names_from=parName, values_from=.value)
  #
  prob_SlargerN <- with(antigenWide, mean(S>N))
  prob_RBDlargerS <- with(antigenWide, mean(RBD>S))
  prob_RBDlargerN <- with(antigenWide, mean(RBD>N))
  prob_RBDlarger0 <- mean(antigenWide$RBD>0)
  comparisonsAntigen <- c("S > N", "RBD > S", "RBD > N", "RBD > 0")
  antigenComps <- data.frame(Model="antigen",
                             Comparison=comparisonsAntigen,
                             Frequency=c(prob_SlargerN,
                                     prob_RBDlargerS,
                                     prob_RBDlargerN,
                                     prob_RBDlarger0))

  ### Technique
  techniqueWide <- slopePosteriors[["technique"]] %>%
    pivot_wider(., id_cols=.draw, names_from=parName, values_from=.value)
  #
  prob_IndirectlargerLFA <- with(techniqueWide, mean(LFA<0))
  prob_IndirectlargerDirect <- with(techniqueWide, mean(Direct<0))
  prob_IndirectlargerComp <- with(techniqueWide, mean(Neutralization<0))
  comparisonsTechnique <- c("Indirect > LFA", "Indirect > Direct",
    "Indirect > Competitive")
  techniqueComps <- data.frame(Model="technique",
                             Comparison=comparisonsTechnique,
                             Frequency=c(prob_IndirectlargerLFA,
                                         prob_IndirectlargerDirect,
                                         prob_IndirectlargerComp))

  ### Antibody
  antibodyWide <- slopePosteriors[["antibody"]] %>%
    pivot_wider(., id_cols=.draw, names_from=parName, values_from=.value)
  #
  prob_IGMlarger0 <- mean(antibodyWide$IgM>0)
  prob_IGAlarger0 <- mean(antibodyWide$IgA>0)
  prob_TOTlarger0 <- mean(antibodyWide$Total>0)
  prob_ALLlarger0 <- mean(with(antibodyWide, (IgM + IgA + Total)>0))
  comparisonsAntibody <- c("IgM > 0", "IgA > 0", "Total > 0",
                            "IgM + IgA + Total > 0")
  antibodyComps <- data.frame(Model="antibody",
                             Comparison=comparisonsAntibody,
                             Frequency=c(prob_IGMlarger0,
                                     prob_IGAlarger0,
                                     prob_TOTlarger0,
                                     prob_ALLlarger0))
}

### Save the comparisons performed
# COMMENT DEPENDING ON WHAT FIT TO USE, ALL DATA OR KNOWN TIMES
if (!knownTimesOnly) {
  allComparisons <- rbind(antigenComps, techniqueComps, antibodyComps, fullModelComps)
  write.csv(allComparisons,
            "../../data/analysis_results/05_characteristics_statistical_significance.csv",
            row.names=FALSE)
} else {
  write.csv(allComparisons,
            "../../data/analysis_results/08_characteristics_statistical_significance_known_times.csv",
            row.names=FALSE)
}


