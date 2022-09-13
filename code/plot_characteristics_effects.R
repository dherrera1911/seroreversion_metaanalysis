####################################################
#
# Plot the density plots of the model parameters related to
# assay characteristics. Generates plots for figures
# 2A, 3A, 4A, 5A, 6A.
# 
# Also, computes some comparisons of the posterior samples
# that are used for significance reporting in the text.
# E.g. the proportion of N < RBD + S.
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

# Some plotting parameters
regLineSize=1
ribbonAlpha=0.2
locLineSize=0.4
locLineAlpha=0.4
locPointSize=1.6
sampleLineSize=0.1
sampleLineAlpha=0.1
dataAlpha <- 0.5

# Vector with model identifiers
modelNames <- c("antigen", "antibody", "technique", "design", "fullModel")

########
### Loop below loads each models posteriors, and makes/saves the density plots.
### Makes figures 2A, 3A, 4A, 5A, 6A.
########

charPlot <- list()
slopePosteriors <- list()
for (m in c(1:length(modelNames))) {
  mN <- modelNames[m]
  posteriorsName <- paste("../data/analysis_results/05_characteristics_",
                          mN, "_posterior_samples.csv", sep="")
  slopePosteriors[[mN]] <- read.csv(posteriorsName, stringsAsFactors=FALSE) %>%
    dplyr::filter(., !is.na(parName))
  if (mN == "antigen") {
    slopePosteriors[[mN]]$parName <- factor(slopePosteriors[[mN]]$parName,
                                            levels=c("N", "S", "RBD"))
  }
  if (mN == "fullModel") {
    slopePosteriors[[mN]]$isAntigen <- slopePosteriors[[mN]]$parName %in%
      c("S", "N", "RBD")
    levelOrder <- c("N", "S", "RBD", "LFA", "sandwich")
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
  plotName <- paste("../data/figures/characteristics_effects_plot_",
                    mN, ".png", sep="")
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

### Antigen
antigenWide <- slopePosteriors[["antigen"]] %>%
  pivot_wider(., id_cols=.draw, names_from=parName, values_from=.value)

prob_SlargerN <- mean(antigenWide$S>antigenWide$N)
prob_RBDlargerS <- mean(antigenWide$RBD>antigenWide$S)
prob_RBDlargerN <- mean(antigenWide$RBD>antigenWide$N)
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

prob_RESTlargerLFA <- mean(techniqueWide$LFA<techniqueWide$Rest)
comparisonsTechnique <- c("Rest > LFA")

techniqueComps <- data.frame(Model="technique",
                           Comparison=comparisonsTechnique,
                           Frequency=c(prob_RESTlargerLFA))

### Design
designWide <- slopePosteriors[["design"]] %>%
  pivot_wider(., id_cols=.draw, names_from=parName, values_from=.value)

prob_SandwichlargerNotSandwich <- mean(designWide$notSandwich<designWide$sandwich)
prob_Sandwichlarger0 <- mean(0<designWide$sandwich)
prob_NotSandwichlarger0 <- mean(0<designWide$notSandwich)
comparisonsTechnique <- c("Sandwich > Not Sandwich", "Sandwich > 0",
  "Not Sandwich > 0")

designComps <- data.frame(Model="design",
                           Comparison=comparisonsTechnique,
                           Frequency=c(prob_SandwichlargerNotSandwich,
                                       prob_Sandwichlarger0,
                                       prob_NotSandwichlarger0))

### Antibody
antibodyWide <- slopePosteriors[["antibody"]] %>%
  pivot_wider(., id_cols=.draw, names_from=parName, values_from=.value)

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

### Full model
fullModelWide <- slopePosteriors[["fullModel"]] %>%
  pivot_wider(., id_cols=.draw, names_from=parName, values_from=.value)

prob_LFAlarger0 <- mean(fullModelWide$LFA>0)
prob_Sandwichlarger0 <- mean(fullModelWide$sandwich>0)
prob_SlargerN <- with(fullModelWide, mean(S>N))
prob_RBDlarger0 <- mean(fullModelWide$RBD>0)
prob_RBDlargerS <- with(fullModelWide, mean(RBD>S))
prob_RBDlargerN <- with(fullModelWide, mean(RBD>N))

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

### Save the comparisons performed
allComparisons <- rbind(antigenComps, techniqueComps, designComps,
                        antibodyComps, fullModelComps)

write.csv(allComparisons,
          "../data/analysis_results/05_characteristics_statistical_significance.csv",
          row.names=FALSE)

