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


modelNames <- c("antigen", "antibody", "technique", "fullModel")

charPlot <- list()
slopePosteriors <- list()
for (m in c(1:length(modelNames))) {
  mN <- modelNames[m]
  posteriorsName <- paste("../data/analysis_results/05_characteristics_",
                          mN, "_posterior_samples.csv", sep="")
  slopePosteriors[[mN]] <- read.csv(posteriorsName, stringsAsFactors=FALSE) %>%
    dplyr::filter(., !is.na(parName))
  charPlot[[m]] <- slopePosteriors[[mN]] %>%
    ggplot(data=., aes(x=.value, fill=parName)) +
    geom_density(alpha=0.3) +
    theme_bw() +
    theme(axis.title.y=element_blank()) +
    guides(fill=guide_legend("Parameter")) +
    xlab("Value")
  plotName <- paste("../data/figures/characteristics_effects_plot_",
                    mN, ".png", sep="")
  ggsave(plotName, charPlot[[m]], units="cm", width=16, height=12)
}

comparisonsDf <- NULL

# Get some statistics comparing parameters

### Antigen
antigenWide <- slopePosteriors[["antigen"]] %>%
  pivot_wider(., id_cols=.draw, names_from=parName, values_from=.value)

prob_SlargerN <- mean(antigenWide$S>antigenWide$N)
prob_RBDlarger0 <- mean((antigenWide$S+antigenWide$RBD)>0)
prob_SNlarger0 <- mean(antigenWide$SN>0)
comparisonsAntigen <- c("S > N", "S+RBD > 0", "SN > 0")

antigenComps <- data.frame(Model="antigen",
                           Comparison=comparisonsAntigen,
                           Frequency=c(prob_SlargerN,
                                   prob_RBDlarger0,
                                   prob_SNlarger0))

### Technique
techniqueWide <- slopePosteriors[["technique"]] %>%
  pivot_wider(., id_cols=.draw, names_from=parName, values_from=.value)

prob_RESTlargerLFA <- mean(techniqueWide$LFA<techniqueWide$Rest)
comparisonsTechnique <- c("Rest > LFA")

techniqueComps <- data.frame(Model="technique",
                           Comparison=comparisonsTechnique,
                           Frequency=c(prob_RESTlargerLFA))



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
prob_SlargerN <- with(fullModelWide, mean(S>N))
prob_RBDlarger0 <- mean(fullModelWide$RBD>0)
prob_SRBDlarger0 <- with(fullModelWide, mean(RBD+S>0))

comparisonsFullModel <- c("LFA > 0", "S > N", "RBD > 0",
                          "S + RBD > 0")

fullModelComps <- data.frame(Model="fullModel",
                           Comparison=comparisonsFullModel,
                           Frequency=c(prob_LFAlarger0,
                                       prob_SlargerN,
                                       prob_RBDlarger0,
                                       prob_SRBDlarger0))

### Save the comparisons performed
allComparisons <- rbind(antigenComps, techniqueComps, antibodyComps,
                        fullModelComps)

write.csv(allComparisons,
          "../data/analysis_results/05_characteristics_statistical_significance.csv",
          row.names=FALSE)

