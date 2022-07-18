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

# Load raw data
seroFitted <- read.csv("../data/processed_data/PCR_to_serotest_all.csv",
                    stringsAsFactors=FALSE)

# Load validation results
validationDf <- read.csv("../data/analysis_results/04_predicted_sensitivities_grouped_CV.csv",
                         stringsAsFactors=FALSE) %>%
dplyr::mutate(., inInterval=(nSeropositives>=predictedPosL) &
              (nSeropositives<=predictedPosH))


#################################
#################################
# Individual assays sensitivity profiles
#################################
#################################

#####################
# Basic model
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

#####################
# Full model
#####################
fullModelProfile <- read.csv("../data/analysis_results/05_characteristics_fullModel_assay_sensitivity_curve.csv",
                           stringsAsFactors=FALSE)

sensitivityCompPlot <- list()
for (np in c(1:nPlots)) {
  testsPlot <- unique(seroFitted$testName)[testLists[[np]]]
  sensitivityCompPlot[[np]] <- sensitivityPlot[[np]] +
    geom_line(data=dplyr::filter(fullModelProfile,
                                 !is.na(testName) & (testName %in% testsPlot)),
              aes(x=time, y=sensitivityMean*100), color="red", linetype="solid",
              size=regLineSize, inherit.aes=FALSE) +
    geom_ribbon(data=dplyr::filter(fullModelProfile,
                                   !is.na(testName) & testName %in% testsPlot),
                aes(x=time, ymin=sensitivityL*100, ymax=sensitivityH*100),
                alpha=ribbonAlpha, colour=NA, show.legend=FALSE,
                inherit.aes=FALSE, fill="red")
ggsave(paste("../data/figures/seroreversion_fit_characteristics", np, ".png", sep=""),
       sensitivityCompPlot[[np]], units="cm", width=24, height=30)
}



#################################
#################################
# Average assay sensitivity profiles
#################################
#################################

############
# Antigen
############

antigenProfile <- read.csv("../data/analysis_results/05_characteristics_antigen_assay_sensitivity_curve.csv",
                           stringsAsFactors=FALSE) %>%
  dplyr::mutate(., averageSens=is.na(testName))

antigenProfile$antigen[with(antigenProfile, N & (!S))] <- "N"
antigenProfile$antigen[with(antigenProfile, N & S & !(RBD))] <- "N-S"
antigenProfile$antigen[with(antigenProfile, !(N) & S & !(RBD))] <- "S"
antigenProfile$antigen[with(antigenProfile, !(N) & S & RBD)] <- "S-RBD"
antigenProfile$antigen[with(antigenProfile, N & S & RBD)] <- "N-S-RBD"

antigenAverages <- filter(antigenProfile, averageSens)
antigenAssays <- filter(antigenProfile, !averageSens)

testAntigens <- dplyr::select(antigenAssays, testName, antigen) %>%
  unique()

antigenPoints <- filter(seroFitted, testName %in% antigenProfile$testName) %>%
  merge(., testAntigens) %>%
  dplyr::mutate(., time=testTime)


antigenProfilePlot <- antigenAverages %>%
  ggplot(., aes(x=time, y=sensitivityMean*100, color=antigen, fill=antigen)) +
  geom_point(data=antigenPoints,
             aes(x=time, y=sensitivityMean, size=sqrt(nSamples)/10),
             color="black", alpha=0.3) +
  geom_line(data=antigenAssays, aes(x=time, y=sensitivityMean*100, group=testName),
            color="black", alpha=0.3) +
  geom_line(size=regLineSize) +
  geom_ribbon(aes(ymin=sensitivityL*100, ymax=sensitivityH*100),
              alpha=ribbonAlpha) +
  facet_wrap(.~antigen, ncol=3, labeller=label_wrap_gen(width=35)) +
  theme_bw() +
  theme(legend.position="top") +
  guides(color=FALSE, size=FALSE, fill=FALSE) +
  xlim(0, 15) +
  xlab("Diagnosis to test (months)") +
  ylab("Sensitivity (%)")

ggsave("../data/figures/characteristics_profiles_antigen.png", antigenProfilePlot,
       units="cm", width=24, height=15)



############
# Technique
############

techniqueProfile <- read.csv("../data/analysis_results/05_characteristics_technique_assay_sensitivity_curve.csv",
                           stringsAsFactors=FALSE) %>%
  dplyr::mutate(., averageSens=is.na(testName))

techniqueProfile$technique[techniqueProfile$LFA] <- "LFA"
techniqueProfile$technique[!techniqueProfile$LFA] <- "Rest"

techniqueAverages <- filter(techniqueProfile, averageSens)
techniqueAssays <- filter(techniqueProfile, !averageSens)

testTechniques <- dplyr::select(techniqueAssays, testName, technique) %>%
  unique()

techniquePoints <- filter(seroFitted, testName %in% techniqueProfile$testName) %>%
  dplyr::select(., -technique) %>%
  merge(., testTechniques) %>%
  dplyr::mutate(., time=testTime)


techniqueProfilePlot <- techniqueAverages %>%
  ggplot(., aes(x=time, y=sensitivityMean*100, color=technique, fill=technique)) +
  geom_point(data=techniquePoints,
             aes(x=time, y=sensitivityMean, size=sqrt(nSamples)/10),
             color="black", alpha=0.3) +
  geom_line(data=techniqueAssays, aes(x=time, y=sensitivityMean*100, group=testName),
            color="black", alpha=0.3) +
  geom_line(size=regLineSize) +
  geom_ribbon(aes(ymin=sensitivityL*100, ymax=sensitivityH*100),
              alpha=ribbonAlpha) +
  facet_wrap(.~technique, ncol=2, labeller=label_wrap_gen(width=35)) +
  theme_bw() +
  theme(legend.position="top") +
  guides(color=FALSE, size=FALSE, fill=FALSE) +
  xlim(0, 15) +
  xlab("Diagnosis to test (months)") +
  ylab("Sensitivity (%)")

ggsave("../data/figures/characteristics_profiles_technique.png", techniqueProfilePlot,
       units="cm", width=19, height=15)



############
# Antibody
############

antibodyProfile <- read.csv("../data/analysis_results/05_characteristics_antibody_assay_sensitivity_curve.csv",
                           stringsAsFactors=FALSE) %>%
  dplyr::mutate(., averageSens=is.na(testName))

antibodyProfile$antibody[with(antibodyProfile, (!IgM) & (!IgA))] <- "IgG"
antibodyProfile$antibody[with(antibodyProfile, IgM & (!IgA))] <- "IgG-IgM"
antibodyProfile$antibody[with(antibodyProfile, (!IgM) & IgA)] <- "IgG-IgA"
antibodyProfile$antibody[with(antibodyProfile, IgM & IgA)] <- "IgG-IgM-IgA"

antibodyAverages <- filter(antibodyProfile, averageSens)
antibodyAssays <- filter(antibodyProfile, !averageSens)

testAntibodies <- dplyr::select(antibodyAssays, testName, antibody) %>%
  unique()

antibodyPoints <- filter(seroFitted, testName %in% antibodyProfile$testName) %>%
  merge(., testAntibodies) %>%
  dplyr::mutate(., time=testTime)

antibodyProfilePlot <- antibodyAverages %>%
  ggplot(., aes(x=time, y=sensitivityMean*100, color=antibody, fill=antibody)) +
  geom_point(data=antibodyPoints,
             aes(x=time, y=sensitivityMean, size=sqrt(nSamples)/10),
             color="black", alpha=0.3) +
  geom_line(data=antibodyAssays, aes(x=time, y=sensitivityMean*100, group=testName),
            color="black", alpha=0.3) +
  geom_line(size=regLineSize) +
  geom_ribbon(aes(ymin=sensitivityL*100, ymax=sensitivityH*100),
              alpha=ribbonAlpha) +
  facet_wrap(.~antibody, ncol=2, labeller=label_wrap_gen(width=35)) +
  theme_bw() +
  theme(legend.position="top") +
  guides(color=FALSE, size=FALSE, fill=FALSE) +
  xlim(0, 15) +
  xlab("Diagnosis to test (months)") +
  ylab("Sensitivity (%)")

ggsave("../data/figures/characteristics_profiles_antibody.png", antibodyProfilePlot,
       units="cm", width=19, height=15)


############
# Full Model
############



