#
# Plot sensitivity across time, both for the 'averages'
# of the different types of assays, and the assay-specific
# sensitivity profiles for associated paper
# "Dynamics of SARS-CoV-2 seroassay sensitivity: a systematic review and modeling study"
# Generates plots for Figure 2 (Main result of the paper), 
# and Figures S1, S3B, S4B, S6.
# 
# The inputs for plotting and statistical analysis are
# the posterior means and intervals of sensitivity
# across time, that are computed when fitting the
# models, the original data points (to overlay with the fits),
# and a dataset with original data and the results of
# cross-validation.
#
# Files with sensitivity profiles are in ../data/processed_data/
# and have names 'XXX_sensitivity_curve.csv', where XXX
# indicates the model fitted.
# 
# The original data points are in file
# ../data/processed_data/PCR_to_serotest_all.csv
#
# The results of crossvalidation are in
# '../data/analysis_results/XXX_grouped_CV', where
# XXX indicates the model that was used for crossvalidation
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
library(viridis)
library(RColorBrewer)
library(gridExtra)
library(ggthemes)
library(lemon)
source("../functions_auxiliary.R")
source("../functions_seroreversion_fit_analysis.R")

regLineSize=1
ribbonAlpha=0.2
locLineSize=0.4
locLineAlpha=0.4
locPointSize=2
sampleLineSize=0.2
sampleLineAlpha=0.1
dataAlpha <- 0.5

knownTimesOnly <- FALSE

# Load raw data
seroFitted <- read.csv("../../data/processed_data/PCR_to_serotest_all.csv",
                    stringsAsFactors=FALSE)

# Put characteristics columns to raw data
seroFitted$RBD <- stringr::str_detect(seroFitted$antigenTarget, "RBD")
seroFitted$S <- stringr::str_detect(seroFitted$antigenTarget, "S") &
  !seroFitted$RBD
seroFitted$N <- stringr::str_detect(seroFitted$antigenTarget, "N") &
  !(seroFitted$RBD | seroFitted$S)
seroFitted$LFA <- stringr::str_detect(seroFitted$technique, "LFIA")
seroFitted$Indirect <- ! seroFitted$LFA &
  stringr::str_detect(seroFitted$design, "indirect")
seroFitted$Direct <- ! seroFitted$LFA &
  stringr::str_detect(seroFitted$design, "sandwich")
seroFitted$Neutralization <- ! seroFitted$LFA &
  stringr::str_detect(seroFitted$design, "competitive")
binaryString <- c("No", "Yes")

if (knownTimesOnly) {
  seroFitted <- dplyr::filter(seroFitted, timeKnown)
}

####################################
####################################
#
# Full Model average sensitivity
# profiles (Fig 2) (MAIN PAPER FIGURE).
# Makes Fig S6 if knownTimesOnly = TRUE
#
####################################
####################################


if (!knownTimesOnly) {
  fullModelProfile <- read.csv("../../data/analysis_results/05_characteristics_fullModel_assay_sensitivity_curve.csv",
                             stringsAsFactors=FALSE) %>%
    dplyr::mutate(., averageSens=is.na(testName))
} else {
  fullModelProfile <- read.csv("../../data/analysis_results/08_characteristics_fullModel_assay_sensitivity_curve_known_times.csv",
                             stringsAsFactors=FALSE) %>%
    dplyr::mutate(., averageSens=is.na(testName))
}

# additional column to organize plots
fullModelProfile$Antigen[with(fullModelProfile, N)] <- "Nucleocapsid"
fullModelProfile$Antigen[with(fullModelProfile, S)] <- "Spike"
fullModelProfile$Antigen[with(fullModelProfile, RBD)] <- "Receptor-binding domain"

levelOrder <- c("Nucleocapsid", "Spike", "Receptor-binding domain")
fullModelProfile$Antigen <- factor(fullModelProfile$Antigen, levels=levelOrder)

binaryString <- c("No", "Yes")
fullModelProfile$LFA <- binaryString[as.integer(fullModelProfile$LFA)+1]
fullModelProfile$Direct <- binaryString[as.integer(fullModelProfile$Direct)+1]
fullModelProfile$Neutralization <- binaryString[as.integer(fullModelProfile$Neutralization)+1]
fullModelProfile$Indirect <- with(fullModelProfile, LFA=="No" & Direct=="No" & Neutralization=="No")
fullModelProfile$Indirect <- binaryString[as.integer(fullModelProfile$Indirect)+1]

fullModelProfile$Design[with(fullModelProfile, LFA=="Yes")] <- "LFA"
fullModelProfile$Design[with(fullModelProfile, LFA=="No" & Direct=="No" &
                             Neutralization=="No")] <- "Quantitative-Indirect"
fullModelProfile$Design[with(fullModelProfile, Direct=="Yes")] <- "Quantitative-Direct"
fullModelProfile$Design[with(fullModelProfile, Neutralization=="Yes")] <- "Quantitative-Competitive"

levelOrder2 <- c("LFA", "Quantitative-Indirect", "Quantitative-Direct", "Quantitative-Competitive")
fullModelProfile$Design <- factor(fullModelProfile$Design, levels=levelOrder2)

fullModelAverages <- filter(fullModelProfile, averageSens)
fullModelAssays <- filter(fullModelProfile, !averageSens)

testfullModel <- dplyr::select(fullModelAssays, testName, Antigen, LFA, Indirect,
                               Direct, Neutralization, Design) %>%
  unique()

fullModelPoints <- seroFitted %>%
  filter(., testName %in% fullModelProfile$testName) %>%
  select(., -LFA, -Indirect, -Direct, -Neutralization) %>%
  merge(., testfullModel) %>%
  dplyr::mutate(., time=testTime)

fullModelProfilePlot <- fullModelAverages %>%
  ggplot(., aes(x=time, y=sensitivityMean*100)) +
  geom_line(size=regLineSize, color="red") +
  geom_ribbon(aes(ymin=sensitivityL*100, ymax=sensitivityH*100),
              alpha=ribbonAlpha, color=NA, fill="red") +
  geom_point(data=fullModelPoints,
             aes(x=time, y=sensitivityMean, size=sqrt(nSamples)/10),
             color="black", alpha=0.2) +
  geom_line(data=fullModelAssays, aes(x=time, y=sensitivityMean*100, group=testName),
            color="black", alpha=0.7) +
  #facet_wrap(~Antigen+Design, labeller=label_both) +
  facet_grid(Antigen~Design, switch="y") +
  theme_bw() +
  #theme(legend.position="top", strip. background = element_blank()) +
  theme(legend.position="top", strip.placement = "outside") +
  guides(color=FALSE, size=FALSE, fill=FALSE) +
  xlim(0, 15) +
  xlab("Diagnosis to test (months)") +
  ylab("Sensitivity (%)")

if (knownTimesOnly) {
  figName <- "../../data/figures/characteristics_profiles_fullModel_known_times.png"
} else {
  figName <- "../../data/figures/characteristics_profiles_fullModel.pdf"
}
ggsave(figName, fullModelProfilePlot, units="cm", width=18, height=15)



#################################
#################################
#
# Individual assays sensitivity profiles (FIGURE S1)
#
#################################
#################################

if (!knownTimesOnly) {
  # Load validation results
  validationDfBasic <- read.csv("../../data/analysis_results/04_predicted_sensitivities_grouped_CV.csv",
                           stringsAsFactors=FALSE) %>%
  dplyr::mutate(., inInterval=(nSeropositives>=predictedPosL) &
                (nSeropositives<=predictedPosH))

  validationDf <- read.csv("../../data/analysis_results/05_characteristics_fullModel_grouped_CV.csv",
                           stringsAsFactors=FALSE) %>%
  dplyr::mutate(., inInterval=(nSeropositives>=predictedPosL) &
                (nSeropositives<=predictedPosH))


  # For tests not in the validation Df of the full model, use the
  # validation result of the basic model
  missingTestsFull <- validationDfBasic[!validationDfBasic$testName %in% validationDf$testName,]

  columns2remove <- names(validationDf)[!names(validationDf) %in% names(validationDfBasic)]
  validationDf <- select(validationDf, -columns2remove)
  validationDf <- rbind(validationDf, missingTestsFull)


  # Load sensitivity profile of basic model
  basicModelProfile <- read.csv("../../data/analysis_results/04_assay_sensitivity_curve.csv",
                                 stringsAsFactors=FALSE)
  fullModelProfile <- read.csv("../../data/analysis_results/05_characteristics_fullModel_assay_sensitivity_curve.csv",
                             stringsAsFactors=FALSE)

  nPlots <- 3 # Number of sheets into which to divide the assay-plots
  nTests <- length(unique(seroFitted$testName))
  testLists <- split(c(1:nTests), cut(seq_along(c(1:nTests)), nPlots,
                                      labels = FALSE)) 

  #####################
  # Assay sensitivity dynamics (Figure S1)
  #####################

  sensitivityCompPlot <- list()
  for (np in c(1:nPlots)) {
    testsPlot <- unique(seroFitted$testName)[testLists[[np]]]
    sensitivityCompPlot[[np]] <- dplyr::filter(validationDf, testName %in% testsPlot) %>%
      ggplot(aes(x=testTime, y=sensitivityMean, color=citationID,
                 shape=inInterval)) +
      geom_line(data=dplyr::filter(basicModelProfile,
                                   !is.na(testName) & (testName %in% testsPlot)),
                aes(x=time, y=sensitivityMean*100), color="black", linetype="solid",
                size=regLineSize, inherit.aes=FALSE) +
      geom_ribbon(data=dplyr::filter(basicModelProfile,
                                     !is.na(testName) & testName %in% testsPlot),
                  aes(x=time, ymin=sensitivityL*100, ymax=sensitivityH*100),
                  alpha=ribbonAlpha, colour=NA, show.legend=FALSE,
                  inherit.aes=FALSE) +
      geom_line(data=dplyr::filter(fullModelProfile,
                                   !is.na(testName) & (testName %in% testsPlot)),
                aes(x=time, y=sensitivityMean*100), color="red", linetype="solid",
                size=regLineSize, inherit.aes=FALSE) +
      geom_ribbon(data=dplyr::filter(fullModelProfile,
                                     !is.na(testName) & testName %in% testsPlot),
                  aes(x=time, ymin=sensitivityL*100, ymax=sensitivityH*100),
                  alpha=ribbonAlpha, colour=NA, show.legend=FALSE,
                  inherit.aes=FALSE, fill="red") +
      geom_pointrange(aes(ymin=sensitivityL, ymax=sensitivityH),
                      position=position_jitter(width=0.15, height=0)) +
      facet_wrap(.~testName, ncol=3, labeller = label_wrap_gen(width=35)) +
      theme_bw() +
      theme(legend.position="top") +
      guides(color=FALSE, shape=guide_legend("In validation interval")) +
      xlim(0, 15) +
      xlab("Diagnosis to test (months)") +
      ylab("Sensitivity (%)")
  ggsave(paste("../../data/figures/seroreversion_fit_characteristics", np, ".png", sep=""),
         sensitivityCompPlot[[np]], units="cm", width=20, height=28)
  }
}



#################################
#################################
#
# Average sensitivity profile plots for alternative models
# Figures S3 and S4
#
#################################
#################################

if (!knownTimesOnly) {
  ############
  # Antigen (Fig 2B) ## Comments in this block of code explain what happens in the
  # following blocks too
  ############

  # Load data, and add a column (averageSens) indicating if a given row of
  # the posterior sensitivity table corresponds to a specific test (indicated
  # in testName), or if it is the average of the test type (is.na(testName))
  antigenProfile <- read.csv("../../data/analysis_results/05_characteristics_antigen_assay_sensitivity_curve.csv",
                             stringsAsFactors=FALSE) %>%
    dplyr::mutate(., averageSens=is.na(testName))

  # Make columns with names for the different combinations of parameters (kinds of tests)
  antigenProfile$antigen[with(antigenProfile, N)] <- "N"
  antigenProfile$antigen[with(antigenProfile, S)] <- "S"
  antigenProfile$antigen[with(antigenProfile, RBD)] <- "RBD"

  # Separate the dataframe into assay-specific sensitivities and average sensitivities
  antigenAverages <- filter(antigenProfile, averageSens)
  antigenAssays <- filter(antigenProfile, !averageSens)

  # Add a column indicating the kind of test to the original data,
  # to plot the data points
  testAntigens <- dplyr::select(antigenAssays, testName, antigen) %>%
    unique()
  antigenPoints <- filter(seroFitted, testName %in% antigenProfile$testName) %>%
    merge(., testAntigens) %>%
    dplyr::mutate(., time=testTime)


  antigenProfilePlot <- antigenAverages %>%
    ggplot(., aes(x=time, y=sensitivityMean*100)) +
    geom_line(size=regLineSize, color="red") +
    geom_ribbon(aes(ymin=sensitivityL*100, ymax=sensitivityH*100),
                alpha=ribbonAlpha, color=NA, fill="red") +
    geom_point(data=antigenPoints,
               aes(x=time, y=sensitivityMean, size=sqrt(nSamples)/10),
               color="black", alpha=0.2) +
    geom_line(data=antigenAssays, aes(x=time, y=sensitivityMean*100, group=testName),
              color="black", alpha=0.7) +
    facet_wrap(.~antigen, ncol=3, labeller=label_wrap_gen(width=35)) +
    theme_bw() +
    theme(legend.position="top") +
    guides(color=FALSE, size=FALSE, fill=FALSE) +
    xlim(0, 15) +
    xlab("Diagnosis to test (months)") +
    ylab("Sensitivity (%)")

  ggsave("../data/figures/characteristics_profiles_antigen.png", antigenProfilePlot,
         units="cm", width=16, height=11)


  ############
  # Technique (Fig 3B)
  ############

  techniqueProfile <- read.csv("../../data/analysis_results/05_characteristics_technique_assay_sensitivity_curve.csv",
                             stringsAsFactors=FALSE) %>%
    dplyr::mutate(., averageSens=is.na(testName))

  # Make columns with names for the different combinations of parameters (kinds of tests)
  techniqueProfile$technique[techniqueProfile$LFA] <- "LFA"
  techniqueProfile$technique[techniqueProfile$Direct] <- "Quantitative-Direct"
  techniqueProfile$technique[techniqueProfile$Neutralization] <- "Quantitative-Competitive"
  techniqueProfile$technique[is.na(techniqueProfile$technique)] <- "Quantitative-Indirect"

  techniqueAverages <- filter(techniqueProfile, averageSens)
  techniqueAssays <- filter(techniqueProfile, !averageSens)

  testTechniques <- dplyr::select(techniqueAssays, testName, technique) %>%
    unique()

  techniquePoints <- filter(seroFitted, testName %in% techniqueProfile$testName) %>%
    dplyr::select(., -technique) %>%
    merge(., testTechniques) %>%
    dplyr::mutate(., time=testTime)

  techniqueProfilePlot <- techniqueAverages %>%
    ggplot(., aes(x=time, y=sensitivityMean*100)) +
    geom_line(size=regLineSize, color="red") +
    geom_ribbon(aes(ymin=sensitivityL*100, ymax=sensitivityH*100),
                alpha=ribbonAlpha, color=NA, fill="red") +
    geom_point(data=techniquePoints,
               aes(x=time, y=sensitivityMean, size=sqrt(nSamples)/10),
               color="black", alpha=0.2) +
    geom_line(data=techniqueAssays, aes(x=time, y=sensitivityMean*100, group=testName),
              color="black", alpha=0.7) +
    facet_wrap(.~technique, ncol=2, labeller=label_wrap_gen(width=35)) +
    theme_bw() +
    theme(legend.position="top") +
    guides(color=FALSE, size=FALSE, fill=FALSE) +
    xlim(0, 15) +
    xlab("Diagnosis to test (months)") +
    ylab("Sensitivity (%)")

  ggsave("../../data/figures/characteristics_profiles_technique.png", techniqueProfilePlot,
         units="cm", width=10, height=12)

  ############
  # Antibody (Fig 5B)
  ############

  antibodyProfile <- read.csv("../../data/analysis_results/05_characteristics_antibody_assay_sensitivity_curve.csv",
                             stringsAsFactors=FALSE) %>%
    dplyr::mutate(., averageSens=is.na(testName))

  # Make columns with names for the different combinations of parameters (kinds of tests)
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
    ggplot(., aes(x=time, y=sensitivityMean*100)) +
    geom_line(size=regLineSize, color="red") +
    geom_ribbon(aes(ymin=sensitivityL*100, ymax=sensitivityH*100),
                alpha=ribbonAlpha, color=NA, fill="red") +
    geom_point(data=antibodyPoints,
               aes(x=time, y=sensitivityMean, size=sqrt(nSamples)/10),
               color="black", alpha=0.2) +
    geom_line(data=antibodyAssays, aes(x=time, y=sensitivityMean*100, group=testName),
              color="black", alpha=0.7) +
    facet_wrap(.~antibody, ncol=2, labeller=label_wrap_gen(width=35)) +
    theme_bw() +
    theme(legend.position="top") +
    guides(color=FALSE, size=FALSE, fill=FALSE) +
    xlim(0, 15) +
    xlab("Diagnosis to test (months)") +
    ylab("Sensitivity (%)")

  ggsave("../../data/figures/characteristics_profiles_antibody.png", antibodyProfilePlot,
         units="cm", width=12, height=12)
}

