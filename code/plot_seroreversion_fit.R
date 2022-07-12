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

#######################
#######################
## 1) Plot the fit to all data overlaid on the raw data
#######################
#######################

### Load data
# Load the fitted serology data
seroFitted <- read.csv("../data/processed_data/PCR_to_serotest_all.csv",
                    stringsAsFactors=FALSE)

# Get the samples of the sensitivity across time posterior for each test
testNames <- unique(seroFitted$testName)
timeVec <- seq(0.5, 12, 0.2)
assayFitDf <- NULL
assaySummaryDf <- NULL
for (tN in testNames) {
  # get posterior 
  testTrace <- dplyr::filter(posteriorTraces, assay==tN)
  sensitivityPosterior <- seroreversion_samples(testTrace, timeVec,
                                                slopeName="assaySlope",
                                                interceptName="assayIntercept",
                                                timeNormalization=4)
  testDf <- data.frame(time=timeVec,
                       sensitivity=sensitivityPosterior$prop_mean,
                       sensitivityL=sensitivityPosterior$prop_L,
                       sensitivityH=sensitivityPosterior$prop_H,
                       testName=tN)
  assayFitDf <- rbind(assayFitDf, testDf)
  summaryStats <- group_by(testTrace, .variable) %>%
    summarize(., meanVal=mean(.value), sdVal=sd(.value))
  assayIntercept <- summaryStats$meanVal[summaryStats$.variable=="assayIntercept"]
  assaySlope <- summaryStats$meanVal[summaryStats$.variable=="assaySlope"]
  assayInterceptSD <- summaryStats$sdVal[summaryStats$.variable=="assayIntercept"]
  assaySlopeSD <- summaryStats$sdVal[summaryStats$.variable=="assaySlope"]
  summaryDf <- data.frame(intercept=assayIntercept,
                          slope=assaySlope,
                          interceptSD=assayInterceptSD,
                          slopeSD=assaySlopeSD, testName=tN)
  assaySummaryDf <- rbind(assaySummaryDf, summaryDf)
}

write.csv(assaySummaryDf, "../data/analysis_results/assay_fit_summary.csv",
          row.names=FALSE)
          

# Plot the fitting results with the raw data
nPlots <- 3
nTests <- length(unique(seroFitted$testName))
testLists <- split(c(1:nTests), cut(seq_along(c(1:nTests)), nPlots,
                                    labels = FALSE)) 
sensitivityPlot <- list()
for (np in c(1:nPlots)) {
  testsPlot <- unique(seroFitted$testName)[testLists[[np]]]
  sensitivityPlot[[np]] <- dplyr::filter(seroFitted, testName %in% testsPlot) %>%
    ggplot(aes(x=testTime, y=sensitivityMean, color=citationID,
               shape=sampleType)) +
    geom_pointrange(aes(ymin=sensitivityL, ymax=sensitivityH),
                    position=position_jitter(width=0.15, height=0)) +
    geom_line(data=dplyr::filter(assayFitDf, testName %in% testsPlot),
              aes(x=time, y=sensitivity*100), color="black", linetype="solid",
              size=regLineSize, inherit.aes=FALSE) +
    geom_ribbon(data=dplyr::filter(assayFitDf, testName %in% testsPlot),
                aes(x=time, ymin=sensitivityL*100, ymax=sensitivityH*100),
                alpha=ribbonAlpha, colour=NA, show.legend=FALSE,
                inherit.aes=FALSE) +
    facet_wrap(.~testName, ncol=3, labeller = label_wrap_gen(width=35)) +
    theme_bw() +
    theme(legend.position="top") +
    guides(color=FALSE, shape=guide_legend("Sample type")) +
    xlim(0, 15) +
    xlab("Diagnosis to test (months)") +
    ylab("Sensitivity (%)")
ggsave(paste("../data/figures/seroreversion_fit", np, ".png", sep=""),
       sensitivityPlot[[np]], units="cm", width=24, height=30)
}


#######################
#######################
## 2) Plot some summaries of parameters
#######################
#######################

# Plot model parameters vs their priors
sigmas <- c("interceptSigma", "slopeSigma", "studySigma")
modelSigmas <- dplyr::filter(posteriorTraces, .variable%in%sigmas)
gammaPars <- list(c(4,4), c(4,4), c(4,4))
priorDf <- NULL
for (ss in c(1:length(sigmas))) {
  x <- seq(0, 2.5, 0.01)
  y <- dgamma(x, gammaPars[[ss]][1], gammaPars[[ss]][2])
  priorDf <- rbind(priorDf, data.frame(x=x, y=y, .variable=sigmas[ss]))
}
sigmasPlot <- ggplot(data=modelSigmas, aes(x=.value)) +
  geom_density() +
  facet_wrap(.~.variable) +
  geom_line(data=priorDf, aes(x=x, y=y), color="red")

# Sigmas credible intervals
meanSlopeSigma <- mean(dplyr::filter(modelSigmas, .variable=="slopeSigma")[[".value"]])
ciSlopeSigma <- quantile(dplyr::filter(modelSigmas, .variable=="slopeSigma")[[".value"]],
probs=c(0.025, 0.975))
meanInterceptSigma <- mean(dplyr::filter(modelSigmas, .variable=="interceptSigma")[[".value"]])
ciInterceptSigma <- quantile(dplyr::filter(modelSigmas, .variable=="interceptSigma")[[".value"]],
probs=c(0.025, 0.975))
meanStudySigma <- mean(dplyr::filter(modelSigmas, .variable=="studySigma")[[".value"]])
ciStudySigma <- quantile(dplyr::filter(modelSigmas, .variable=="studySigma")[[".value"]],
probs=c(0.025, 0.975))


# Plot parameters
slope_intercept <- ggplot(data=assaySummaryDf, aes(x=intercept, y=slope)) +
  geom_point() +
  theme_bw()
corTest <- cor.test(assaySummaryDf$slope, assaySummaryDf$intercept)
ggsave("../data/figures/seroreversion_assays_correlation.png", slope_intercept,
  units="cm", width=20, height=20)

# Plot histogram of parameters
# Make short labels
assaySummaryDf$testName2 <-
  stringr::str_replace(assaySummaryDf$testName, " \\s*\\([^\\)]+\\)", "") %>%
  stringr::str_replace(., "Anti-SARS-CoV-2", "") %>%
  stringr::str_replace(., "SARS-CoV-2", "") %>%
  stringr::str_trunc(., width=25, ellipsis="")

sortedAssays <- assaySummaryDf$testName2[order(assaySummaryDf$slope, decreasing=T)]
assaySummaryDf$testName2 <- factor(assaySummaryDf$testName2, levels=sortedAssays)
slopeHist <- assaySummaryDf %>%
  ggplot(., aes(x=slope, y=testName2)) +
  geom_pointrange(aes(xmin=slope-slopeSD*1.96, xmax=slope+slopeSD*1.96)) +
  theme_bw() +
  theme(axis.title.y=element_blank()) +
  xlab("Assay slope")

sortedAssays <- assaySummaryDf$testName2[order(assaySummaryDf$intercept, decreasing=T)]
assaySummaryDf$testName2 <- factor(assaySummaryDf$testName2, levels=sortedAssays)
interceptHist <- assaySummaryDf %>%
  ggplot(., aes(x=intercept, y=testName2)) +
  geom_pointrange(aes(xmin=intercept-interceptSD*1.96,
                      xmax=intercept+interceptSD*1.96)) +
  theme_bw() +
  theme(axis.title.y=element_blank()) +
  xlab("Assay intercept")

paramPlot <- ggarrange(plotlist=list(slopeHist, interceptHist),
                       ncol=2)

ggsave("../data/figures/assayPars_general.png", paramPlot,
  units="cm", width=22, height=14)



#######################
#######################
## 3) Plot parameters by assay characteristics
#######################
#######################
#assaySummaryDf <- read.csv("../data/analysis_results/assay_fit_summary.csv",
#                           stringsAsFactors=FALSE)
#
## Load assay characteristics dataframe
#assayChars <- read.csv("../data/raw_data/assay_characteristics.csv",
#                       stringsAsFactors=FALSE, na.strings=c("", "--")) %>%
#  dplyr::rename(., testName=test_name_long, antigen_target=antigen.target,
#  assay_type=assay.type, isotype_target=isotype.target)
#
## Merge assay characteristics with assay fit
#assaySummaryDf <- dplyr::filter(assayChars, testName %in% assaySummaryDf$testName) %>%
#  merge(assaySummaryDf, ., all.x=TRUE)
#assaySummaryDf$antigen_target[is.na(assaySummaryDf$antigen_target)] <- "NA"
#assaySummaryDf$assay_type[is.na(assaySummaryDf$assay_type)] <- "NA"
#assaySummaryDf$isotype_target[is.na(assaySummaryDf$isotype_target)] <- "NA"
#
# Summarize the slope distribution grouped by assay characteristics
#antigenSummary <- group_by(assaySummaryDf, antigen_target) %>%
#  summarize(., meanSlope=mean(slope), sdSlope=sd(slope))
#antigenSlopes <- assaySummaryDf %>%
#  ggplot(., aes(x=antigen_target, y=slope, color=antigen_target)) +
#  geom_pointrange(aes(ymin=slope-slopeSD*1.96,
#                      ymax=slope+slopeSD*1.96),
#                  position=position_jitter(width=0.4, height=0),
#                  alpha=0.3) +
#  geom_pointrange(data=antigenSummary,
#                  aes(y=meanSlope, ymin=meanSlope-sdSlope*1.96,
#                      ymax=meanSlope+sdSlope*1.96), size=1) + 
#  theme_bw() +
#  theme(legend.position="none") +
#  xlab("Antigen") +
#  ylab("Slope")
#
#typeSummary <- group_by(assaySummaryDf, assay_type) %>%
#  summarize(., meanSlope=mean(slope), sdSlope=sd(slope))
#typeSlopes <- assaySummaryDf %>%
#  ggplot(., aes(x=assay_type, y=slope, color=assay_type)) +
#  geom_pointrange(aes(ymin=slope-slopeSD*1.96,
#                      ymax=slope+slopeSD*1.96),
#                  position=position_jitter(width=0.4, height=0),
#                  alpha=0.3) +
#  geom_pointrange(data=typeSummary,
#                  aes(y=meanSlope, ymin=meanSlope-sdSlope*1.96,
#                      ymax=meanSlope+sdSlope*1.96), size=1) + 
#  theme_bw() +
#  theme(legend.position="none") +
#  xlab("Type") +
#  ylab("Slope")
#
#isotypeSummary <- group_by(assaySummaryDf, isotype_target) %>%
#  summarize(., meanSlope=mean(slope), sdSlope=sd(slope))
#isotypeSlopes <- assaySummaryDf %>%
#  ggplot(., aes(x=isotype_target, y=slope, color=isotype_target)) +
#  geom_pointrange(aes(ymin=slope-slopeSD*1.96,
#                      ymax=slope+slopeSD*1.96),
#                  position=position_jitter(width=0.4, height=0),
#                  alpha=0.3) +
#  geom_pointrange(data=isotypeSummary,
#                  aes(y=meanSlope, ymin=meanSlope-sdSlope*1.96,
#                      ymax=meanSlope+sdSlope*1.96), size=1) + 
#  theme_bw() +
#  theme(legend.position="none") +
#  xlab("Isotype") +
#  ylab("Slope")
#
#
#antigenSummary <- group_by(assaySummaryDf, antigen_target) %>%
#  summarize(., meanSlope=mean(intercept), sdSlope=sd(intercept))
#antigenIntercepts <- assaySummaryDf %>%
#  ggplot(., aes(x=antigen_target, y=intercept, color=antigen_target)) +
#  geom_pointrange(aes(ymin=intercept-interceptSD*1.96,
#                      ymax=intercept+interceptSD*1.96),
#                  position=position_jitter(width=0.4, height=0),
#                  alpha=0.3) +
#  geom_pointrange(data=antigenSummary,
#                  aes(y=meanSlope, ymin=meanSlope-sdSlope*1.96,
#                      ymax=meanSlope+sdSlope*1.96), size=1) + 
#  theme_bw() +
#  theme(legend.position="none") +
#  xlab("Antigen") +
#  ylab("Intercept")
#
#typeSummary <- group_by(assaySummaryDf, assay_type) %>%
#  summarize(., meanSlope=mean(intercept), sdSlope=sd(intercept))
#typeIntercepts <- assaySummaryDf %>%
#  ggplot(., aes(x=assay_type, y=intercept, color=assay_type)) +
#  geom_pointrange(aes(ymin=intercept-interceptSD*1.96,
#                      ymax=intercept+interceptSD*1.96),
#                  position=position_jitter(width=0.4, height=0),
#                  alpha=0.3) +
#  geom_pointrange(data=typeSummary,
#                  aes(y=meanSlope, ymin=meanSlope-sdSlope*1.96,
#                      ymax=meanSlope+sdSlope*1.96), size=1) + 
#  theme_bw() +
#  theme(legend.position="none") +
#  xlab("Type") +
#  ylab("Intercept")
#
#isotypeSummary <- group_by(assaySummaryDf, isotype_target) %>%
#  summarize(., meanIntercept=mean(intercept), sdIntercept=sd(intercept))
#isotypeIntercepts <- assaySummaryDf %>%
#  ggplot(., aes(x=isotype_target, y=intercept, color=isotype_target)) +
#  geom_pointrange(aes(ymin=intercept-interceptSD*1.96,
#                      ymax=intercept+interceptSD*1.96),
#                  position=position_jitter(width=0.4, height=0),
#                  alpha=0.3) +
#  geom_pointrange(data=isotypeSummary,
#                  aes(y=meanIntercept, ymin=meanIntercept-sdIntercept*1.96,
#                      ymax=meanIntercept+sdIntercept*1.96), size=1) + 
#  theme_bw() +
#  theme(legend.position="none") +
#  xlab("Isotype") +
#  ylab("Intercept")
#
#
#ggsave("../data/figures/assayPars_slopes_by_antigen.png", antigenSlopes,   
#       units="cm", width=20, height=7)
#ggsave("../data/figures/assayPars_slopes_by_type.png", typeSlopes,   
#       units="cm", width=20, height=7)
#ggsave("../data/figures/assayPars_slopes_by_isotype.png", isotypeSlopes,   
#       units="cm", width=20, height=7)
#ggsave("../data/figures/assayPars_intercepts_by_antigen.png", antigenIntercepts,   
#       units="cm", width=20, height=7)
#ggsave("../data/figures/assayPars_intercepts_by_type.png", typeIntercepts,   
#       units="cm", width=20, height=7)
#ggsave("../data/figures/assayPars_intercepts_by_isotype.png", isotypeIntercepts,   
#       units="cm", width=20, height=7)
#

#######################
#######################
## 4) Plot fitted vs predicted values
#######################
#######################

validationDf <- read.csv("../data/analysis_results/sensitivity_decay_validation.csv",
                         stringsAsFactors=FALSE) %>%
  dplyr::mutate(., sensitivityMeanPred=sensitivityMeanPred*100,
                sensitivityLPred=sensitivityLPred*100,
                sensitivityHPred=sensitivityHPred*100)

nPlots <- 2
nTests <- length(unique(validationDf$testName))
testLists <- split(c(1:nTests), cut(seq_along(c(1:nTests)), nPlots,
                                    labels = FALSE)) 
validationPlot <- list()
for (np in c(1:nPlots)) {
  testsPlot <- unique(validationDf$testName)[testLists[[np]]]
  validationPlot[[np]] <- dplyr::filter(validationDf, testName %in% testsPlot) %>%
    ggplot(aes(x=sensitivityMean, y=sensitivityMeanPred, color=citationID,
               shape=sampleType)) +
    geom_pointrange(aes(ymin=sensitivityLPred, ymax=sensitivityHPred),
                    position=position_jitter(width=0.15, height=0)) +
    geom_pointrange(aes(xmin=sensitivityL, xmax=sensitivityH),
                    position=position_jitter(width=0.15, height=0)) +
    geom_abline(slope=1, intercept=0) +
    facet_wrap(.~testName, ncol=3) +
    theme_bw() +
    theme(legend.position="top") +
    guides(color=FALSE, shape=guide_legend("Sample type")) +
    xlab("Observed sensitivity (%)") +
    ylab("Predicted sensitivity (%)")
ggsave(paste("../data/figures/sensitivity_validation", np, ".png", sep=""),
       validationPlot[[np]], units="cm", width=24, height=30)
}

validationTimePlot <- list()
for (np in c(1:nPlots)) {
  testsPlot <- unique(validationDf$testName)[testLists[[np]]]
  validationTimePlot[[np]] <- dplyr::filter(validationDf, testName %in% testsPlot) %>%
    ggplot(aes(x=testTime, y=sensitivityMean, shape=sampleType)) +
    geom_pointrange(aes(ymin=sensitivityL, ymax=sensitivityH),
                    color="blue", alpha=0.4) +
    geom_segment(aes(x=testTime, xend=testTime+0.2, y=sensitivityMean,
                     yend=sensitivityMeanPred), color="red") +
    facet_wrap(.~testName, ncol=3) +
    theme_bw() +
    theme(legend.position="top") +
    guides(color=FALSE, shape=guide_legend("Sample type")) +
    xlim(0, 15) +
    xlab("Diagnosis to test (months)") +
    ylab("Sensitivity (%)")
ggsave(paste("../data/figures/sensitivity_validation_time", np, ".png", sep=""),
       validationTimePlot[[np]], units="cm", width=24, height=30)
}

# test how many times the observed sensitivity falls within the
# uncertainty interval
sensSamples <- 500
studySamples <- 10
predL <- NULL
predH <- NULL
for (r in c(1:nrow(validationDf))) {
  sdSens <- with(validationDf, (sensitivityHPred[r]-sensitivityLPred[r])/3.92)
  nTests <- validationDf$nSamples[r]
  testSens <- NULL
  for (s in c(1:sensSamples)) {
    sampleSens <- rnorm(1, mean=validationDf$sensitivityMeanPred[r], sd=sdSens)
    sampleSens <- pmin(pmax(sampleSens, 0), 100)
    testSens <- c(rbinom(studySamples, size=validationDf$nSamples[r],
                         prob=sampleSens/100)/nTests, testSens)
  }
  predL[r] <- quantile(testSens, 0.025) * 100
  predH[r] <- quantile(testSens, 0.975) * 100
}

inCI <- (validationDf$sensitivityMean < predH) & (validationDf$sensitivityMean > predL)
mean(inCI)


#######################
#######################
## 5) Plot assay characteristics parameters
#######################
#######################

############
# Antigen effect
############
antigenSamplesDf <- read.csv("../data/analysis_results/7_assay_characteristics_antigenTarget.csv",
                             stringsAsFactors=FALSE)

############## Remove after fixing name in data
antigenPlot <- dplyr::filter(antigenSamplesDf, !is.na(parName)) %>%
  ggplot(data=., aes(x=.value, fill=parName)) +
  geom_density(alpha=0.3) +
#  facet_wrap(.~parName, ncol=1) +
  theme_bw() +
  theme(axis.title.y=element_blank()) +
  guides(fill=guide_legend("Parameter")) +
  xlab("Value")
ggsave("../data/figures/characteristicValue_Antigen.png", antigenPlot,
       units="cm", width=16, height=12)

# Get how many times a parameter is larger than the other, as statistical
# test
antigenWide <- dplyr::filter(antigenSamplesDf, .variable=="charSlope") %>%
  pivot_wider(., id_cols=.draw, names_from=parName, values_from=.value)
s_nMean <- mean(antigenWide$S - antigenWide$N)
s_nCI <- quantile(antigenWide$S - antigenWide$N, probs=c(0.025, 0.975))
rbdMean <- mean(antigenWide$S + antigenWide$RBD)
rbdCI <- quantile(antigenWide$S + antigenWide$RBD, probs=c(0.025, 0.975))
snMean <- mean(antigenWide$SN)
snCI <- quantile(antigenWide$SN, probs=c(0.025, 0.975))

############
# test effect
############
techniqueSamplesDf <- read.csv("../data/analysis_results/7_assay_characteristics_technique.csv",
                             stringsAsFactors=FALSE)

techniquePlot <- dplyr::filter(techniqueSamplesDf, !is.na(parName)) %>%
  ggplot(data=., aes(x=.value, fill=parName)) +
  geom_density(alpha=0.3) +
#  facet_wrap(.~parName, ncol=1) +
  theme_bw() +
  theme(axis.title.y=element_blank()) +
  guides(fill=guide_legend("Parameter")) +
  xlab("Value")
ggsave("../data/figures/characteristicValue_Technique.png", techniquePlot,
       units="cm", width=14, height=10)

techniqueWide <- dplyr::filter(techniqueSamplesDf, .variable=="charSlope") %>%
  pivot_wider(., id_cols=.draw, names_from=parName, values_from=.value)
#lfa_chem <- mean(techniqueWide$LFA > techniqueWide$Chemoluminiscense)
#lfa_elisa <- mean(techniqueWide$LFA > techniqueWide$ELISA)
#chem_elisa <- mean(techniqueWide$Chemoluminiscense > techniqueWide$ELISA)
lfa_restMean <- mean(techniqueWide$LFA - techniqueWide$Rest)
lfa_restCI <- quantile(techniqueWide$LFA-techniqueWide$Rest, probs=c(0.025, 0.975))

############
# antibody effect
############
antibodySamplesDf <- read.csv("../data/analysis_results/7_assay_characteristics_antibody.csv",
                             stringsAsFactors=FALSE)

antibodyPlot <- dplyr::filter(antibodySamplesDf, !is.na(parName)) %>%
  ggplot(data=., aes(x=.value, fill=parName)) +
  geom_density(alpha=0.3) +
#  facet_wrap(.~parName, ncol=1) +
  theme_bw() +
  theme(axis.title.y=element_blank()) +
  guides(fill=guide_legend("Parameter")) +
  xlab("Value")
ggsave("../data/figures/characteristicValue_Antibody.png", antibodyPlot,
       units="cm", width=14, height=10)

antibodyWide <- dplyr::filter(antibodySamplesDf, .variable=="charSlope") %>%
  pivot_wider(., id_cols=.draw, names_from=parName, values_from=.value)
#lfa_chem <- mean(antibodyWide$LFA > antibodyWide$Chemoluminiscense)
#lfa_elisa <- mean(antibodyWide$LFA > antibodyWide$ELISA)
#chem_elisa <- mean(antibodyWide$Chemoluminiscense > antibodyWide$ELISA)
igmMean <- mean(antibodyWide$IgM)
igmCI <- quantile(antibodyWide$IgM, probs=c(0.025, 0.975))
igaMean <- mean(antibodyWide$IgA)
igaCI <- quantile(antibodyWide$IgA, probs=c(0.025, 0.975))
totMean <- mean(antibodyWide$Total)
totCI <- quantile(antibodyWide$Total, probs=c(0.025, 0.975))

############
# antibody22 effect
############
antibody2SamplesDf <- read.csv("../data/analysis_results/7_assay_characteristics_antibody2.csv",
                             stringsAsFactors=FALSE)

antibody2Plot <- dplyr::filter(antibody2SamplesDf, !is.na(parName)) %>%
  ggplot(data=., aes(x=.value, fill=parName)) +
  geom_density(alpha=0.3) +
#  facet_wrap(.~parName, ncol=1) +
  theme_bw() +
  theme(axis.title.y=element_blank()) +
  guides(fill=guide_legend("Parameter")) +
  xlab("Value")

ggsave("../data/figures/characteristicValue_Antibody2.png", antibody2Plot,
       units="cm", width=14, height=10)

antibody2Wide <- dplyr::filter(antibody2SamplesDf, .variable=="charSlope") %>%
  pivot_wider(., id_cols=.draw, names_from=parName, values_from=.value)
#lfa_chem <- mean(antibody2Wide$LFA > antibody2Wide$Chemoluminiscense)
#lfa_elisa <- mean(antibody2Wide$LFA > antibody2Wide$ELISA)
#chem_elisa <- mean(antibody2Wide$Chemoluminiscense > antibody2Wide$ELISA)
otherMean <- mean(antibody2Wide$Others)
otherCI <- quantile(antibody2Wide$Others, probs=c(0.025, 0.975))



#######################
#######################
## 6) Average sensitivities by assay characteristics
#######################
#######################
assayChars <- read.csv("../data/raw_data/assay_characteristics.csv",
                       stringsAsFactors=FALSE, na.strings=c("", "--")) %>%
  dplyr::rename(., testName=test_name_long, antigen_target=antigen.target,
  assay_type=assay.type, isotype_target=isotype.target)

# First, add assay characteristics to the dataframe with the fitted data
antigenVec <- NULL
antibodyVec <- NULL
techniqueVec <- NULL
for (r in c(1:nrow(seroFitted))) {
  charRow <- assayChars$testName == seroFitted$testName[r]
  if (any(charRow )) {
    antigenVec <- c(antigenVec, assayChars$antigen_target[charRow])
    antibodyVec <- c(antibodyVec, assayChars$isotype_target[charRow])
    techniqueVec <- c(techniqueVec, assayChars$assay_type[charRow])
  } else {
    antigenVec <- c(antigenVec, NA)
    antibodyVec <- c(antibodyVec, NA)
    techniqueVec <- c(techniqueVec, NA)
  }
}

seroFitted$antigenTarget <- antigenVec
seroFitted$antibodyTarget <- antibodyVec
seroFitted$technique <- techniqueVec

# Add binary columns indicating presence or absence of characteristics
seroFitted$spikeAntigen <- (str_detect(seroFitted$antigenTarget, "S") |
  str_detect(seroFitted$antigenTarget, "RBD"))
seroFitted$rbdAntigen <- str_detect(seroFitted$antigenTarget, "RBD")
seroFitted$nuceloAntigen <- str_detect(seroFitted$antigenTarget, "N")
seroFitted$LFA <- str_detect(seroFitted$technique, "LFIA")
seroFitted$ELISA <- str_detect(seroFitted$technique, "ELISA")
seroFitted$chemlum <- (str_detect(seroFitted$technique, "CMIA") |
                       str_detect(seroFitted$technique, "CLIA"))
seroFitted$IgM <- str_detect(seroFitted$antibodyTarget, "IgM")
seroFitted$IgA <- str_detect(seroFitted$antibodyTarget, "IgA")

###########
# Antigen target sensitivity
###########

# match labels of fitted data to labels of model levels
# The new labels indicate the name of the parameters that
# are "activated" for a given assay, separated by "-".
# It is important that these labels match the labels
# generated by "charCombs", in order for the code to work
spike_RBD <- with(seroFitted, !is.na(antigenTarget) &
                  spikeAntigen & rbdAntigen & !(nuceloAntigen))
seroFitted$antigenTarget[spike_RBD] <- "S-RBD"
spike_N <- with(seroFitted, !is.na(antigenTarget) &
                  spikeAntigen & !(rbdAntigen) & nuceloAntigen)
seroFitted$antigenTarget[spike_N] <- "S-N-SN"
all_antigen <- with(seroFitted, !is.na(antigenTarget) &
                  spikeAntigen & rbdAntigen & nuceloAntigen)
seroFitted$antigenTarget[all_antigen] <- "S-N-SN-RBD"

antigenRef <- dplyr::select(seroFitted, testName, antigenTarget) %>%
  unique(.)

######
# Mean sensitivity across time for each param combination
######
timeVec <- seq(0.5, 12, 0.2)
charCombs <- list("S", "N", c("S", "N", "SN"), c("S", "RBD"),
                  c("S", "N", "SN", "RBD"))
antigenSens <- mean_param_seroreversion(antigenSamplesDf, timeVec,
                            slopeName="charSlope",
                            interceptName="intercept",
                            charCombs=charCombs,
                            timeNormalization=4)

# Get sensitivity across time for each assay fitted
assayRefDf <- unique(select(seroFitted, testName, antigenTarget)) %>%
  dplyr::rename(., chars=antigenTarget)
assayAntigenSerorev <- assay_char_sensitivity_posterior(antigenSamplesDf,
                                 timeVec=timeVec,
                                 charCombs=charCombs,
                                 assayRefDf=assayRefDf,
                                 timeNormalization=4)

regLineSize <- 2

seroFitted$chars <- seroFitted$antigenTarget
antigenSensPlot <- antigenSens %>%
  ggplot(., aes(x=time, y=sensitivity*100, color=chars, fill=chars)) +
  geom_point(data=dplyr::filter(seroFitted, !is.na(antigenTarget)),
             aes(x=testTime, y=sensitivityMean, size=sqrt(nSamples)/10),
             color="black", alpha=0.3) +
  geom_line(data=assayAntigenSerorev, aes(x=time, y=sensitivity*100, group=testName),
            color="black", alpha=0.3) +
  geom_line(size=regLineSize) +
  geom_ribbon(aes(ymin=sensitivityL*100, ymax=sensitivityH*100),
              alpha=ribbonAlpha) +
  facet_wrap(.~chars, ncol=2, labeller=label_wrap_gen(width=35)) +
  theme_bw() +
  theme(legend.position="top") +
  guides(color=FALSE, size=FALSE, fill=FALSE) +
  xlim(0, 15) +
  xlab("Diagnosis to test (months)") +
  ylab("Sensitivity (%)")

ggsave("../data/figures/antigen_sensitivity_time.png", antigenSensPlot,
       units="cm", width=20, height=25)

###########
# Test technique
###########
timeVec <- seq(0.5, 14, 0.2)
charCombs <- list("LFA", "Rest")
techniqueSens <- mean_param_seroreversion(techniqueSamplesDf, timeVec,
                            slopeName="charSlope",
                            interceptName="intercept",
                            charCombs=charCombs,
                            timeNormalization=4)

# Get sensitivity across time for each assay fitted
seroFitted$technique[which(seroFitted$LFA)] <- "LFA"
seroFitted$technique[which(!seroFitted$LFA)] <- "Rest"

assayRefDf <- unique(select(seroFitted, testName, technique)) %>%
  dplyr::rename(., chars=technique)
assayTechniqueSerorev <- assay_char_sensitivity_posterior(techniqueSamplesDf,
                                 timeVec=timeVec,
                                 charCombs=charCombs,
                                 assayRefDf=assayRefDf,
                                 timeNormalization=4)

seroFitted$chars <- seroFitted$technique
techniqueSensPlot <- techniqueSens %>%
  ggplot(., aes(x=time, y=sensitivity*100, color=chars, fill=chars)) +
  geom_point(data=dplyr::filter(seroFitted, !is.na(chars)),
             aes(x=testTime, y=sensitivityMean, size=sqrt(nSamples)/10),
             color="black", alpha=0.3) +
  geom_line(data=assayTechniqueSerorev, aes(x=time, y=sensitivity*100, group=testName),
            color="black", alpha=0.3) +
  geom_line(size=regLineSize) +
  geom_ribbon(aes(ymin=sensitivityL*100, ymax=sensitivityH*100),
              alpha=ribbonAlpha) +
  facet_wrap(.~chars, ncol=2, labeller=label_wrap_gen(width=35)) +
  theme_bw() +
  theme(legend.position="top") +
  guides(color=FALSE, size=FALSE, fill=FALSE) +
  xlim(0, 15) +
  xlab("Diagnosis to test (months)") +
  ylab("Sensitivity (%)")

ggsave("../data/figures/technique_sensitivity_time.png", techniqueSensPlot,
       units="cm", width=20, height=14)

###########
# Antibody
###########
timeVec <- seq(0.5, 14, 0.2)
charCombs <- list("IgG", c("IgG", "IgM"), c("IgG", "IgA"),
                  c("IgG", "IgM", "IgA", "Total"))
antibodySens <- mean_param_seroreversion(antibodySamplesDf, timeVec,
                            slopeName="charSlope",
                            interceptName="intercept",
                            charCombs=charCombs,
                            timeNormalization=4)

# Get sensitivity across time for each assay fitted
seroFitted$antibodyTarget <- stringr::str_replace(seroFitted$antibodyTarget,
                                                  " etc.", "") %>%
  stringr::str_replace(., "IgG IgM IgA", "IgG-IgM-IgA-Total") %>%
  stringr::str_replace(., " ", "-")

assayRefDf <- unique(select(seroFitted, testName, antibodyTarget)) %>%
  dplyr::rename(., chars=antibodyTarget)
assayAntibodySerorev <- assay_char_sensitivity_posterior(antibodySamplesDf,
                                 timeVec=timeVec,
                                 charCombs=charCombs,
                                 assayRefDf=assayRefDf,
                                 timeNormalization=4)

seroFitted$chars <- seroFitted$antibodyTarget
antibodySensPlot <- antibodySens %>%
  ggplot(., aes(x=time, y=sensitivity*100, color=chars, fill=chars)) +
  geom_point(data=dplyr::filter(seroFitted, !is.na(chars)),
             aes(x=testTime, y=sensitivityMean, size=sqrt(nSamples)/10),
             color="black", alpha=0.3) +
  geom_line(data=assayAntibodySerorev, aes(x=time, y=sensitivity*100, group=testName),
            color="black", alpha=0.3) +
  geom_line(size=regLineSize) +
  geom_ribbon(aes(ymin=sensitivityL*100, ymax=sensitivityH*100),
              alpha=ribbonAlpha) +
  facet_wrap(.~chars, ncol=2, labeller=label_wrap_gen(width=35)) +
  theme_bw() +
  theme(legend.position="top") +
  guides(color=FALSE, size=FALSE, fill=FALSE) +
  xlim(0, 15) +
  xlab("Diagnosis to test (months)") +
  ylab("Sensitivity (%)")

ggsave("../data/figures/antibody_sensitivity_time.png", antibodySensPlot,
       units="cm", width=20, height=12)




#######################
#######################
## 7) Full model
#######################
#######################

validationDf <- read.csv("../data/analysis_results/7_assay_characteristics_fullModel_validation.csv",
                         stringsAsFactors=FALSE)

good <- filter(validationDf, (nSeropositives>=predictedPosL) &
               (nSeropositives<=predictedPosH))
bad <- filter(validationDf, !((nSeropositives>=predictedPosL) &
               (nSeropositives<=predictedPosH)))




