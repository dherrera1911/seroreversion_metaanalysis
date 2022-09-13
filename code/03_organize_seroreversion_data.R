library(dplyr)
library(tidyr)
library(lubridate)
source("./functions_auxiliary.R")

############
# 1) Load data with known diagnosis to serosurvey time
############
seroKnown <- read.csv("../data/raw_data/PCR_to_serotest_known_times.csv",
                     stringsAsFactors=FALSE) %>%
  dplyr::mutate(., test.name=str_replace(test.name, " $", ""),
                number.of.seropositives.among.prior.positives=
                  as.integer(number.of.seropositives.among.prior.positives),
                number.of.prior.positives=
                  as.integer(number.of.prior.positives),
                timeKnown=TRUE) %>%
  as_tibble(.)

# Change variable names to be more friendly
newVarNames <- c("phase_id", "region", "country", "location",
                 "sampleType", "startDate",
                 "endDate", "midpointDate", "testTime", "testName",
                 "citationID", "nSeropositives", "nSamples",
                 "sensitivityMean", "sensitivityL", "sensitivityH",
                 "includedTable", "notes", "timeKnown")
names(seroKnown) <- newVarNames

# Relabel test times that give intervals into single times
seroKnown$testTime[seroKnown$testTime=="0.5 - 2"] <- "1"
seroKnown$testTime[seroKnown$testTime=="1 - 2"] <- "1.5"
seroKnown$testTime[seroKnown$testTime=="1 - 3"] <- "2"
seroKnown$testTime[seroKnown$testTime=="2 - 3"] <- "2.5"
seroKnown$testTime[seroKnown$testTime=="2 - 4"] <- "3"
seroKnown$testTime[seroKnown$testTime=="3 - 4"] <- "3.5"
seroKnown$testTime[seroKnown$testTime=="4 - 5"] <- "4.5"
seroKnown$testTime[seroKnown$testTime=="4 - 6"] <- "5"
seroKnown$testTime[seroKnown$testTime=="5 - 7"] <- "6"
seroKnown$testTime[seroKnown$testTime==">6"] <- "7"
# Convert test times to integer
seroKnown$testTime <- as.numeric(seroKnown$testTime)

############
# 2) Data with estimated diagnosis to serosurvey time
############
seroEstimated <- read.csv("../data/processed_data/PCR_to_serotest_estimated_times.csv",
                          stringsAsFactors=FALSE) %>%
  dplyr::mutate(., timeKnown=FALSE)

seroEstimated$includedTable <- NA
seroEstimated <- dplyr::filter(seroEstimated, !(is.na(nSamples) & is.na(sensitivityL)))

#############
# 3) Put together datasets with known & estimated times
#############

seroAll <- rbind(seroKnown, seroEstimated)
# Select which dataset to fit
seroFitted <- seroAll
# Remove mixed assays
removedWords <- c("OR", "re-test", "\\+", "AND", "after")
for (w in c(1:length(removedWords))) {
  seroFitted <- dplyr::filter(seroFitted,
                              !stringr::str_detect(testName, removedWords[w]))
}

#############
# 4) Estimate some missing data
#############

# Give confidence intervals to datapoints lacking them (only N samples and N positive)
for (r in c(1:nrow(seroFitted))) {
  if (is.na(seroFitted[[r,"sensitivityL"]])) {
    seroprevConfint <- binomial_confint(seroFitted[[r,"nSamples"]],
                                        seroFitted[[r, "nSeropositives"]])
    seroFitted[[r,"sensitivityL"]] <- signif(seroprevConfint$lower*100, digits=4)
    seroFitted[[r,"sensitivityH"]] <- signif(seroprevConfint$upper*100, digits=4)
    seroFitted[[r,"sensitivityMean"]] <- signif(seroFitted[[r,"nSeropositives"]]/
                                              seroFitted[[r,"nSamples"]]*100, digits=4)
  }
}

# Fit beta distribution to studies without raw data, to estimate raw data
nonRaw <- is.na(seroFitted$nSamples)
fittedBetas <- fit_beta_ci(meanEstimate=seroFitted$sensitivityMean[nonRaw]/100,
                            lower=seroFitted$sensitivityL[nonRaw]/100,
                            upper=seroFitted$sensitivityH[nonRaw]/100)
totalCount <- with(fittedBetas, shape1+shape2)
meanSensitivity <- with(fittedBetas, shape1/(shape1+shape2)*100)
confintSensitivity <- binomial_confint(round(totalCount), round(fittedBetas$shape1))
# Add estimated raw data to the dataframe
seroFitted$nSamples[nonRaw] <- with(fittedBetas, round(shape1+shape2))
seroFitted$nSeropositives[nonRaw] <- round(fittedBetas$shape1)


############
# 5) Add assay characteristics to the dataframe
############
assayChars <- read.csv("../data/raw_data/assay_characteristics.csv",
                       stringsAsFactors=FALSE, na.strings=c("", "--")) %>%
  dplyr::rename(., testName=test_name_long, antigen_target=antigen.target,
  isotype_target=isotype.target, assay_type=assay.type, 
  test_design=test.design) %>%
  dplyr::filter(., testName %in% seroFitted$testName)

# add assay characteristics column to the data on sensitivity
antigenVec <- NULL
antibodyVec <- NULL
techniqueVec <- NULL
designVec <- NULL
for (r in c(1:nrow(seroFitted))) {
  charRow <- assayChars$testName == seroFitted$testName[r]
  if (any(charRow )) {
    antigenVec <- c(antigenVec, assayChars$antigen_target[charRow])
    antibodyVec <- c(antibodyVec, assayChars$isotype_target[charRow])
    techniqueVec <- c(techniqueVec, assayChars$assay_type[charRow])
    designVec <- c(designVec, assayChars$test_design[charRow])
  } else {
    antigenVec <- c(antigenVec, NA)
    antibodyVec <- c(antibodyVec, NA)
    techniqueVec <- c(techniqueVec, NA)
    designVec <- c(designVec, NA)
  }
}

seroFitted$antigenTarget <- antigenVec
seroFitted$antibodyTarget <- antibodyVec
seroFitted$technique <- techniqueVec
seroFitted$design <- designVec

# Filter tests that don't belong
seroFitted <- dplyr::filter(seroFitted, !stringr::str_detect(testName, "Obolensk"))

# save data to fit, for easier access
write.csv(seroFitted, "../data/processed_data/PCR_to_serotest_all.csv",
          row.names=FALSE)


