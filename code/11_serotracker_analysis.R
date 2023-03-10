##################################
# 
# This script analyzes how many of the data points in
# the SeroTracker dataset are at high risk of seroreversion
# bias, as indicated by the assay used.
# 
# We note that we do not account for other factors, such as whether
# the authors in the study actually correct for seroreversion
# in some way, or the delay between epidemic wave and the
# serological survey.
# 
# At the end of the script, we obtain the values that are reported in
# Table 2 in the associated article:
# https://www.medrxiv.org/content/10.1101/2022.09.08.22279731v3
# now in press at Eurosurveillance
# 
# Script authored by Daniel Herrera-Esposito.
# For questions, contact me at dherrera1911[at]gmail.com
# 
# Final version revised 10/03/2023
# 
##################################

library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)

seroTrack <- read.csv("../data/serotracker_dataset.csv",
                       stringsAsFactors=FALSE) %>%
  tidyr::as_tibble(.) %>%
  dplyr::filter(., is_unity_aligned=="Unity-Aligned")

seroTrack$sampling_end_date <- lubridate::ymd(seroTrack$sampling_end_date)
seroTrack$sampling_start_date <- lubridate::ymd(seroTrack$sampling_start_date)

agTargets <- seroTrack$antibody_target
agIsotype <- seroTrack$isotypes
seroTrack$sampling_mid_date <- lubridate::today()
for (r in c(1:nrow(seroTrack))){
  agTargets[r] <- str_replace_all(agTargets[r], pattern="'", replacement="")
  agTargets[r] <- str_replace_all(agTargets[r], pattern="\\[", replacement="")
  agTargets[r] <- str_replace_all(agTargets[r], pattern="\\]", replacement="")
  agTargets[r] <- str_replace_all(agTargets[r], pattern=",", replacement="")
  agTargets[r] <- str_replace_all(agTargets[r], pattern="  ", replacement="_")
  agTargets[r] <- str_replace_all(agTargets[r], pattern=" ", replacement="")
  agTargets[r] <- str_replace_all(agTargets[r], pattern="__", replacement="_")
  agIsotype[r] <- str_replace_all(agIsotype[r], pattern="'", replacement="")
  agIsotype[r] <- str_replace_all(agIsotype[r], pattern="\\[", replacement="")
  agIsotype[r] <- str_replace_all(agIsotype[r], pattern="\\]", replacement="")
  agIsotype[r] <- str_replace_all(agIsotype[r], pattern=",", replacement="")
  agIsotype[r] <- str_replace_all(agIsotype[r], pattern="  ", replacement="_")
  agIsotype[r] <- str_replace_all(agIsotype[r], pattern=" ", replacement="")
  agIsotype[r] <- str_replace_all(agIsotype[r], pattern="__", replacement="_")
  seroTrack$sampling_mid_date[r] <- with(seroTrack,
                                         mean(c(sampling_end_date[r], sampling_start_date[r])))
}

seroTrack$antibody_target <- agTargets
seroTrack$isotypes <- agIsotype

seroTrack <- dplyr::filter(seroTrack, sampling_mid_date < "2022-05-01")

seroTrack <- dplyr::select(seroTrack, antibody_target, sampling_mid_date,
  test_name, test_manufacturer, source_name, test_type) %>%
  unique()

tests <- unique(seroTrack$test_name)
testType <- unique(seroTrack$test_type)
testTarget <- seroTrack$antibody_target

### Indicate if assay targets only nucleocapsid
seroTrack$onlyN <- grepl("^Nucleocapsid.*\\(N-protein)$", seroTrack$antibody_target)
### Indicate if assay targeting only N is a known sandwich assay
# Get the names of assays targeting only nucleocapsid
assayNamesN <- sort(unique(seroTrack$test_name[seroTrack$onlyN]))
assayManN <- sort(unique(seroTrack$test_manufacturer[seroTrack$onlyN]))
directAssayNames_N <- c("*Elecsys*", "*Neutralization*", "*Platelia*")
isNDirect <- FALSE
for (dn in directAssayNames_N){
  isNDirect <- (grepl(dn, seroTrack$test_name) | isNDirect) & seroTrack$onlyN
}
seroTrack$directN <- isNDirect

seroTrack$highRisk <- with(seroTrack, (onlyN & !directN) | (test_type=="LFIA"))

# Get an estimate of all studies in meta-analysis with high-risk of
# seroreversion
# Remove repeated data points from same study
unityStudies <- dplyr::select(seroTrack, source_name, highRisk) %>%
  unique(.)

# Get an estimate of data points in 1st half 2020 with high seroreversion potential
sero2020_1 <- dplyr::filter(seroTrack, sampling_mid_date<"2020-07-01") %>%
  dplyr::select(., source_name, highRisk)
mean(sero2020_1$highRisk)
sum(sero2020_1$highRisk)
length(sero2020_1$highRisk)

# Get an estimate of data points in 2nd half 2020 with high seroreversion potential
sero2020_2 <- dplyr::filter(seroTrack, sampling_mid_date>="2020-07-01" &
                            sampling_mid_date<"2021-01-01")
  dplyr::select(., source_name, highRisk)
mean(sero2020_2$highRisk)
sum(sero2020_2$highRisk)
length(sero2020_2$highRisk)

# Get an estimate of data points in 1st half 2021 with high seroreversion potential
sero2021_1 <- dplyr::filter(seroTrack, sampling_mid_date >= "2021-01-01" &
                            sampling_mid_date<"2021-07-01") %>%
  dplyr::select(., source_name, highRisk)
mean(sero2021_1$highRisk)
sum(sero2021_1$highRisk)
length(sero2021_1$highRisk)

# Get an estimate of data points in 2nd half 2021 with high seroreversion potential
sero2021_2 <- dplyr::filter(seroTrack, sampling_mid_date>="2021-07-01" &
                            sampling_mid_date<"2022-01-01")
  dplyr::select(., source_name, highRisk)
mean(sero2021_2$highRisk)
sum(sero2021_2$highRisk)
length(sero2021_2$highRisk)

