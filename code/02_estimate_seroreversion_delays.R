library(dplyr)
library(tidyr)
library(lubridate)
library(rstan)
library(bayesplot)
library(tidybayes)

source("./functions_auxiliary.R")

############
# Data of death and case dynamics
############
dynamics <- read.csv("../data/raw_data/death_dynamics.csv", stringsAsFactors=FALSE) %>%
  dplyr::mutate(., date=date(date)) %>%
  as_tibble(.)

############
# Data of PCR test, serological retest, unknown time
############
seroUnknown <- read.csv("../data/raw_data/PCR_to_serotest_unknown_times.csv",
                     stringsAsFactors=FALSE) %>%
  dplyr::mutate(., test.name=str_replace(test.name, " $", ""),
                number.of.seropositives.among.prior.positives=
                  as.integer(number.of.seropositives.among.prior.positives),
                number.of.prior.positives=as.integer(number.of.prior.positives)) %>%
  as_tibble(.)
# Change variable names to be more friendly
newVarNames <- c("phase_id", "region", "country", "location", "sampleType",
                 "startDate", "endDate", "midpointDate", 
                 "testName", "citationID", "nSeropositives",
                 "nSamples", "sensitivityMean", "sensitivityL",
                 "sensitivityH", "notes")
names(seroUnknown) <- newVarNames
# Make dates dates
seroUnknown$midpointDate <- lubridate::mdy(seroUnknown$midpointDate)

# make names of serotest table compatible with dynamics table
newLocs <- c("Ischgl", "Spain", "Salt Lake",
             "Ethiopia", "Tamil Nadu", "India", "Tamil Nadu", "Kashmir",
             "Kashmir", "India", "Ireland", "Lithuania", "Tirschenreuth",
             "Arkansas", "Arkansas", "Arkansas",
             "United States", "Houston", "Canada", "Madryn",
             "Puducherry", "Puducherry", "Puducherry", "Indonesia",
             "Hyderabad", "Japan", "Japan",
             "Delhi", "Delhi", "Ischgl",
             "France", "France", "Kupferzell", "Feilnbach",
             "Pune", "Reutlingen",
             "Freiburg", "Straubing", "Aachen",
             "Reutlingen", "Osnabruck", "Mitte", "Freiburg", "Magdeburg",
             "Aachen", "Osnabruck", "Magdeburg", "Chemnitz",
             "Vorpommern", "England",
             "England", "England", "England", "Cantabria", "Gauteng",
             "Wuhan", "Wuhan", "Wuhan", "Tunisia", "Delhi",
             "Andorra", "Chile", "Los Angeles", "Holyoke",
             "Hillsborough", "Tyumen", "Khabarovsk",
             "Saint Petersburg", "Leningrad", "Sverdlovsk",
             "Tatarstan", "Moscow", "Chelyabinsk", "Irkutsk",
             "Saratov", "Kaliningrad", "Murmansk", "Krasnoyarsk",
             "Novosibirsk", "Stavropol",
             "Spain", "Spain",  "Spain",  "Spain",
             "Portugal", "India", "India", "Lebanon", 
             "Connecticut", "Delhi", "Delhi")

new <- c("Arkansas", "Japan", "Feilnbach")

seroUnknown$newLoc <- newLocs
seroUnknown$pcrReferenceDate <- seroUnknown$midpointDate

# Notes:
# 1) For Ischgl, approximate median infection date from paper = 2020-03-16
# 2) For Pollán, Spain, add results of <= 14 days
# 3) Added symptoms date to Castiglioni D'Adda. Not sure if to include,
#   possible bias due to confirming testing of rapid test positive.

# Remove last 14 days from Chennai, since they exclude those in their sample
# both studies?
chennaiInd1 <- which(seroUnknown$location == "Tamil Nadu: Chennai")[1]
seroUnknown$pcrReferenceDate[chennaiInd1] <- seroUnknown$pcrReferenceDate[chennaiInd1] - 14
# Check india source
indiaInd1 <- which(seroUnknown$country == "India" & (seroUnknown$location %in%
                    c("National", "National (unvaccinated population)")))[1]
seroUnknown$pcrReferenceDate[indiaInd1] <- seroUnknown$pcrReferenceDate[indiaInd1] - 14
# put a midpoint date to Hillsborough
hillStr <- "Hillsborough County, unvaccinated population"
hillsboroughInd <- which(seroUnknown$location == hillStr)
seroUnknown$pcrReferenceDate[hillsboroughInd] <- date("2021-01-01")
seroUnknown$midpointDate[hillsboroughInd] <- date("2021-01-01")
# put a midpoint date to Tunisia
tunisiaInd <- which(seroUnknown$country == "Tunisia")
seroUnknown$pcrReferenceDate[tunisiaInd] <- date("2021-04-01")
seroUnknown$midpointDate[tunisiaInd] <- date("2021-04-01")

# Estimate the delays for each row
estimatedTime <- NULL
for (nr in c(1:nrow(seroUnknown))) {
  seroLoc <- seroUnknown$newLoc[nr]
  locDynamics <- dplyr::filter(dynamics, Location==seroLoc)
  locDynamics <- dplyr::filter(locDynamics, date<=seroUnknown$pcrReferenceDate[nr])
  maxCases <- locDynamics$cases[nrow(locDynamics)]
  medianInd <- which(locDynamics$cases>=(maxCases/2))[1]
  medianDay <- locDynamics$date[medianInd]
  estimatedTime[nr] <- (lubridate::interval(medianDay, seroUnknown$midpointDate[nr])/
    months(1))
}

if (any(estimatedTime==0)) {
  estimatedTime[estimatedTime==0] <- 1
}

# Change date for Castiglioni D'Adda
#castInd <- which(seroUnknown$newLoc=="Adda")
#meanSymptomDate <- date("2020-03-01")
#estimatedTime[castInd] <- (lubridate::interval(meanSymptomDate,
#                                               seroUnknown$midpointDate[castInd])/
#    months(1))
# Change date for Odisha Citation 93
odiInd <- which(seroUnknown$newLoc=="Odisha" & seroUnknown$phase_id==2)
estimatedTime[odiInd] <- 1
# Change date for Ischgl citation 149. Paper reports more geographically
# precise dates
ischglInd <- which(seroUnknown$newLoc=="Ischgl")
estimatedTime[ischglInd] <- lubridate::interval(date("2020-03-15"),
                                                  seroUnknown$midpointDate[ischglInd])/months(1)
 
seroUnknown$testTime <- pmax(round(estimatedTime), 1)
seroUnknown <- dplyr::select(seroUnknown, -newLoc, -pcrReferenceDate)

write.csv(seroUnknown, "../data/analysis_results/PCR_to_serotest_estimated_times.csv",
          row.names=FALSE)



############
# Data of serological test and serological retest
############
sero2sero <- read.csv("../data/raw_data/serotest_to_serotest.csv",
                     stringsAsFactors=FALSE) %>%
  dplyr::mutate(., test.name.for.later.sampling=
                str_replace(test.name.for.later.sampling, " $", ""),
                number.of.later.seropositives.among.initial.seropositives=
                  as.integer(number.of.later.seropositives.among.initial.seropositives),
                number.of.initial.seropositives=as.integer(number.of.initial.seropositives)) %>%
  as_tibble(.)

# Change variable names to be more friendly
newVarNames <- c("phase_id", "country", "location", "sampleType",
                 "midpointInitial", "midpointFollowup", "testTime",
                 "testName", "citationID", "nSeropositives",
                 "nSamples", "sensitivityMean", "sensitivityL",
                 "sensitivityH", "includedTable", "notes")

names(seroUnknown) <- newVarNames

# Remove mixed assays
seroUnknown <- dplyr::filter(seroUnknown, !stringr::str_detect(testName, "OR"))

# Give confidence intervals to datapoints lacking them (only N samples and N positive)
seroUnknown$sensitivityL <- as.numeric(seroUnknown$sensitivityL)
seroUnknown$sensitivityH <- as.numeric(seroUnknown$sensitivityH)
for (r in c(1:nrow(seroUnknown))) {
  if (is.na(seroUnknown[[r,"sensitivityL"]])) {
    seroprevConfint <- binomial_confint(seroUnknown[[r,"nSamples"]],
                                        seroUnknown[[r, "nSeropositives"]])
    seroUnknown[[r,"sensitivityL"]] <- signif(seroprevConfint$lower*100, digits=4)
    seroUnknown[[r,"sensitivityH"]] <- signif(seroprevConfint$upper*100, digits=4)
    seroUnknown[[r,"sensitivityMean"]] <- signif(seroUnknown[[r,"nSeropositives"]]/
                                              seroUnknown[[r,"nSamples"]]*100, digits=4)
  }
}


# Relabel test times that give intervals into single times
seroUnknown$testTime[seroUnknown$testTime=="1 - 2"] <- "1.5"
seroUnknown$testTime[seroUnknown$testTime=="3 - 4"] <- "3.5"
seroUnknown$testTime[seroUnknown$testTime=="≥5"] <- "6"
# Convert test times to integer
seroUnknown$testTime <- as.numeric(seroUnknown$testTime)




