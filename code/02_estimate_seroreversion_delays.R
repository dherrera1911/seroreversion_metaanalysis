###############################
# 
# This script estimates the time delays between
# diagnosis and serological testing, using the time series of
# case reports collected in script 01.
# 
# The script generates as output the file
# "../data/processed_data/PCR_to_serotest_estimated_times.csv", which
# contains the reported data together with the times of
# diagnosis-to-testing estimated here.
#
# Script authored by Daniel Herrera-Esposito.
# For questions, contact me at dherrera1911[at]gmail.com
# 
# Final version revised 10/03/2023
# 
###############################

library(dplyr)
library(tidyr)
library(lubridate)

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

# make names of serotest table compatible with case dynamics table
newLocs <- c("Ischgl", "Spain", "Salt Lake",
             "Ethiopia", "India", "Kashmir", "Kashmir", "India",
             "Ireland", "Lithuania", "Tirschenreuth",
             "Arkansas", "Arkansas", "Arkansas",
             "United States", "Houston", "Canada", "Madryn",
             "Puducherry", "Puducherry", "Puducherry", "Indonesia",
             "Hyderabad", "Japan", "Japan",
             "Delhi", "Delhi", "Ischgl",
             "France", "France regions", "Kupferzell", "Feilnbach",
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
             "11 Cantons", "Portugal", "India", "India",
             "Peru", "Lebanon", 
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

write.csv(seroUnknown, "../data/processed_data/PCR_to_serotest_estimated_times.csv",
          row.names=FALSE)

