library(dplyr)
library(lubridate)
library(stringr)
library(ggplot2)
library(ggpubr)
library(matrixStats)

# return the 95% CI of a binomial coefficient
binomial_confint <- function(countTotal, occurrences, input="count"){
  lower <- NULL
  upper <- NULL
  if (input=="count") {
    countOccurrences <- occurrences
  } else {
    countOccurrences <- round(occurrences * countTotal)
  }
  for (i in c(1:length(countTotal))) {
    confint <- binom.test(countOccurrences[i], countTotal[i])$conf.int
    lower[i] <- confint[1]
    upper[i] <- confint[2]
  }
  confintList <- list(lower=lower, upper=upper)
}

# function used for plotting log axis
scaleFun <- function(x) sprintf("%1g", x)

# Return the parameters of a beta distribution, as
# estimated from moment matching its mean and 95% CI bounds
fit_beta_ci <- function(meanEstimate, lower, upper) {
  var <- ((upper-lower)/4)^2
  commonFactor <- meanEstimate*(1-meanEstimate)/var - 1
  shape1 <- meanEstimate*commonFactor
  shape2 <- (1-meanEstimate)*commonFactor
  fittingDf <- data.frame(meanEstimate=meanEstimate, lower=lower, upper=upper,
                    shape1=NA, # stores estimates of shape1
                    shape2=NA, # same for shape 2
                    priorMean=NA, # the mean of the fitted prior
                    priorLB=NA,  # the 0.025 quantile of the fitted prior
                    priorUB=NA) # the 0.975 quantile of the fitted prior
  fittingDf$shape1 <- shape1
  fittingDf$shape2 <- shape2
  fittingDf$priorMean <- shape1/(shape1+shape2)
  for (i in c(1:nrow(fittingDf))) {
    shape1i <- fittingDf$shape1[i]
    shape2i <- fittingDf$shape2[i]
    fittingDf$priorLB[i] <- qbeta(p=.025, shape1=shape1i, shape2=shape2i)
    fittingDf$priorUB[i] <- qbeta(p=.975, shape1=shape1i, shape2=shape2i)
  }
  return(fittingDf)
}


# Make the batches for grouped-in-time cross-validation of the
# seroreversion fit.
#
# Each batch contains, for nAssayPerGroup assays, either
# all datapoints with times shorter than the average time
# for that test, or datapoints with times larger. That is,
# either the first half or the second half of the points for a test
# are left out in each batch. This way, this crossvalidation
# reflects extrapolation across time and assays
#
# Also, the function tries to balance the batches by pairing
# assays with high number of datapoints to those with low number
# of datapoints.
grouped_val_batches <- function(seroFitted, testsPerGroup) {
  # Get number of data points in each test, and sort to
  # match together high N tests with low N tests for CV
  testSamples <- group_by(seroFitted, testName) %>%
    summarize(., multiTime=length(unique(testTime)) > 1,
              nTimes=length(unique(testTime)),
              meanTime=mean(testTime))
  testSamples <- arrange(testSamples, nTimes)
  halfRows <- ceiling(nrow(testSamples)/2)
  testSamples[c(1:halfRows),] <- testSamples[c(halfRows:1),]
  # Make the CV group indexing
  nGroups <- ceiling(nrow(testSamples)/testsPerGroup)
  testGroup <- rep(c(1:nGroups), testsPerGroup)[1:nrow(testSamples)]
  validationGroup <- rep(NA, nrow(seroFitted))
  for (ng in c(1:nGroups)) {
    # tests belonging to this group
    tempTests <- testSamples$testName[testGroup==ng]
    for (tt in tempTests) {
      testInds <- which(seroFitted$testName == tt)
      meanTime <- testSamples$meanTime[testSamples$testName==tt]
      # Group data points by first or second half
      if (testSamples$multiTime[testSamples$testName==tt]) {
        timeGroup <- as.integer(seroFitted$testTime[testInds] >= meanTime) + 1
      } else {
        timeGroup <- 1
      }
      # Add group index (ng-1)*2 to the timeGroup
      validationGroup[testInds] <- (ng-1)*2 + timeGroup
    }
  }
  return(validationGroup)
}


# Make matrix with indicator variable for assay characteristics, to use for
# model fit
make_characteristics_matrix <- function(seroFitted, charsName) {
  nChars <- length(charsName)
  assayVec <- as.factor(seroFitted$testName)
  charLong <- NULL
  for (cn in charsName) {
    charLong <- c(charLong, as.integer(seroFitted[[cn]]))
  }
  charLong <- c(charLong, as.integer(assayVec))
  charMatrix <- matrix(charLong, ncol=nChars+1)
  charMatrix <- unique(charMatrix)
  charMatrix <- charMatrix[order(charMatrix[,nChars+1]),]
  charMatrix <- charMatrix[,c(1:nChars)]
  return(charMatrix)
}



