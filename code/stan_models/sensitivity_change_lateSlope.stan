data {
  int<lower=0> N;                           // number of observations
  int<lower=0> K;                           // number of assays
  int<lower=0> M;                           // number of studies:assays
  int<lower=0> C;                           // number of characteristics
  int<lower=1, upper=K> assay[N];           // assay index vector
  int<lower=1, upper=M> study[N];           // study:assay index vector
  matrix[K,C] characteristics;              // assay characteristics that we analyze
  vector[N] timeVecEarly;                   // predictor
  vector[N] timeVecLate;                    // predictor
  int<lower=0> nPositive[N];                // Number of positive samples
  int<lower=0> nTested[N];                  // Number of samples tested
}

parameters {
  // 
  vector[C] charSlope;                // slope effets of assay characteristics
  real meanSlopeLate;
  real intercept;                     // mean intercept
  real<lower=0> slopeSigmaEarly;           // sd of the slope
  real<lower=0> slopeSigmaLate;           // sd of the slope
  real<lower=0> interceptSigma;       // sd of the intercept
  real<lower=0> studySigma;           // sd of the study effect
  vector[K] assayIntercept;           // intercept of each assay 
  vector[K] assaySlopeLate;               // intercept change of each assay 
  vector[K] assaySlopeEarly;               // intercept change of each assay 
  vector[M] studyIntercept;           // intercept of each study:assay
}

transformed parameters{
  vector<lower=0, upper=1>[N] sensitivity;  // variable to contain % outcome
  sensitivity = inv_logit(assayIntercept[assay] + studyIntercept[study] +
    (assaySlopeEarly[assay]) .* timeVecEarly +
    (assaySlopeLate[assay]) .* timeVecLate);  // logistic regression model 
}

model {
  charSlope ~ normal(0, 100);
  meanSlopeLate ~ normal(0, 100);
  intercept ~ normal(0, 100);
  slopeSigmaEarly ~ gamma(4, 4);
  slopeSigmaLate ~ gamma(4, 4);
  interceptSigma ~ gamma(4, 4);
  studySigma ~ gamma(4, 4);
  assaySlopeEarly ~ normal(characteristics * charSlope, slopeSigmaLate);
  assaySlopeLate ~ normal(meanSlopeLate, slopeSigmaLate);
  assayIntercept ~ normal(intercept, interceptSigma);
  studyIntercept ~ normal(0, studySigma);
  nPositive ~ binomial(nTested, sensitivity);
}


