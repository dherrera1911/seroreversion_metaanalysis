data {
  int<lower=0> N;                           // number of observations
  int<lower=0> K;                           // number of assays
  int<lower=0> M;                           // number of studies:assays
  int<lower=0> C;                           // number of characteristics
  int<lower=1, upper=K> assay[N];           // assay index vector
  int<lower=1, upper=M> study[N];           // study:assay index vector
  matrix[K,C] characteristics;              // assay characteristics that we analyze
  int<lower=0> nNegative[N];             // Number of positive samples
  int<lower=0> nTested[N];             // Number of samples tested
}

parameters {
  // 
  vector[C] charIntercept;            // intercept effets of assay characteristics
  real<lower=0> interceptSigma;       // sd of the intercept among assays of same class
  real<lower=0> studySigma;           // sd of the study effect
  vector[K] assayIntercept;           // intercept of each assay 
  vector[M] studyIntercept;           // intercept of each study:assay
}

transformed parameters{
  vector<lower=0, upper=1>[N] specificity;  // variable to contain % outcome
  specificity = inv_logit(assayIntercept[assay] + studyIntercept[study]);  // logistic regression model 
}

model {
  charIntercept ~ normal(0, 100);
  interceptSigma ~ gamma(4, 4);
  studySigma ~ gamma(4, 4);
  assayIntercept ~ normal(characteristics * charIntercept, interceptSigma);
  studyIntercept ~ normal(0, studySigma);
  nNegative ~ binomial(nTested, specificity);
}


