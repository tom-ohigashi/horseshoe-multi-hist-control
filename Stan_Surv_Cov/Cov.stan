data {
  int I_E;                                 // number of patients for event data
  int I_C;                                 // number of patients for censorring data
  int<lower=0,upper=1> TRT_E[I_E];         // treatment indicator for event data
  int<lower=0,upper=1> TRT_C[I_C];         // treatment indicator for censoring data
  real<lower=0>  Time_E[I_E];
  real<lower=0>  Time_C[I_C];
  int Ncov;
  matrix[I_E, Ncov] X_E;
  matrix[I_C, Ncov] X_C;
}

parameters {
  real eta_CC;
  real<lower=0> lambda_T;
  vector[Ncov] beta_cov;
}

transformed parameters{
  vector[I_E] hazard_E; // part of linear regression
  vector[I_C] hazard_C; // part of linear regression
  
  for(i in 1:I_E){
    hazard_E[i] = lambda_T * exp(eta_CC * TRT_E[i] + X_E[i]*beta_cov);
  }
  for(i in 1:I_C){
    hazard_C[i] = lambda_T * exp(eta_CC * TRT_C[i] + X_C[i]*beta_cov);
  }
}

model {
  target += normal_lpdf(beta_cov|0,10);
  target += normal_lpdf(eta_CC|0,10);
  target += gamma_lpdf(lambda_T|1, 1);
  target += exponential_lpdf(Time_E|hazard_E);
  target += exponential_lccdf(Time_C|hazard_C);
}
