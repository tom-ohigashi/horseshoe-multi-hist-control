data {
  int I_E;                                 // number of patients for event data
  int I_C;                                 // number of patients for censorring data
  int<lower=0,upper=1> TRT_E[I_E];         // treatment indicator for event data
  int<lower=0,upper=1> TRT_C[I_C];         // treatment indicator for censoring data
  int SID_E[I_E];                          // 
  int SID_C[I_C];
  real<lower=0>  Time_E[I_E];
  real<lower=0>  Time_C[I_C];
  int H;                                   // Numver of Historical study
  real<lower=0> betascale;
  int<lower=0> nu;
  int Ncov;
  matrix[I_E, Ncov] X_E;
  matrix[I_C, Ncov] X_C;
}

parameters {
  real eta_CC;
  real<lower=0> lambda_T;
  vector[H] beta_raw;
  real<lower=0> tau;
  vector<lower=0>[H] lambda;
  vector[Ncov] beta_cov;
}

transformed parameters{
  vector[I_E] hazard_E; // part of linear regression
  vector[I_C] hazard_C; // part of linear regression
  
  vector[H] eta_h;
  vector[H] beta;

  beta = tau * lambda .* beta_raw;
  eta_h = eta_CC + beta;

  for(i in 1:I_E){
    if(SID_E[i] <= H){
      hazard_E[i] = lambda_T * exp(eta_h[SID_E[i]] + X_E[i]*beta_cov);
    }else{
      hazard_E[i] = lambda_T * exp(eta_CC * TRT_E[i] + X_E[i]*beta_cov);
    }
  }
    for(i in 1:I_C){
    if(SID_C[i] <= H){
      hazard_C[i] = lambda_T * exp(eta_h[SID_C[i]] + X_C[i]*beta_cov);
    }else{
      hazard_C[i] = lambda_T * exp(eta_CC * TRT_C[i] + X_C[i]*beta_cov);
    }
  }
}

model {
  target += normal_lpdf(beta_cov|0,10);
  target += normal_lpdf(eta_CC|0,10);
  target += normal_lpdf(beta_raw|0,1);
  target += student_t_lpdf(tau|nu,0,betascale);
  target += student_t_lpdf(lambda|1,0,1);
  target += gamma_lpdf(lambda_T|1, 1);
  target += exponential_lpdf(Time_E|hazard_E);
  target += exponential_lccdf(Time_C|hazard_C);
}
