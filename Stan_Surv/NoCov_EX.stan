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
  real Tau;
}

parameters {
  real<lower=0> lambda_T;
  real mu;
  real<lower=0> tau;
  vector[H+1] eta_raw;
}

transformed parameters{
  vector[I_E] hazard_E; // part of linear regression
  vector[I_C] hazard_C; // part of linear regression
  vector[H] eta_h;
  real eta_CC;

  for(h in 1:H){
    eta_h[h] = mu + tau * eta_raw[h];
  }
  eta_CC = mu + tau * eta_raw[H+1];

  for(i in 1:I_E){
    if(SID_E[i] <= H){
      hazard_E[i] = lambda_T * exp(eta_h[SID_E[i]]);
    }else{
      hazard_E[i] = lambda_T * exp(eta_CC * TRT_E[i]);
    }
  }
    for(i in 1:I_C){
    if(SID_C[i] <= H){
      hazard_C[i] = lambda_T * exp(eta_h[SID_C[i]]);
    }else{
      hazard_C[i] = lambda_T * exp(eta_CC * TRT_C[i]);
    }
  }
}

model {
  target += normal_lpdf(mu|0, 10);
  target += normal_lpdf(tau|0,Tau);
  target += normal_lpdf(eta_raw|0,1);
  target += gamma_lpdf(lambda_T|0.01, 0.01);
  target += exponential_lpdf(Time_E|hazard_E);
  target += exponential_lccdf(Time_C|hazard_C);
}
