data {
  int<lower=0> H;
  int<lower=0> x_h[H];
  int<lower=0> n_h[H];
  int<lower=0> x_CC;
  int<lower=0> n_CC;
  int<lower=0> x_CT;
  int<lower=0> n_CT;
}

parameters {
  real<lower=-1, upper=1> g;
  real mu;
  real<lower=0> tau;
  vector[H] eta_raw;
  real eta_C;
}

transformed parameters{
  real<lower=0, upper=1> theta_T;
  real<lower=0, upper=1> theta_C;
  vector<lower=0, upper=1>[H] theta_h;
  vector[H] eta_h;

  for(h in 1:H){
    eta_h[h] = mu + tau*eta_raw[h];
  }
  theta_h = inv_logit(eta_h);
  theta_C = inv_logit(eta_C);
  theta_T = theta_C + g;
}

model {
  target += normal_lpdf(mu|0, 100);
  target += normal_lpdf(tau|0, 1);
  target += normal_lpdf(eta_raw|0,1);
  target += log_sum_exp(log(0.9) + normal_lpdf(eta_C|mu,tau), log(0.1) + normal_lpdf(eta_C|mu,100));
  target += normal_lpdf(g|0, 100);
  target += binomial_lpmf(x_h|n_h,theta_h);
  target += binomial_lpmf(x_CC|n_CC,theta_C);
  target += binomial_lpmf(x_CT|n_CT,theta_T);
}
