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
  real<lower=0, upper=1> theta_C;
  real<lower=-1, upper=1> g;
}

transformed parameters{
  real<lower=0, upper=1> theta_T;
  theta_T = theta_C + g;
}

model {
  target += normal_lpdf(g|0, 100);
  target += beta_lpdf(theta_C|1,1);
  target += binomial_lpmf(x_h|n_h,theta_C);
  target += binomial_lpmf(x_CC|n_CC,theta_C);
  target += binomial_lpmf(x_CT|n_CT,theta_T);
}

