data {
  int<lower=0> H;
  int<lower=0> x_h[H];
  int<lower=0> n_h[H];
  int<lower=0> x_CC;
  int<lower=0> n_CC;
  int<lower=0> x_CT;
  int<lower=0> n_CT;
  real<lower=0> betascale;
  int<lower=0> nu;
}

parameters {
  real<lower=-1, upper=1> g;
  real<lower=0> tau;
  real eta_C;
  vector<lower=0>[H] lambda;
  vector[H] beta_raw;
}

transformed parameters{
  real<lower=0, upper=1> theta_C;
  real<lower=0, upper=1> theta_T;
  vector[H] eta_h;
  vector[H] beta;
  beta = tau * lambda .* beta_raw;
  theta_C = inv_logit(eta_C);
  theta_T = theta_C + g;
  eta_h = eta_C + beta;
}

model {
  target += normal_lpdf(eta_C|0,100);
  target += normal_lpdf(beta_raw|0,1);
  target += student_t_lpdf(tau|nu,0,betascale);
  target += student_t_lpdf(lambda|1,0,1);
  target += binomial_logit_lpmf(x_h|n_h,eta_h);
  target += binomial_lpmf(x_CC|n_CC,theta_C);
  target += normal_lpdf(g|0, 100);
  target += binomial_lpmf(x_CT|n_CT,theta_T);
}
