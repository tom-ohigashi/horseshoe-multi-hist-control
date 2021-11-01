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
  real g;
  real<lower=0, upper=1> delta[H];
  real<lower=0> a;
  real<lower=0> b;
}

transformed parameters{
  real<lower=0, upper=1> theta_T;
  real<lower=0, upper=1> mu;
  real<lower=0, upper=1> sigmasq;
  real lik_C;
  real lik_T;
  real SC;
  real delta_sum1[H];
  real delta_sum2[H];
  real delta_sum3[H];
  real sum_log_delta[H];

  theta_T = theta_C + g;
  mu = a / (a + b);
  sigmasq = (mu * (1 - mu)) / (a + b + 1);
  
  for(h in 1:H){
    delta_sum1[h] = delta[h] * x_h[h];
    delta_sum2[h] = delta[h] * (n_h[h] - x_h[h]);
    delta_sum3[h] = delta[h] * n_h[h];
    sum_log_delta[h] = log_sum_exp(log(0.9) + beta_lpdf(delta[h]|a,b), log(0.1) + (normal_lpdf(delta[h]|0,sqrt(sigmasq / 6.25)) - 1*normal_lccdf(0|0,sqrt(sigmasq / 6.25))));
  }
  
  lik_C = (sum(delta_sum1) + x_CC) * log(theta_C) + (sum(delta_sum2) + (n_CC - x_CC)) * log(1 - theta_C);
  lik_T = x_CT * log(theta_T) + (n_CT - x_CT) * log(1 - theta_T);
  SC = lgamma(sum(delta_sum1) + 1) + lgamma(sum(delta_sum2) + 1) - lgamma(sum(delta_sum3) + 1 + 1);
}

model {
  target += beta_lpdf(theta_C|1,1);
  target += normal_lpdf(g|0,100);
  target += uniform_lpdf(mu|0,1);
  target += gamma_lpdf(1 / sigmasq|0.01,0.01);
  target += -log(sigmasq ^ 2);
  target += sum(sum_log_delta);
  target += lik_C;
  target += lik_T;
  target += -SC;
}
