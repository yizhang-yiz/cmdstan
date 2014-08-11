data {
  int<lower=0> N;
  int<lower=0> n;
  matrix[N,n] x;
  int<lower=0,upper=1> z[N];
}
parameters {
  vector[n] beta_raw;
  real<lower=0> tau;
}
transformed parameters {
  vector[n] beta;
  beta <- tau * beta_raw;
}
model {
  tau ~ cauchy(0, 1);
  beta_raw ~ normal(0, 1);
  z ~ bernoulli_logit(x*beta);
}