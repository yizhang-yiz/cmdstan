data { 
  int N;
  vector[N] QoI_d;
} 

parameters {
  real<lower = 0> k[2];
} 

transformed parameters{
  real theta[2];
  real QoI[1];
  theta[1] = k[1];
  theta[2] = k[2];
  QoI = laplace_pde_forward_pde(theta);
}

model {
  k ~ normal(1.0, 0.01);
  for(i in 1:N) QoI_d[i] ~ normal(QoI[1], 0.01);
}
