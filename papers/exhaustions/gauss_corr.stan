transformed data {
  vector[100] mu;
  real rho;
  cov_matrix[100] Sigma;
  cholesky_factor_cov[100] L_Sigma;

  mu <- rep_vector(0, 100);

  //vector[100] e;

  rho <- 0.9;
  
  for (i in 1:100)
    for (j in 1:100)
      Sigma[i, j] <- pow(rho, abs(i - j));

  L_Sigma <- cholesky_decompose(Sigma);

  /*
  e <- eigenvalues_sym(Sigma);

  for (n in 1:50)
    print(2 * pi() *  sqrt(e[n]));
  */

}

parameters {
  vector[100] x;
}

model {
  x ~ multi_normal_cholesky(mu, L_Sigma);
}
