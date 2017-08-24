functions {
  // create an integer sequence
  int[] seq_int(int start, int end) {
    int N = end - start + 1;
    int seq[N];
    for(i in 1:N) seq[i] = i + start - 1;
    return(seq);
  }

  real[] twoCmtOral_ode(real t,
			real[] y,
			real[] theta,
			real[] x_r,
			int[] x_i) {
    real dydt[3];
    real ka;
    real k12;
    real k21;
    real k10;
    real dose = x_r[1];

    ka = theta[1];
    k10 = theta[2];
    k12 = theta[3];
    k21 = theta[4];

    ## we add a dose at t=24 within about 0.5 time units
    dydt[1] = -ka * y[1]; ## + dose * exp(normal_lpdf(t | 24, 0.5/4.));
    dydt[2] =  ka * y[1] - (k10 + k12) * y[2] + k21 * y[3];
    dydt[3] =  k12 * y[2] - k21 * y[3];
    return(dydt);
  }

  
  real[] mpi_function(real[] theta, real[] x_r, int[] x_i) {
    int T = x_i[1];
    real run[T,3] = integrate_ode_bdf(twoCmtOral_ode,
                                      theta[1:3], 
                                      0, x_r[2:(T+1)],
                                      theta[4:7],
                                      x_r[1:1], x_i,
                                      1E-7, 1E-7, 1000);

    return(run[:,2]);
  }

  real[] run_mpi_function(real[,] Theta, real[,] X_r, int[,] X_i);

  real[] integrate_oral_2cmt_serial(real[,] Theta, 
                                    int[] M,
                                    real[] time,
                                    real[,] x_r, int[,] x_i) {
    int J = size(M);
    real res[sum(M)];
    int cj;
    cj = 1;
    for(j in 1:J) {
      real run[M[j],3];
      run = integrate_ode_bdf(twoCmtOral_ode,
                              Theta[j,1:3],
                              0, time,
                              Theta[j,4:7],
                              x_r[j], x_i[j],
                              1E-7, 1E-7, 1000);

      for(m in 1:M[j])
        res[cj + m - 1] = run[m,2];

      cj = cj + M[j];
    }

    return(res);
  }
}
data {
  int<lower=1> T;
  int<lower=1> J;
  real<lower=0> theta[4];
  real<lower=0> theta_sd[4];
  real<lower=0> dose;
  int<lower=1> worker;
  real<lower=1> scale;
  int<lower=0,upper=1> use_mpi;
}
transformed data {
  real state0[J,3];
  real t0[J];
  real time[T*J];
  int M[J];
  real x_r[J,1+T];
  int x_i[J,1];
  vector[J] yobs_T;
  int<lower=0,upper=1> parallel = worker > 1 ? 1 : 0;
  real Theta_0[J,7];

  for(j in 1:J) {
    state0[j,1] = dose + j - 1;
    state0[j,2] = 0;
    state0[j,3] = 0;

    x_r[j,1] = dose + j - 1;
    x_r[j,2:T+1] = to_array_1d(seq_int(1, T));
    x_r[j,2:T+1] = to_array_1d(to_vector(x_r[j,2:T+1]) * scale);
    x_i[j,1] = T;

    t0[j] = 0;
    M[j] = T;

    Theta_0[j] = rep_array(1.0, 7);

    time[(j-1) * T + 1 : j * T] = to_array_1d(to_vector(seq_int(1, T)) * scale);
  }

  yobs_T = rep_vector(100., J);

  // obsolete as we distribute the data once
  //world_map = setup_mpi_function(Theta_0, x_r, x_i);

  if(parallel) {
    print("Parallel ODE integration using MPI.");
  } else {
    print("Serial ODE integration.");
  }
  
  if(use_mpi) {
    print("Using MPI.");
  } else {
    print("Not using MPI.");
  }

}
parameters {
  real<lower=0> theta_v[4];
  real<lower=0> dose0_v[J];
}
transformed parameters {
  
  //print("yhat = ", yhat);
  
  //reject("OK, we are good for now");
}
model {
  vector[J] conc_T;
  real yhat[sum(M)];
  real Theta[J,7];

  for(j in 1:J) {
    Theta[j,1] = dose0_v[j];
    Theta[j,2] = 0;
    Theta[j,3] = 0;
    Theta[j,4:7] = theta_v;
  }
  
  if(use_mpi) {
    yhat = run_mpi_function(Theta, x_r, x_i);
  } else {
    yhat = integrate_oral_2cmt_serial(Theta, M, time[1:T], x_r, x_i);
  }

  for(j in 1:J)
    conc_T[j] = yhat[(j-1) * T + T];
  //print("conc_T = ", conc_T);

  theta_v ~ lognormal(log(theta), theta_sd);

  dose0_v ~ lognormal(log(10.), 0.2);

  target += normal_lpdf(yobs_T | conc_T, 5);
}
generated quantities {
  /*
  real single_res[sum(M),3];
  real parallel_res[sum(M),3];

  single_res = integrate_oral_2cmt_serial(state0, t0, M, time[1:T], theta_v, x_r, x_i);

  parallel_res = integrate_oral_2cmt_parallel(state0, t0, M,
  time[1:T], theta_v, x_r, x_i);
  */
}
