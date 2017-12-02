functions {
  // create an integer sequence
  int[] seq_int(int start, int end) {
    int N = end - start + 1;
    int seq[N];
    for(i in 1:N) seq[i] = i + start - 1;
    return(seq);
  }

  real log_diff_exp_abs(real la, real lb) {
    return(0.5 * log_diff_exp(log_sum_exp(2*la, 2*lb), log(2) + la + lb));
  }

  /*
   * transform from log of lka, CL, V1, Q, V2 to log of lka, alpha,
   * beta, A, B
   */
  vector trans_oral2cmt_macro2micro(real lka, real lCL, real lV1, real lQ, real lV2) {
    vector[3] lk;
    real lkSum;
    vector[5] mm;
    // first project onto "k" parametrization and then onto micro
    // and macro constants

    // lk
    lk[1] = lCL - lV1;
    // lk12
    lk[2] = lQ - lV1;
    //lk21
    lk[3] = lQ - lV2;

    // log(k+k12+k21)
    lkSum = log_sum_exp(lk);

    mm[1] = lka;

    // check that discriminat is for all patients real, i.e. that 
    // (k10 + k12 + k21)^2 - 4 k10 k21 > 0
    // otherwise the eigenvalues would be imaginary leading to oscillations
    if(2*lkSum < log(4.0) + lk[1] + lk[3])
      reject("System discriminant must be real!");

    // log of second rate constant beta
    mm[3] = log(0.5) + log_diff_exp(lkSum, 0.5 * log_diff_exp(2*lkSum, log(4.0) + lk[1] + lk[3]) );

    // log of first rate constant alpha
    mm[2] = lk[1] + lk[3] - mm[3];

    // macro constants
    mm[4] = log_diff_exp_abs(lk[3], mm[2]) - log_diff_exp_abs(mm[2], mm[3]);
    mm[5] = log_diff_exp_abs(lk[3], mm[3]) - log_diff_exp_abs(mm[3], mm[2]);

    return(mm);
  }

  matrix pk_oral_2cmt(vector state0, vector Dt,
                      real lka, real lalphaR, real lbetaR, real lA, real lB) {
    real lstateRefOral; // ref state for the 2-cmt with oral cmt (only the oral cmt)
    real lstateRef[2];  // ref state for the 2-cmt without oral cmt
    int N;
    real alphaR;
    real betaR;
    real ka;
    real lAt;
    real lBt;
    matrix[num_elements(Dt),3] lstate;
    real lk12;
    real lk21;
    real lD;
    real lad2;
    real lbd2;
    real ltemp;

    N = num_elements(Dt);
    
    ka = exp(lka);
    alphaR = exp(lalphaR);
    betaR  = exp(lbetaR);
    
    // Bateman coefficients
    lAt = lA + lka - log_diff_exp_abs(lka, lalphaR);
    lBt = lB + lka - log_diff_exp_abs(lka, lbetaR );
    
    // needed constant for the unobserved peripheral cmt C
    lD = log_diff_exp(lalphaR, lbetaR);   // discriminant which is always positive
    ltemp = log_sum_exp(lB + lD, lbetaR);
    lk12 = log_diff_exp(log_sum_exp(lalphaR, lbetaR), log_sum_exp(2*ltemp, lalphaR + lbetaR) - ltemp );
    lk21 = log_diff_exp(lalphaR, lA + lD);

    lad2 = 2 * log_diff_exp_abs(lalphaR, lka);
    lbd2 = 2 * log_diff_exp_abs(lbetaR , lka);
    
    // by convention time starts just at the first observation
    lstateRefOral = state0[1];
    lstateRef[1]  = state0[2];
    lstateRef[2]  = state0[3];
    for(i in 1:N) {
      lstate[i,1] = lstateRefOral - ka * Dt[i];
      // solution for the concentration which is in the central and
      // peripheral cmt
      lstate[i,2] = lstateRef[1] + log_sum_exp(lA - alphaR * Dt[i], lB - betaR * Dt[i]);
      lstate[i,3] = lstateRef[2] + log_sum_exp(lB - alphaR * Dt[i], lA - betaR * Dt[i]);

      // other changes in the state can only meaningful be calculated
    // if Dt[i] is large enough to allow for diffusion to occur
      if(Dt[i] > 0.) {
        lstate[i,2] = log_sum_exp(lstate[i,2], lstateRef[2] + lD - lk12 + lA + lB + log_diff_exp(- betaR * Dt[i], - alphaR * Dt[i]) );
        lstate[i,3] = log_sum_exp(lstate[i,3], lstateRef[1] + lk12 - lD + log_diff_exp(-betaR * Dt[i], -alphaR * Dt[i]));


        // add in the part which stems from oral cmt which results in
        // the superposition of Bateman functions in the main cmt
        lstate[i,2] = log_sum_exp(lstate[i,2], lstateRefOral + log_sum_exp(lAt + log_diff_exp_abs( -alphaR * Dt[i], -ka * Dt[i]),
                                                                           lBt + log_diff_exp_abs( -betaR  * Dt[i], -ka * Dt[i]))
                                  );
        // last, we add into the peripheral cmt the effect of the oral
        // cmt dosing
        //lstate[i,3] = log_sum_exp(lstate[i,3], lstateRefOral + lk12 + lka - lD - ka * Dt[i] + log( D*A2i*B2i + A2i * exp(- (alphaR-ka) * Dt[i]) - B2i * exp(-(betaR-ka)*Dt[i])   ) );
        // the huge expression below is a sign-sorted version of (brute force)
        // k12 * ka /[ (alpha-ka) * (beta-ka) ] * [ exp(-ka * t) - (ka-beta)/(alpha-beta) * exp(-alpha*t) + (ka-alpha)/(alpha-beta) * exp(-beta*t) ]
        lstate[i,3] = log_sum_exp(lstate[i,3], lstateRefOral + lk12 + lka - lD - lad2 - lbd2 +
                                  log_diff_exp(log_sum_exp(log_sum_exp(lD + log_sum_exp(lalphaR + lbetaR, 2*lka) - ka * Dt[i], lalphaR + lbd2 - alphaR * Dt[i]), lka    + lad2 - betaR * Dt[i] ),
                                               log_sum_exp(log_sum_exp(lD + lka + log_sum_exp(lalphaR,   lbetaR) - ka * Dt[i], lka     + lbd2 - alphaR * Dt[i]), lbetaR + lad2 - betaR * Dt[i] )
                                               )
                                  );
      }
    }
    
    return lstate;
  }

  real[] twoCmtOral_analytic(real dose, real[] Dt, real[] theta) {
    vector[5] lmicro = trans_oral2cmt_macro2micro(log(theta[1]), log(theta[2]), 0, log(theta[3]),
                                                  log(theta[3]) - log(theta[4]));
    vector[3] state0;
    matrix[num_elements(Dt),3] lstate;
    state0[1] = log(dose); state0[2] = -35; state0[3] = -35;
    lstate = pk_oral_2cmt(state0, to_vector(Dt), lmicro[1], lmicro[2], lmicro[3], lmicro[4], lmicro[5]);
    return(to_array_1d(exp(lstate[:,2])));
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

    ka = theta[1];
    k10 = theta[2];
    k12 = theta[3];
    k21 = theta[4];

    dydt[1] = -ka * y[1];
    dydt[2] =  ka * y[1] - (k10 + k12) * y[2] + k21 * y[3];
    dydt[3] =  k12 * y[2] - k21 * y[3];
    return(dydt);
  }

  
  vector mpi_function(vector eta, vector theta, real[] x_r, int[] x_i) {
    int T = x_i[1];
    int use_ode = x_i[2];
    int use_shared = x_i[3];
    vector[1] run;
    real state0[3];
    real theta_run[4];
    
    state0[1] = theta[1];
    state0[2] = 0;
    state0[3] = 0;

    if(use_shared) { 
      theta_run = to_array_1d(eta);
    } else {
      theta_run = to_array_1d(theta[2:5]);
    }

    if(use_ode) {
      real run2[T,3] = integrate_ode_bdf(twoCmtOral_ode,
                                         state0, 
                                         0, x_r[2:(T+1)],
                                         theta_run,
                                         x_r[1:0], x_i,
                                         1E-5, 1E-7, 1000);
      run[1] = run2[T,2];
    } else {
      real run2[T] = twoCmtOral_analytic(theta[1], x_r[2:(T+1)], theta_run);
      run[1] = run2[T];
    }
 
    return(run);
  }

  vector map_rect(vector eta, vector[] Theta, real[,] X_r, int[,] X_i);

  vector map_rect_serial(vector eta,
                         vector[] Theta, 
                         int[] M,
                         real[,] x_r, int[,] x_i) {
    int J = size(M);
    vector[sum(M)] res;
    int cj;
    cj = 1;
    for(j in 1:J) { 
      vector[M[j]] run;

      run = mpi_function(eta, Theta[j], x_r[j], x_i[j]);
      
      for(m in 1:M[j])
        res[cj + m - 1] = run[m];

      cj = cj + M[j];
    }

    return(res);
  }
}
data {
  int<lower=1> T;
  int<lower=1> J;
  real<lower=0> eta[4];
  real<lower=0> eta_sd[4];
  real<lower=0> dose;
  real<lower=1> scale;
  int<lower=0,upper=1> use_mpi;
  int<lower=0,upper=1> use_ode;
  int<lower=0,upper=1> use_shared;
}
transformed data {
  real t0[J];
  real time[T*J];
  int M[J];
  real x_r[J,1+T];
  int x_i[J,3];
  vector[J] yobs_T;

  yobs_T = rep_vector(100., J);

  for(j in 1:J) {
    yobs_T[j] = 1 + dose + j - 1;

    x_r[j,1] = yobs_T[j];
    x_r[j,2:T+1] = to_array_1d(seq_int(1, T));
    x_r[j,2:T+1] = to_array_1d(to_vector(x_r[j,2:T+1]) * scale);

    x_i[j,1] = T;
    x_i[j,2] = use_ode;
    x_i[j,3] = use_shared;

    t0[j] = 0;
    M[j] = 1;

    time[(j-1) * T + 1 : j * T] = to_array_1d(to_vector(seq_int(1, T)) * scale);
  }

  // obsolete as we distribute the data once
  //world_map = setup_mpi_function(Theta_0, x_r, x_i);

  print("Problem size J = ", J);
  
  if(use_mpi) {
    print("Using MPI.");
  } else {
    print("Not using MPI.");
  }
  if(use_ode) {
    print("Using ODE integration.");
  } else {
    print("Using analytic solution.");
  }
  if(use_shared) {
    print("Using shared parameters.");
  } else {
    print("Not using shared parameters.");
  }
}
parameters {
  vector[4] log_eta_v;
  vector[J] log_dose0_v;
}
transformed parameters {
  vector[4] eta_v;
  vector[J] dose0_v;
  eta_v = exp(log_eta_v);
  dose0_v = exp(log_dose0_v);
  //print("yhat = ", yhat);
  
  //reject("OK, we are good for now");
}
model {

  if(use_shared) {
    vector[1] Theta[J];
    vector[J] conc_T;

    for(j in 1:J) {
      Theta[j,1] = dose0_v[j];
    }

    if(use_mpi) {
      conc_T = map_rect(eta_v, Theta, x_r, x_i);
    } else {
      conc_T = map_rect_serial(eta_v, Theta, M, x_r, x_i);
    }
    target += normal_lpdf(yobs_T | conc_T, 5);     
  } else {
    vector[0] eta0;
    vector[5] Theta[J];
    vector[J] conc_T;

    for(j in 1:J) {
      Theta[j,1] = dose0_v[j];
      Theta[j,2:5] = eta_v;
    }
  
    if(use_mpi) {
      conc_T = map_rect(eta0, Theta, x_r, x_i);
    } else {
      conc_T = map_rect_serial(eta0, Theta, M, x_r, x_i);
    }
    target += normal_lpdf(yobs_T | conc_T, 5);     
  }

  //for(j in 1:J)
  //conc_T[j] = yhat[(j-1) * T + T];
  //print("conc_T = ", conc_T);

  log_eta_v ~ normal(log(eta), eta_sd);

  log_dose0_v ~ normal(log(10.), 0.2);


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
