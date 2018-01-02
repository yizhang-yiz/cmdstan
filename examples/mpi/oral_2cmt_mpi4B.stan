functions {
  // count number times elem appears in test set
  int count_relem(real[] test, real elem) {
    int count = 0;
    for(i in 1:num_elements(test))
      if(test[i] == elem)
        count = count + 1;
    return(count);
  }

  // find elements in test which are equal to elem
  int[] which_relem(real[] test, real elem) {
    int res[count_relem(test, elem)];
    int ci = 1;
    for(i in 1:num_elements(test))
      if(test[i] == elem) {
        res[ci] = i;
        ci = ci + 1;
      }
    return(res);
  }
  
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
   * NOTE: A & B are not scaled by any volume!
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
    if(2*lkSum < log(4.0) + lk[1] + lk[3]) {
      reject("System discriminant must be real!");
    }

    // first and second rate constant are roots for pq-quadratic form
    // with p = -kSum, q = k21 * k = k[3] * k[1]
    
    // log of first rate constant alpha
    mm[2] = log(0.5) + log_sum_exp(lkSum, 0.5 * log_diff_exp(2*lkSum, log(4.0) + lk[1] + lk[3]) );

    // log of second rate constant beta (Vieta)
    mm[3] = lk[1] + lk[3] - mm[2];

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

  real[] twoCmtOral_analytic(real ldose, real[] Dt, real[] theta) {
    //            (lka, lCL, lV1, lQ, lV2)
    // theta = log(c(ka, CL, CLr, V, Vr))
    real lV1 = theta[4] + theta[5];
    vector[5] lmicro = trans_oral2cmt_macro2micro(theta[1], // ka
                                                  theta[2] + theta[3], // CL
                                                  lV1, // V1
                                                  theta[2] - theta[3], // Q1
                                                  theta[4] - theta[5]); // V2
    // lmicro: lka, log(alpha), log(beta), lA, lB
    vector[3] state0 = [ ldose, -25, -25 ]';
    matrix[num_elements(Dt),3] lstate = pk_oral_2cmt(state0, to_vector(Dt), lmicro[1], lmicro[2], lmicro[3], lmicro[4], lmicro[5]);
    return(to_array_1d(lstate[:,2] - lV1));
  }

  real[] twoCmtOral_ode(real t,
			real[] y,
			real[] theta,
			real[] x_r,
			int[] x_i) {
    real dydt[3];
    real ka = theta[1];
    real CL = theta[2];
    real V1 = theta[3];
    real Q  = theta[4];
    real V2 = theta[5];
    
    dydt[1] = -ka * y[1];
    dydt[2] =  ka * y[1] - (CL + Q) * y[2]/V1 + Q * y[3]/V2;
    dydt[3] =  Q * y[2]/V1 - Q * y[3]/V2;
    return(dydt);
  }

  
  vector mpi_function(vector mu, vector eta, real[] x_r, int[] x_i) {
    int T = x_i[1];
    int use_ode = x_i[2];
    int num_omega = x_i[3];
    int ind_omega[num_omega] = x_i[4:(4+num_omega-1)];
    int known_sigma_y = x_i[3+num_omega+1];
    int known_omega = x_i[3+num_omega+2];
    real run[T];
    real ldose = x_r[2*T+1];
    real sigma_y = mu[1];
    vector[6] theta = mu[2:7];
    vector[num_omega] omega = mu[8:(8+num_omega-1)];
    //vector[num_omega] true_omega = x_r[2*T+2:(2*T+2+num_omega-1)];
    real eta_a[6] = to_array_1d(theta);
    vector[1] res;

    // take over the parameters which are per subject
    eta_a[ind_omega] = to_array_1d(eta);
    
    if(use_ode) {
      real state0[3] = { exp(ldose + eta_a[1]), 0, 0 } ;
      real lV1 = eta_a[1+4] + eta_a[1+5];
      real theta_macro[5] = {
        exp(eta_a[1+1]), // ka
        exp(eta_a[1+2] + eta_a[1+3]), // CL
        exp(lV1), // V1
        exp(eta_a[1+2] - eta_a[1+3]), // Q
        exp(eta_a[1+4] - eta_a[1+5])  // V2
      };
      real run2[T,3] = integrate_ode_rk45(twoCmtOral_ode,
                                          state0, 
                                          0, x_r[1:T],
                                          theta_macro,
                                          x_r[1:0], x_i,
                                          1E-5, 1E-7, 1000);
      run = log(run2[:,2]);
      for(i in 1:T)
        run[i] = run[i] - lV1;
    } else {
      run = twoCmtOral_analytic(ldose + eta_a[1], x_r[1:T], eta_a[2:6]);
    }

    /*
    for(i in 1:T) {
      if(is_nan(run[i])) {
        print("Found NaN, i = ", i, "; Dt = ", x_r[i]);
        print("eta_a = ", eta_a);
      }
    }
    */

    res[1] = normal_lpdf(x_r[(T+1):2*T] | run, known_sigma_y ? 0.05 : sigma_y )
             + normal_lpdf( eta | theta[ind_omega], known_omega ? to_vector(x_r[2*T+2:(2*T+2+num_omega-1)]) : omega );

    return( res );
  }

  /**/
  vector map_rect_mpi(vector eta, vector[] Theta, real[,] X_r, int[,] X_i);
  vector map_rect_serial(vector eta, vector[] Theta, real[,] X_r, int[,] X_i);
  /**/
  /*
  vector map_rect_mpi(vector eta, vector[] Theta, real[,] X_r, int[,] X_i) {
    vector[1] res;
    return(res);
  }
  vector map_rect_serial(vector eta, vector[] Theta, real[,] X_r, int[,] X_i) {
    vector[1] res;
    return(res);
  }
  */
  
  vector map_rect_stan(vector eta,
                       vector[] Theta, 
                       int[] M,
                       real[,] x_r, int[,] x_i) {
    int J = size(M);
    vector[sum(M)] res;
    int cj = 1;
    for(j in 1:J) { 
      vector[M[j]] run = mpi_function(eta, Theta[j], x_r[j], x_i[j]);
      
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
  real<lower=0> time[T];
  real true_theta[6];
  real<lower=0> true_omega[6];
  int<lower=0,upper=2> use_map_rect;
  int<lower=0,upper=1> use_ode;
  real<lower=0> dose;
  int<lower=0,upper=1> known_sigma_y;
  int<lower=0,upper=1> known_omega;
}
transformed data {
  real t0[J];
  int M[J];
  real ldose = log(dose);
  int num_omega = 6 - count_relem(true_omega, 0);
  int x_i[J,3+num_omega+2];
  real x_r[J,2*T+1+num_omega];
  int ind_omega[num_omega];
  matrix[J,6] Eta_sim;
  row_vector[6] eta_mean = rep_row_vector(0, 6);
  cholesky_factor_cov[num_omega] L_Omega_known;
  int ind_active_theta[5] = {2, 3, 4, 5, 6};
  vector[5] true_theta_macro = [
    exp(true_theta[1+1]), // ka
    exp(true_theta[1+2] + true_theta[1+3]), // CL
    exp(true_theta[1+4] + true_theta[1+5]), // V1
    exp(true_theta[1+2] - true_theta[1+3]), // Q
    exp(true_theta[1+4] - true_theta[1+5])  // V2
  ]';

  {
    int ci = 1;
    for(i in 1:6) {
      if(true_omega[i] > 0) {
        ind_omega[ci] = i;
        ci = ci + 1;
      }
    }
  }

  L_Omega_known = cholesky_decompose(diag_matrix(square(to_vector(true_omega[ind_omega]))));

  // we fix the population frel to be 1
  if(true_theta[1] != 0)
    reject("Population frel must be 1!");

  print("True population parameters ", round(to_vector(true_theta)*100)/100);
  print("True population parameters (ka, CL, V1, Q, V2) ", round(100*true_theta_macro)/100);
  print("Number of random effects ", num_omega);
  print("Random effects on parameters ", ind_omega);
  print("Random effects omegas ", true_omega[ind_omega]);

  if(known_omega) {
    print("Assuming known omegas.");
  } else {
    print("Assuming unknown omegas.");
  }

  if(known_sigma_y) {
    print("Assuming known sigma_y.");
  } else {
    print("Assuming unknown sigma_y.");
  }

  // ensure that the mean over the individualized parameters is
  // exactly equal the true mean value
  for(j in 1:J) {
    Eta_sim[j] = to_row_vector(true_theta);
    for(i in 1:6) {
      if(true_omega[i] > 0)
        Eta_sim[j,i] = normal_rng(true_theta[i], true_omega[i]);
      eta_mean[i] = eta_mean[i] + Eta_sim[j,i];
    }
  }

  eta_mean = eta_mean / J;
  
  for(j in 1:J) {
    real theta_j[6] = to_array_1d(Eta_sim[j] - eta_mean + to_row_vector(true_theta));
    real ldose_j = normal_rng(ldose, log(4)/1.96);

    x_r[j,1:T] = time;
    x_r[j,(T+1):2*T] = twoCmtOral_analytic(ldose_j + theta_j[1], time, theta_j[2:6]);
    x_r[j,2*T+1] = ldose_j;

    for(k in 1:T)
      x_r[j,T+k] = x_r[j,T+k] + normal_rng(0, 0.05);

    x_r[j,(2*T+2):(2*T+2+num_omega-1)] = true_omega[ind_omega];

    x_i[j,1] = T;
    x_i[j,2] = use_ode;
    x_i[j,3] = num_omega;
    x_i[j,4:(4+num_omega-1)] = ind_omega;
    x_i[j,3+num_omega+1] = known_sigma_y;
    x_i[j,3+num_omega+2] = known_omega;

    t0[j] = 0;
    M[j] = 1;
  }
 
  print("Problem size J = ", J);
  
  if(use_map_rect == 0) {
    print("Using map_rect_mpi.");
  } else if(use_map_rect == 1) {
    print("Using map_rect_serial.");
  } else if(use_map_rect == 2) {
    print("Using map_rect_stan.");
  }
  if(use_ode) {
    print("Using ODE integration.");
  } else {
    print("Using analytic solution.");
  }
}
parameters {
  vector[num_elements(ind_active_theta)] theta_raw;
  vector[num_omega] Eta_v[J];
  //vector[J] Eta_v[num_omega];
  vector<lower=0>[num_omega] omega_raw;
  real<lower=0> sigma_y_raw;
}
transformed parameters {
  vector[6] theta_v = rep_vector(0, 6);
  vector<lower=0>[num_omega] omega_v;
  real<lower=0> sigma_y;
  // we shift by the true theta's such that init=0 will start with the
  // true values.
  theta_v[ind_active_theta] = to_vector(true_theta[ind_active_theta]) + theta_raw;
  //theta_v[5] = 0;
  //theta_v = to_vector(true_theta) + theta_raw;
  omega_v = known_omega ? to_vector(true_omega[ind_omega]) : omega_raw;
  sigma_y = known_sigma_y ? 0.05 : sigma_y_raw;
}
model {
  vector[1+6+num_omega] mu;
  //vector[num_omega] Eta_a[J];

  mu[1] = sigma_y;
  mu[2:7] = theta_v;
  mu[8:(8+num_omega-1)] = omega_v;

  //theta_v ~ normal(to_vector(true_theta), 1);
  theta_raw ~ normal(0, 1);
  
  // CP parametrization
  /*
  if(known_omega) {
    Eta_v ~ multi_normal_cholesky(theta_v[ind_omega], L_Omega_known);
  } else {
    Eta_v ~ multi_normal_cholesky(theta_v[ind_omega], cholesky_decompose(diag_matrix(square(omega_v))));
  }
  */

  /*
  if(known_omega) {
    for(i in 1:num_omega)
      Eta_v[i] ~ normal(theta_v[ind_omega[i]], true_omega[ind_omega[i]]);
  } else {
    for(i in 1:num_omega)
      Eta_v[i] ~ normal(theta_v[ind_omega[i]], omega_v[i]);
  }
  */
  
  omega_raw ~ normal(0, 1);
  sigma_y_raw ~ normal(0, 0.25);

  /*
  for(j in 1:J) {
    for(i in 1:num_omega)
      Eta_a[j,i] = Eta_v[i,j];
  }
  */

  if(use_map_rect == 0) {
    target += map_rect_mpi(mu, Eta_v, x_r, x_i);
  } else if(use_map_rect == 1) {
    target += map_rect_serial(mu, Eta_v, x_r, x_i);
  } else {
    target += map_rect_stan(mu, Eta_v, M, x_r, x_i);
  }
}
generated quantities {
  vector[num_elements(ind_active_theta)] bias_theta = theta_v[ind_active_theta] - to_vector(true_theta[ind_active_theta]);
  vector[num_elements(ind_active_theta)] mse_theta = square(bias_theta);
  vector[num_omega] bias_omega = omega_v - to_vector(true_omega[ind_omega]);
  vector[num_omega] mse_omega = square(bias_omega);
  vector[1] bias_sigma_y = [ sigma_y - 0.05 ]';  
}
