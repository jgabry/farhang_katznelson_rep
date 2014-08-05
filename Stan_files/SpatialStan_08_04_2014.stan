data {
// Dimensions of the data
  int<lower=1>            N ; # number of obs in data
  int<lower=0>            R ; # number of groups (i.e. region/period combos)
  int<lower=0>            J ; # number of betas (number of vars in X plus 1 for intercept)
  int<lower=0>            C ; # number of deltas (number of vars in Z)

// Outcome variable 
  int<lower=0,upper=1>    Y[N] ; # indicator for prolabor vote
  
// Predictors w/ effects assumed constant over regional/temporal dims ##
  matrix[N, C]            Z ;  # dem party & labor committee indicators 
  
// Predictors to get smoothing priors
  matrix[N, R]		        X[J-1] ;  # urbanpct, aapct, unionpop
  
// Stuff for spatial-smoothing priors
  matrix[N, R]		        M ;         # incidence matrix
  matrix[R, R] 		        A ;         # adjacency matrix
  vector[R]               d ;         # diagonal vector of D (degree matrix)
  int                     A_N ;       # number of adjacent region pairs
  int                     A1[A_N] ;   # first half of adjacency pairs
  int                     A2[A_N] ;   # second half of adjacency pairs
}

transformed data {
  real<lower=0>           sigma_delta ; # sd for priors on delta
  real<lower=0>           dispersion ;  # to be used making prior for rho
  vector[C]               mu_delta ;    # mean for priors on delta
  real                    mu_beta ;     # prior mean for betas 1 through 3
  vector[R]               sqrt_d ;      # diagonal elements of degree matrix

  dispersion  <- 1 ;
  sigma_delta <- 100.0 ;
  mu_beta     <- 0.0 ;
  mu_delta    <- rep_vector(0.0, C) ;

  for (r in 1:R) {
    sqrt_d[r] <- sqrt(d[r]) ;
  }

}

parameters {  
  vector[C]  				          delta_noise ; # for coefs on Z
  vector[R]					          beta[J] ;     # for coefs on X & intercept/region-period effect

  real<lower=0>	              tau_mean ;
  vector<lower=0>[J]          tau_unit ;
  real<lower=0,upper=1>       rho_mean ;
  vector<lower=0,upper=1>[J]  rho ;
  real                        alpha0 ; # intercept -> hyperparameter
}

transformed parameters {
  vector[C]  	delta ;
  vector[J]   tau ;
  
  delta <- mu_delta + delta_noise*sigma_delta ;
  tau <- tau_mean*tau_unit ;
}

model {
// declare temp variables
  real shape_1 ;
  real shape_2 ;
  matrix[R,R] K[J] ;                  # for constructing prec matrix 
  vector[R] beta_m_mean[J] ;          # beta - mu_beta for j=1,2,3 & beta - alpha_0 for j = 4.
  row_vector[R] beta_m_mean_t_A[J];   # beta_m_mean' * A


  
// prior on delta_noise
  delta_noise  ~ normal(0, 1) ; # implies delta ~ N(mu_delta, sigma_delta)

// prior on tau_unit (implies tau ~ exponential(1.0 / tau_mean))
  tau_unit ~ exponential(1.0) ;  
  
// prior for rhos 
  shape_1 <- rho_mean * dispersion ;
  shape_2 <- (1.0 - rho_mean) * dispersion ;
  rho ~ beta(shape_1, shape_2) ;

  
// fill in beta_m_mean and beta_m_mean_t_A vectors
  for (j in 1:(J-1) ) {
    beta_m_mean[j] <- beta[j] - mu_beta ;
    beta_m_mean_t_A[j] <- rep_vector(0.0, R)' ;
  }
    beta_m_mean[J] <- beta[J] - alpha0 ;
    beta_m_mean_t_A[J] <- rep_vector(0.0, R)' ;

  for (j in 1:J) {
    for (i in 1:A_N) {
      beta_m_mean_t_A[j, A1[i]] <- (beta_m_mean_t_A[j, A1[i]] 
                                      + beta_m_mean[j, A2[i]]) ;
      beta_m_mean_t_A[j, A2[i]] <- (beta_m_mean_t_A[j, A2[i]] 
                                      + beta_m_mean[j, A1[i]]) ;
    }
  }

// increment log probability manually
  for (j in 1:J) { 
    increment_log_prob(-0.5*tau[j]*dot_self(sqrt_d .* beta_m_mean[j]) ) ;
    increment_log_prob( 0.5*tau[j]*rho[j]*dot_product(beta_m_mean_t_A[j],
                                                      beta_m_mean[j]) ) ;
    K[j] <- -rho[j] * A ;

    for(r in 1:R) {
      K[j,r,r] <- K[j,r,r] + d[r] ;
    }

    increment_log_prob( 0.5 * log_determinant( K[j] ) ) ;
    increment_log_prob( 0.5 * R * log( tau[j] ) ) ;
  }

  
// construct predictor 
 { # begin local block

  vector[N]  	eta_log ; 
          
  eta_log <- Z*delta + M*beta[J] ;

  for (j in 1:(J-1) ) {
    eta_log <- eta_log + X[j]*beta[j] ;
  }
  
  Y ~ bernoulli_logit(eta_log) ;

 } # end local block

}
