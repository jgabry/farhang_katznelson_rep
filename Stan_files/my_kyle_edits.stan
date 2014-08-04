data {
// Dimensions of the data
  int<lower=1>            N ; # number of obs in data
  int<lower=0>  	        R ; # number of groups (i.e. region/period combos)


// Outcome variable 
  int<lower=0,upper=1>    Y[N] ; # indicator for prolabor vote
  
// Predictors w/ effects assumed constant over regional/temporal dims ##
  matrix[N, 2]		        Z ;  # dem party & labor committee indicators 
  
// Predictors to get smoothing priors
  matrix[N, R]		        X[3] ;  # urbanpct, aapct, unionpop
  
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
  vector[2]               mu_delta ;    # mean for priors on delta
  vector[R]               sqrt_d ;      # diagonal elements of degree matrix
  

  sigma_delta <- 100.0 ;
  mu_delta    <- rep_vector(0.0, 2) ;

  for (r in 1:R) {
    sqrt_d[r] <- sqrt(d[r]) ;
  }

}

parameters {  
  vector[2]  				          delta_noise ; # for coefs on Z
  vector[R]					          beta[4] ;     # for coefs on X & reg-per effect

  vector<lower=0.000001>[4]	  tau ;
  real<lower=0,upper=1>       rho[4] ;
  real                        alpha0 ; # intercept -> hyperparameter
}

transformed parameters {
  vector[2]  	delta ;
  
  delta <- mu_delta + delta_noise*sigma_delta ;
}

model {
// declare temp variables
  matrix[R,R] K[4] ;                  # tau*K = prec matrix 
  vector[R] beta_m_alpha0[4] ;        # beta - alpha0
  row_vector[R] beta_m_alpha0t_A[4];  # (beta - alpha0)' * A
  
// prior on delta_noise
  delta_noise  ~ normal(0, 1) ; # implies delta ~ N(mu_delta, sigma_delta)

// prior on taus
  tau ~ gamma(0.5, 0.0005);  
  
// no explicit prior on rhos: implied uniform prior over support 


  
// find (beta - alpha0)' * A
  for (j in 1:4) {
    beta_m_alpha0[j] <- beta[j] - alpha0 ;
    beta_m_alpha0t_A[j] <- rep_vector(0.0, R)' ;
  }

  for (j in 1:4) {
    for (i in 1:A_N) {
      beta_m_alpha0t_A[j, A1[i]] <- (beta_m_alpha0t_A[j, A1[i]] 
                                      + beta_m_alpha0[j, A2[i]]) ;
      beta_m_alpha0t_A[j, A2[i]] <- (beta_m_alpha0t_A[j, A2[i]] 
                                      + beta_m_alpha0[j, A1[i]]) ;
    }
  }

// increment log probability manually
  for (j in 1:4) { 
    increment_log_prob(-0.5*tau[j]*dot_self(sqrt_d .* beta_m_alpha0[j]) ) ;
    increment_log_prob( 0.5*tau[j]*rho[j]*dot_product(beta_m_alpha0t_A[j],
                                                      beta_m_alpha0[j]) ) ;
    K[j] <- -rho[j] * A ;

    for(r in 1:R) {
      K[j,r,r] <- K[j,r,r] + d[r] ;
    }

    increment_log_prob( 0.5 * log_determinant( K[j] ) ) ;
    increment_log_prob( 0.5 * R * log( tau[j] ) ) ;
  }

  
// construct linear predictor eta 
 { # local block
  vector[N]  	eta ; 

  eta <- (
          Z*delta        # dem & laborcomm
          + X[1]*beta[1] # urbanpct
          + X[2]*beta[2] # aapct
          + X[3]*beta[3] # unionpop
          + M*beta[4]    # region-period effect
          ) ;
  
  Y ~ bernoulli_logit(eta) ;
 }

}
