### Changes: moved alpha0 from params to transformed data and set it to 0 ###


## BEGIN data BLOCK ##
data {
  // Dimensions 
  int<lower=1>            N ; # number of obs in data
  int<lower=0>            R ; # number of groups (i.e. region/period combos)
  int<lower=0>            J ; # number of betas (including intercept/region-period effect)
  int<lower=0>            C ; # number of deltas (number of vars in Z)
  
  // Outcome variable 
  int<lower=0,upper=1>    Y[N] ; # indicator for prolabor vote
  
  // Predictors w/ effects assumed constant over regional/temporal dims ##
  matrix[N, C]            Z ;  # dem party & labor committee indicators 
  
  // Predictors to get smoothing priors
  matrix[N, R]            X[J] ;  # urbanpct, aapct, unionpop
  
  // Stuff for spatial-smoothing priors
  matrix[R, R]            A ;         # adjacency matrix
  vector[R]               d ;         # diagonal vector of degree matrix
  int                     A_N ;       # number of adjacent region pairs
  int                     A1[A_N] ;   # first half of adjacency pairs (rows)
  int                     A2[A_N] ;   # second half of adjacency pairs (cols)
}
## END data BLOCK ##

## BEGIN transformed data BLOCK ##
transformed data {
  real                    alpha0 ;      # intercept -> hyperparameter
  real<lower=0>           sigma_delta ; # prior sd for deltas
  vector[C]               mu_delta ;    # prior means for deltas
  real                    mu_beta ;     # prior mean for betas (except beta[4], the intercept)
  vector[R]               sqrt_d ;      # sqrt of diagonal elements of degree matrix
  
  alpha0 <- 0.0 ;
  sigma_delta <- 100.0 ;
  mu_beta     <- 0.0 ;
  mu_delta    <- rep_vector(0.0, C) ;
  
  for (r in 1:R) {
    sqrt_d[r] <- sqrt(d[r]) ;
  }
  
}
## END transformed data BLOCK ##

## BEGIN parameters BLOCK ##
parameters {  
  vector[C]  				          delta_noise ; 
  vector[R]					          beta[J] ;     # includes coefs on X & intercept
  
  real<lower=0>	              tau_mean ;
  vector<lower=0>[J]          tau_unit ;
  real<lower=0,upper=1>       rho_mean ;
  real<lower=0>               rho_dispersion ;  
  vector<lower=0,upper=1>[J]  rho ;
}
## END parameters BLOCK ##

## BEGIN transformed parameters BLOCK ##
transformed parameters {
  vector[C]  	delta ;
  vector[J]   tau ;
  
  delta <- mu_delta + delta_noise*sigma_delta ;
  tau <- tau_mean*tau_unit ;
}
## END transformed parameters BLOCK ##

## BEGIN model BLOCK ##
model {
  // declare temp variables
  real shape_1 ;                      # for beta prior on rhos
  real shape_2 ;                      # for beta prior on rhos
  matrix[R,R] K[J] ;                  # for constructing prec matrices
  vector[R] beta_m_mean[J] ;          # beta - mu_beta for j=1,2,3 & beta - alpha_0 for j = 4.
  row_vector[R] beta_m_mean_t_A[J];   # beta_m_mean' * A
  
  
  // prior on alpha0 is implied improper flat prior on (-inf, +inf)
  
  // prior on delta_noise (implies delta ~ normal(mu_delta, sigma_delta))
  delta_noise  ~ normal(0, 1) ; 
  
  // prior on tau_unit (implies tau ~ exponential(1.0 / tau_mean) & taus a priori independent)
  # tau_mean ~ uniform(0, +inf) implied by declared support in parameters block
  tau_unit ~ exponential(1.0) ; 

  
  // prior for rhos
  # rho_mean ~ uniform(0,1) implied by declared support in parameters block
  rho_dispersion ~ exponential(1.0) ;
  shape_1 <- rho_mean * rho_dispersion ;
  shape_2 <- (1.0 - rho_mean) * rho_dispersion ;
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
      increment_log_prob( 0.5*tau[j]*rho[j]*dot_product(beta_m_mean_t_A[j],beta_m_mean[j]) ) ;
      
      K[j] <- -rho[j] * A ;
      
      for(r in 1:R) {
        K[j,r,r] <- K[j,r,r] + d[r] ;
      }
      
      increment_log_prob( 0.5 * log_determinant( K[j] ) ) ;
      increment_log_prob( 0.5 * R * log( tau[j] ) ) ;
    }
    
    
    
{ # begin local block
  
  // construct linear predictor 
  vector[N]  	eta_log ; 
  
  eta_log <- Z*delta;
  
  for (j in 1:J) {
    eta_log <- eta_log + X[j]*beta[j] ;
  }
  
  Y ~ bernoulli_logit(eta_log) ;
  
} # end local block

}
## END model BLOCK ##
