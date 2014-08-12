### Changes ###

# 1. removed M incidence matrix. M is now the last (i.e. Jth) matrix in the X array

# 2. dispersion now named rho_dispersion and moved from transformed data to parameters block  
# and given exponential(1.0) prior instead of fixed at 1


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
//  matrix[N, R]            X[J] ;  # urbanpct, aapct, unionpop
  
  // Stuff for spatial-smoothing priors
  matrix[R, R]   	        A ;         # adjacency matrix
  matrix[N, R]            M ;
  vector[R]               d ;         # diagonal vector of degree matrix
  int                     A_N ;       # number of adjacent region pairs
  int                     A1[A_N] ;   # first half of adjacency pairs (rows)
  int                     A2[A_N] ;   # second half of adjacency pairs (cols)
}
## END data BLOCK ##

## BEGIN transformed data BLOCK ##
transformed data {
  real<lower=0>           sigma_delta ; # prior sd for deltas
  vector[C]               mu_delta ;    # prior means for deltas
  vector[R]               sqrt_d ;      # sqrt of diagonal elements of degree matrix
  
  sigma_delta <- 100.0 ;
  mu_delta    <- rep_vector(0.0, C) ;
  
  for (r in 1:R) {
    sqrt_d[r] <- sqrt(d[r]) ;
  }
  
}
## END transformed data BLOCK ##

## BEGIN parameters BLOCK ##
parameters {  
  vector[C]  				          delta_noise ; 
  vector[R]					          alpha ;
  
  real<lower=0>               tau ;
  real<lower=0,upper=1>       rho ;
}
## END parameters BLOCK ##

## BEGIN transformed parameters BLOCK ##
transformed parameters {
  vector[C]  	delta ;
  
  delta <- mu_delta + delta_noise*sigma_delta ;
}
## END transformed parameters BLOCK ##

## BEGIN model BLOCK ##
model {
  // declare temp variables
  matrix[R,R] K ;                  # for constructing prec matrices
  vector[R] alpha_m_mean ;          # beta - mu_beta for j=1,2,3 & beta - alpha_0 for j = 4.
  row_vector[R] alpha_m_mean_t_A;   # alpha_m_mean' * A
  
  
  // prior on alpha0 is implied improper flat prior on (-inf, +inf)
  
  // prior on delta_noise (implies delta ~ normal(mu_delta, sigma_delta))
  delta_noise  ~ normal(0, 1) ; 
  
  tau ~ gamma(0.5, 0.005) ;
  
  
  // fill in alpha_m_mean and alpha_m_mean_t_A vectors

    alpha_m_mean <- alpha - 0.0 ;
    alpha_m_mean_t_A <- rep_vector(0.0, R)' ;
    
      for (i in 1:A_N) {
        alpha_m_mean_t_A[A1[i]] <- (alpha_m_mean_t_A[A1[i]] 
                                      + alpha_m_mean[A2[i]]) ;
        alpha_m_mean_t_A[A2[i]] <- (alpha_m_mean_t_A[A2[i]] 
                                      + alpha_m_mean[A1[i]]) ;
      }
    
    // increment log probability manually 
      increment_log_prob(-0.5*tau*dot_self(sqrt_d .* alpha_m_mean) ) ;
      increment_log_prob( 0.5*tau*rho*dot_product(alpha_m_mean_t_A,alpha_m_mean) ) ;
      
      K <- -rho * A ;
      
      for(r in 1:R) {
        K[r,r] <- K[r,r] + d[r] ;
      }
      
      increment_log_prob( 0.5 * log_determinant( K ) ) ;
      increment_log_prob( 0.5 * R * log( tau ) ) ;
    }
    
    
    
{ # begin local block
  
  // construct linear predictor 
  vector[N]  	eta_log ; 
  
  eta_log <- Z*delta + M*alpha ;
  
  Y ~ bernoulli_logit(eta_log) ;
  
} # end local block

}
## END model BLOCK ##
