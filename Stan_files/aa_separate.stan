### Changes ###

# 1. rhos and taus iid unif(0,1) and half-cauchy(0,1), respectively

## BEGIN data BLOCK ##
data {
  // Dimensions 
  int<lower=1>            N ; # number of obs in data
  int<lower=0>            R ; # number of groups (i.e. region/period combos)
  int<lower=0>            C ; # number of deltas (number of vars in Z)
  int<lower=0>            J ; # number of beta vectors (including intercepts)
  int<lower=0>            P ; # number of spatial regions (for aa) 
  
  // Outcome variable 
  int<lower=0,upper=1>    Y[N] ; # indicator for prolabor vote
  
  // Predictors w/ effects assumed constant over regional/temporal dims ##
  matrix[N, C]            Z ;  # dem party & labor committee indicators 

  // Predictors to get smoothing priors
  matrix[N, R]            X[J] ;  # urban proportion, union pop & intercepts

  matrix[N, P]            X_aa ;  # aa proportion 
  
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
  real<lower=0>           sigma_alpha0 ;  # prior sd for alpha0
  real<lower=0>           sigma_delta ;   # prior sd for deltas
  vector[C]               mu_delta ;      # prior means for deltas
  real                    mu_alpha0 ;     # prior mean for alpha0
  real                    mu_beta ;       # prior mean for betas 1 through J-1
  vector[R]               sqrt_d ;        # sqrt of diagonal elements of degree matrix
  
  sigma_delta   <- 100.0 ;
  mu_delta      <- rep_vector(0.0, C) ;
  mu_alpha0     <- 0.0 ;
  sigma_alpha0  <- 3.0 ;
  mu_beta       <- 0.0 ;
  
  for (r in 1:R) {
    sqrt_d[r] <- sqrt(d[r]) ;
  }
  
}
## END transformed data BLOCK ##

## BEGIN parameters BLOCK ##
parameters {  
  cholesky_factor_corr[P]     L_aa ;
  vector<lower=0>[P]          scale_aa_noise ;
  vector[P]                   beta_aa_noise ;
  vector[R]                   beta[J] ;
  vector[C]    			          delta_noise ; 
  real                        alpha0_noise ;
  vector<lower=0>[J]          tau_noise ; 
  vector<lower=0,upper=1>[J]  rho ;
}
## END parameters BLOCK ##

## BEGIN transformed parameters BLOCK ##
transformed parameters {
  vector[P]           beta_aa ;
  vector[P]           scale_aa ;
  vector[C]  	        delta ;
  vector<lower=0>[J]  tau ;
  real                alpha0 ;
  
  # note period before star to indicate element-wise multiplication
  for(p in 1:P) {
    scale_aa[p] <- 3 * tan(pi() * (Phi_approx(scale_aa_noise[p]) - 0.5)) ;
  }
  beta_aa <- scale_aa .* (L_aa * beta_aa_noise) ;

  alpha0    <- mu_alpha0 + alpha0_noise*sigma_alpha0 ;
  delta     <- mu_delta + delta_noise*sigma_delta ;

  for (j in 1:J) {
      tau[j] <- tan( pi() * (Phi_approx(tau_noise[j]) - 0.5) ) ;
  }

}
## END transformed parameters BLOCK ##

## BEGIN model BLOCK ##
model {
  // declare temp variables
  matrix[R,R] K[J] ;                  # for constructing precision matrices
  vector[R] beta_m_mean[J] ;          # beta - mu_beta (or beta - alpha0 for beta[4,]) 
  row_vector[R] beta_m_mean_t_A[J] ;  # beta_m_mean' * A 
  
  // standard normal distributions for the "noises" 
  delta_noise   ~ normal(0, 1) ;    // implies delta ~ normal(mu_delta, sigma_delta)
  tau_noise     ~ normal(0, 1) ;    // implies tau ~ half-cauchy(0,1) 
  alpha0_noise  ~ normal(0, 1) ;    // implies alpha0 ~ normal(mu_alpha0, sigma_alpha0)

  // prior for rhos is implied uniform(0,1) 


  // priors for beta_aa
  L_aa ~ lkj_corr_cholesky(1) ;
  beta_aa_noise ~ normal(0,1) ;
      # implies beta_aa ~ multi_normal(rep_vector(0, P), 
      #           diag_matrix(scale_aa) * L * L' diag_matrix(scale_aa))
  scale_aa_noise ~ normal(0,1) ;    // implies scale_aa ~ cauchy(0,3) ;

  
  // fill in beta_m_mean and beta_m_mean_t_A vectors
  for (j in 1:(J-1) ) {
    beta_m_mean[j] <- beta[j] - mu_beta ;
    beta_m_mean_t_A[j] <- rep_vector(0.0, R)' ;
  }
    beta_m_mean[J] <- beta[J] - alpha0 ;
    beta_m_mean_t_A[J] <- rep_vector(0.0, R)' ;
    
  for (j in 1:J) {
    for (i in 1:A_N) {
      beta_m_mean_t_A[j, A1[i]] <- (beta_m_mean_t_A[j, A1[i]] + beta_m_mean[j, A2[i]]) ;
      beta_m_mean_t_A[j, A2[i]] <- (beta_m_mean_t_A[j, A2[i]] + beta_m_mean[j, A1[i]]) ;
    }
  }
    
// increment log probability manually for multi_normal_prec  
  for (j in 1:J) { 
    increment_log_prob(-0.5*tau[j]*dot_self(sqrt_d .* beta_m_mean[j]) ) ;
    increment_log_prob( 0.5*tau[j]*rho[j]*dot_product(beta_m_mean_t_A[j], beta_m_mean[j]) ) ;
      
    K[j] <- -rho[j] * A ;
      
    for(r in 1:R) {
      K[j,r,r] <- K[j,r,r] + d[r] ;
    }
      
    increment_log_prob( 0.5 * log_determinant( K[j] ) ) ;
    increment_log_prob( 0.5 * R * log( tau[j] ) ) ;
  }
    

{ # begin local block for likelihood
  
// construct linear predictor 
  vector[N]    eta_log ; 
  
  eta_log <- Z*delta + X_aa*beta_aa;
  
  for (j in 1:J) {
    eta_log <- eta_log + X[j]*beta[j] ;
  }
  
  Y ~ bernoulli_logit(eta_log) ;
  
} # end local block

}
## END model BLOCK ##
