### changes ###
# 1. rhos get beta(1.5, 1.5) priors
# 2. unionpop converted to proportion

## BEGIN data ##
data {
  // Dimensions 
  int<lower=1>            N ; # number of obs in data
  int<lower=0>            R ; # number of groups (i.e. region/period combos)
  int<lower=0>            C ; # number of deltas (number of vars in Z)
  int<lower=0>            J ; # number of beta vectors + 1 for alpha vector
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
}
## END data ##

## BEGIN transformed data ##
transformed data {
  real<lower=0>   sigma_alpha0 ;  # prior sd for alpha0
  real<lower=0>   sigma_delta ;   # prior sd for deltas
  vector[C]       mu_delta ;      # prior means for deltas
  real            mu_alpha0 ;     # prior mean for alpha0
  
  sigma_alpha0  <- 3.0 ;
  sigma_delta   <- 100.0 ;
  mu_delta      <- rep_vector(0.0, C) ;
  mu_alpha0     <- 0.0 ;
}
## END transformed data ##

## BEGIN parameters ##
parameters {  
  real                        alpha0_noise ;
  row_vector[R]               alpha_noise ;
  row_vector[R]               beta_noise[J-1] ;
  vector<lower=0>[J]          tau_noise ; 
  vector<lower=0,upper=1>[J]  rho ;
  vector[P]                   beta_aa_noise ;
  cholesky_factor_corr[P]     L_aa ;
  vector<lower=0>[P]          scale_aa_noise ;
  vector[C]        	          delta_noise ; 
}
## END parameters ##

## BEGIN transformed parameters ##
transformed parameters {
  real                alpha0 ;
  vector[C]           delta ;
  vector[P]           scale_aa ;
  vector[P]           beta_aa ;
  vector<lower=0>[J]  tau ;
  vector[R]           beta[J-1] ;
  vector[R]           alpha ;


  alpha0    <- mu_alpha0 + alpha0_noise*sigma_alpha0 ;
  delta     <- mu_delta + delta_noise*sigma_delta ;
  
  for (p in 1:P) {
    scale_aa[p] <- 3 * tan(pi() * (Phi_approx(scale_aa_noise[p]) - 0.5)) ;
  }
  beta_aa <- scale_aa .* (L_aa * beta_aa_noise) ; # note: .* is element-wise mult.
   
  for (j in 1:J) {
    tau[j] <- tan( pi() * (Phi_approx(tau_noise[j]) - 0.5) ) ;
  }
  
  for (j in 1:(J-1)) {
    matrix[R,R] precision ;
    precision <- -rho[j] * A ;
    for (r in 1:R) {
      precision[r,r] <- precision[r,r] + d[r] ;
    }
    beta[j] <- mdivide_right_tri_low(beta_noise[j],
               cholesky_decompose(precision))' / sqrt(tau[j]) ;
  }

  { # BEGIN local
    matrix[R,R] precision;
    precision <- -rho[J] * A;
    for (r in 1:R) {
      precision[r,r] <- precision[r,r] + d[r] ;
    }
    alpha <- alpha0 + mdivide_right_tri_low(alpha_noise,
             cholesky_decompose(precision))' / sqrt(tau[J]);
  } # END local

}
## END transformed parameters ##

## BEGIN model ##
model { 
  // prior for rhos
  rho ~ beta(1.5, 1.5) ;

  // std. normal dist. for the "noises" 
  delta_noise   ~ normal(0, 1) ;    // implies delta ~ normal(mu_delta, sigma_delta)
  tau_noise     ~ normal(0, 1) ;    // implies tau ~ half-cauchy(0,1) 
  alpha0_noise  ~ normal(0, 1) ;    // implies alpha0 ~ normal(mu_alpha0, sigma_alpha0)
  alpha_noise   ~ normal(0, 1) ;    // implies alpha ~ MVN(alpha0, precision) 
  for (j in 1:(J-1)) {
    beta_noise[j] ~ normal(0, 1) ;  // implies beta[j] ~ MVN(0, precision) 
  }  

  
  // stuff for beta_aa
  L_aa ~ lkj_corr_cholesky(1) ;
  scale_aa_noise ~ normal(0,1) ;    // implies scale_aa ~ cauchy(0,3) ;
  beta_aa_noise ~ normal(0,1) ;
    # implies beta_aa ~ multi_normal(rep_vector(0, P), 
    #           diag_matrix(scale_aa) * L * L' diag_matrix(scale_aa))

        
  { # BEGIN local (for likelihood)
  
  // construct linear predictor 
    vector[N]    eta_log ; 
  
    eta_log <- Z*delta + X_aa*beta_aa + X[J]*alpha;
  
    for (j in 1:(J-1)) {
      eta_log <- eta_log + X[j]*beta[j] ;
    }
  
    Y ~ bernoulli_logit(eta_log) ;
  
  } # END local 

}
## END model ##