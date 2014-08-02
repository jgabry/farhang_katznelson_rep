data { 
  ## constants ##
  int<lower=1>                      N ; # number of observations
  int<lower=0>                      R ; # number of regions
  
  ## predictors with effects assumed constant over regional/temporal dims ##
  vector[N]                         DEM ;
  vector[N]                         LAB ;
  
  ## design matrices for the other predictors ##
  matrix[N, R]                      XX ;
  matrix[N, R]                      XX_urb ;
  matrix[N, R]                      XX_aa ;
  matrix[N, R]                      XX_union ;
  
  ## penalty matrix ##
  matrix[R, R]                      KK ;
  matrix[R, R]                      D ;
  
  ## outcome variable ##
  int<lower=0,upper=1>              VOTE[N] ; 
}

transformed data {
  real<lower=0> sigma_delta ;
  sigma_delta <- 100.0 ;
}


parameters {
  vector[R]                         beta[4] ;
  vector<lower=0.00000001>[4]       tau_sq ;
  vector[2]                         delta_noise ;
  real<lower=0, upper=1>            rho[4] ; // strength of spatial correlation
}

transformed parameters {
  vector[2] delta ;
  vector[4] sigma_sq;
  
  delta <- sigma_delta*delta_noise ;  
  for(j in 1:4) sigma_sq[j] <- 1/tau_sq[j] ;
}

model {
  vector[N] eta ;
  
  eta <- (XX*beta[1]
          + XX_urb*beta[2]
          + XX_aa*beta[3]
          + XX_union*beta[4]
          + delta[1]*DEM 
          + delta[2]*LAB) ;
  
  VOTE ~ bernoulli_logit(eta) ;
  
  
  for(j in 1:4){
    increment_log_prob(0.5 * R * log( sigma_sq[j] ) 
                       + 0.5 * log( determinant(D - rho[j] * KK) ) );
  }
  
  delta_noise ~ normal(0, 1) ;
  
  tau_sq ~ inv_gamma(0.001, 0.001) ;
  
}
