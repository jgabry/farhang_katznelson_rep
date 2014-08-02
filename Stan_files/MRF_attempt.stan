data { # stuff to be passed to Stan from R
  
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
//  matrix[R,R] D;      // diagonal matrix with d_ii = sum(c_ij)
  real<lower=0> sigma_delta ;

  sigma_delta <- 100.0 ;

//  for (i in 1:R) {
//    for (j in 1:R) {
//      D[i,j] <- if_else(i==j, sum(row(KK, i)), 0.0);
//      }
//    }
}


parameters {
  vector[R] beta[4] ;
  vector<lower=0.00000001>[4] tau_sq ;
  vector[2] delta_noise ;
  real<lower=0,upper=1> rho[4]; // strength of spatial correlation
}

transformed parameters {
  vector[2] delta ;
  delta <- sigma_delta*delta_noise ;  
}

model {
  vector[N] eta ;
  matrix[R,R] TAU[4];

  for(j in 1:4){
    TAU[j] <- (1/tau_sq[j]) * (D - rho[j]*KK);
  }
  eta <- (XX*beta[1]
         + XX_urb*beta[2]
         + XX_aa*beta[3]
         + XX_union*beta[4]
         + delta[1]*DEM 
         + delta[2]*LAB) ;

  VOTE ~ bernoulli_logit(eta) ;
  
  for(j in 1:4){
    // beta[j] ~ multi_normal_prec(rep_vector(0, R), KK / tau_sq[j]) ;
      beta[j] ~ multi_normal_prec(rep_vector(0, R), TAU[j]) ;
  }

  delta_noise ~ normal(0, 1) ;

  tau_sq ~ inv_gamma(0.001, 0.001) ;

}