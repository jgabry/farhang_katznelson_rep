data { # stuff to be passed to Stan from R
  
  ## constants ##
  int<lower=1>                      N ; # number of observations
  int<lower=0>                      R ; # number of regions
  int<lower=1>                      Krank ; 
  
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

  ## outcome variable ##
  int<lower=0,upper=1>              VOTE[N] ; 

}

transformed data {
  real<lower=0> sigma_delta ;
  sigma_delta <- 100.0 ;
}

parameters {
  vector[Krank] beta_free[4] ;
  vector<lower=0.00000001>[4] tau_sq ;
  vector[2] delta_noise ;
}

transformed parameters {  
  vector[R - Krank] beta_pinned[4] ;
  vector[2] delta ;

  delta <- sigma_delta*delta_noise ;  
  
  for(j in 1:4){ 
    for(r in 1:(R - Krank)){
      beta_pinned[j,r] <- (transpose(beta_free[j]) * block(KK, 1, 1, Krank, Krank) 
                          * beta_free[j]) ;
    }
  }
   

}

model {
  vector[N] eta ;
  vector[R] beta[4];

  for(j in 1:4){
    for(r in 1:Krank) {
      beta[j, r] <- beta_free[j,r] ;
    }
    for(i in 1:(R - Krank)) {
      beta[j, i] <- beta_pinned[j, i] ;
    }
  }
  
  eta <- (XX*beta[1]
         + XX_urb*beta[2]
         + XX_aa*beta[3]
         + XX_union*beta[4]
         + delta[1]*DEM 
         + delta[2]*LAB) ;

  VOTE ~ bernoulli_logit(eta) ;
  
//  for(j in 1:4){
  //  beta_free[j] ~ multi_normal(rep_vector(0, Krank), M*tau_sq[j]) ;
//  }

  for(j in 1:4){
    beta[j] ~ multi_normal_prec(rep_vector(0, R), KK/tau_sq[j]) ;
  }

  delta_noise ~ normal(0, 1) ;

  tau_sq ~ inv_gamma(0.001, 0.001) ;

}