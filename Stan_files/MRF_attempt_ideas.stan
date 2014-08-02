data { # stuff to be passed to Stan from R
  
  ## constants ##
  int<lower=1>                      I ; # number of observations
  int<lower=0>                      R ; # number of regions
  
  ## independent variables ##
  vector[I]                         DEM ;
  vector[I]                         LAB ;

  ## x design matrices ##
  matrix[I, R]                      XX ;
  matrix[I, R]                      XX_urb ;
  matrix[I, R]                      XX_aa ;
  matrix[I, R]                      XX_union ;

  ## penalty matrix ##
  matrix[R, R]                      KK ;

  ## outcome variable ##
  int<lower=0,upper=1>              VOTE[I] ; 

}

transformed data {

  vector[R] zeros ;
  matrix[R,R] K_inv ;
  matrix[R,R] K_chol ;

  for (r in 1:R) {
    zeros[r] <- 0 ;
  }

  K_inv <- inverse(KK) ;
  K_chol <- cholesky_decompose(KK) ;
}

parameters {
//  vector[R] beta[4] ;
  vector[R] beta_1 ;
  vector[R] beta_2 ;
  vector[R] beta_3 ;
  vector[R] beta_4 ;
  vector[2] delta ;
  vector<lower=0.00000001>[4] tau_sq ;
}

transformed parameters {
  vector<lower=0>[4] sigma_sq ;

  for (j in 1:4) {
    sigma_sq[j] <- 1/tau_sq[j] ;
  }
  
}

model {
//  vector[I] eta ;

  tau_sq ~ inv_gamma(0.001, 0.001) ;
  
//  for(j in 1:4){
//    beta[j] ~ multi_normal_prec(zeros, KK / tau_sq[j]) ;
//  }

# beta_1 ~ multi_normal_prec(zeros, KK / tau_sq[1]) ;
# beta_2 ~ multi_normal_prec(zeros, KK / tau_sq[2]) ;
# beta_3 ~ multi_normal_prec(zeros, KK / tau_sq[3]) ;
# beta_4 ~ multi_normal_prec(zeros, KK / tau_sq[4]) ;

#beta_1 ~ multi_normal_cholesky(zeros, sigma_sq[1]*K_chol) ;
#beta_2 ~ multi_normal_cholesky(zeros, sigma_sq[2]*K_chol) ;
#beta_3 ~ multi_normal_cholesky(zeros, sigma_sq[3]*K_chol) ;
#beta_4 ~ multi_normal_cholesky(zeros, sigma_sq[4]*K_chol) ;

 beta_1 ~ multi_normal(rep_vector(0, R), KK / tau_sq[1]) ;
 beta_2 ~ multi_normal(zeros, KK / tau_sq[2]) ;
 beta_3 ~ multi_normal(zeros, KK / tau_sq[3]) ;
 beta_4 ~ multi_normal(zeros, KK / tau_sq[4]) ;


// for(j in 1:4){
//  increment_log_prob( - 0.5 * R * log(sigma_sq[j]) + 0.5 * log(determinant(KK)) );
// }

# eta <- (XX*beta_1
#         + XX_urb*beta_2
#         + XX_aa*beta_3
#         + XX_union*beta_4
#         + delta[1]*DEM 
#         + delta[2]*LAB) ;

//  VOTE ~ bernoulli_logit(eta) ;
  
  VOTE ~ bernoulli_logit(XX*beta_1
        + XX_urb*beta_2
        + XX_aa*beta_3
        + XX_union*beta_4
        + delta[1]*DEM 
        + delta[2]*LAB) ;
}