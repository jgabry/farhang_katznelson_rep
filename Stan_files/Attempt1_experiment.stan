data { # stuff to be passed to Stan from R

## constants ##
  int<lower=1>                      I ; # number of observations
  int<lower=1>                      N ; # length of map vector (total number of neighbors)
  int<lower=1>                      R ; # number of regions
  int<lower=0>                      K ; # number of parameters with intrinsic autoregressive priors

## data ##
  real<lower=0, upper=100>          URB[I] ;
  real<lower=0, upper=100>          AA[I] ;
  real<lower=0>                     UNION[I] ;
  int<lower=0, upper=1>             DEM[I] ;
  int<lower=0, upper=1>             LABORCOMM[I] ;
  int<lower=0, upper=1>             VOTE[I] ; 
  
## stuff for spatial smoothing ##
  int<lower=0, upper=R>             region[I] ; 
  int<lower=0, upper=R>             map[N] ;
  int<lower=0, upper=N>             off[R+1] ;
  real<lower=0>                     Nneighs[R] ;  
}

transformed data {
## shape parameter for gamma prior on the taus ##
    real  shape ; 
    shape <- 0.0001 + R / 2.0 ;
}

parameters {
  real          b_DEM ;
  real          b_LABORCOMM ;
  vector[R]     B[K] ;
  real<lower=0> tau[K] ;
}

transformed parameters {
  vector[N]   b_neigh[K] ;
  vector[R]   B_bar[K] ;
  vector[R]   tau_lik[K] ;
  vector[R]   dev[K] ;
  real        d[K] ;

## fill in b_neigh ##
  for (k in 1:K) { 
    for (n in 1:N) {
      b_neigh[k, n] <- B[k, map[n]] ;
    }
  }
 
## fill in B_bar ##
  for (r in 1:R) { 
    real b_neigh_seg[K, off[r + 1] - off[r]] ;
    for (k in 1:K) {
      for (l in (off[r] + 1):off[r+1]) {
        b_neigh_seg[k,l] <- b_neigh[k,l] ; 
      }
      B_bar[k, r] <- mean(b_neigh_seg[k]) ;
    }
  }  



## fill in dev, tau_lik & d ##
  for (k in 1:K) { 
    for (r in 1:R) {
      dev[k, r] <- 1.0 / sqrt(Nneighs[r] * tau[k]) ;
      tau_lik[k, r] <- Nneighs[r] * B[k,r] * (B[k,r] - B_bar[k,r]) ;
    }
    
    d[k] <- 0.0001 + sum(tau_lik[k]) / 2.0;
  }

}

model {

## priors ##
  tau           ~ gamma(shape, d) ;     # gamma(shape, inverse-scale) 
  b_DEM         ~ normal(0.0, 100.0) ;  # normal(mean, sd) 
  b_LABORCOMM   ~ normal(0.0, 100.0) ;

## stochastic functions of hyperparameters ##
  for (k in 1:K) {
    for (r in 1:R) {
      B[k,r] ~ normal(B_bar[k,r], dev[k,r]) ;
    }
  }

  
## likelihood ##
  for(i in 1:I) {
    real eta ; 
    eta <- ( 
            B[1, region[i]] 
          + B[2, region[i]] * URB[i] 
          + B[3, region[i]] * AA[i] 
          + B[4, region[i]] * UNION[i] 
          + b_DEM * DEM[i] 
          + b_LABORCOMM * LABORCOMM[i] 
            ) ;

    VOTE[i] ~ bernoulli_logit(eta) ;
  }

}