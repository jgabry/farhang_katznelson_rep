data { # stuff to be passed to Stan from R

## constants ##
  int<lower=1>                      I ; # number of observations
  int<lower=1>                      N ; # length of map vector
  int<lower=0>                      R ; # number of regions

## independent variables ##
  real<lower=0,upper=100>           URB[I] ;
  real<lower=0,upper=100>           AA[I] ;
  real<lower=0>                     UNION[I] ;
  int<lower=0,upper=1>              DEM[I] ;
  int<lower=0,upper=1>              LABORCOMM[I] ;

## outcome variable ##
  int<lower=0,upper=1>              VOTE[I] ; 

## stuff for spatial smoothing ##
  int<lower=0,upper=R>              region[I] ; 
  int<lower=0,upper=R>              map[N] ;
  int<lower=0,upper=N>              off[R + 1] ;
  real<lower=0>                     Nneighs[R] ;  
}

transformed data {
## shape parameter for gamma priors for the taus ##
  real              shape ; 
  shape <- 0.0001 + R / 2.0 ;
}

parameters {
## coefficients on predictors with effects assumed constant over regional/temporal dims ##
  real  b_DEM ;
  real  b_LABORCOMM ;

## parameters to be given the instrinsic autoregressive priors ## 
  real  b_0[R] ;
  real  b_URB[R] ;
  real  b_AA[R] ;
  real  b_UNION[R] ;

## hyperparameters ##
  real<lower=0> tau_b_0 ;
  real<lower=0> tau_b_URB ;
  real<lower=0> tau_b_AA ;
  real<lower=0> tau_b_UNION ;
}

transformed parameters {

  real  b_0_neigh[N] ;
  real  b_URB_neigh[N] ;
  real  b_AA_neigh[N] ;
  real  b_UNION_neigh[N] ;
  
  real  b_0_bar[R] ;
  real  b_URB_bar[R] ;
  real  b_AA_bar[R] ;
  real  b_UNION_bar[R] ;
  
  real  tau_b_0_lik[R] ;
  real  tau_b_URB_lik[R] ;
  real  tau_b_AA_lik[R] ;
  real  tau_b_UNION_lik[R] ;
  real  b_0_dev[R] ;
  real  b_URB_dev[R] ;
  real  b_AA_dev[R] ;
  real  b_UNION_dev[R] ;
  
  real<lower=0>  d_b_0 ;
  real<lower=0>  d_b_URB ;
  real<lower=0>  d_b_AA ;
  real<lower=0>  d_b_UNION ;
  

  for(n in 1:N){
    b_0_neigh[n] <- b_0[ map[n] ] ;
    b_URB_neigh[n] <- b_URB[ map[n] ] ;
    b_AA_neigh[n] <- b_AA[ map[n] ] ;
    b_UNION_neigh[n] <- b_UNION[ map[n] ] ;
  }
  
  for(r in 1:R){
    real seg_0[ off[r+1] - off[r] ] ;
    real seg_URB[ off[r+1] - off[r] ] ;
    real seg_AA[ off[r+1] - off[r] ] ;
    real seg_UNION[ off[r+1] - off[r] ] ;
    seg_0               <- segment(b_0_neigh, off[r] + 1, off[r+1] - off[r]) ; 
    seg_URB             <- segment(b_URB_neigh, off[r] + 1, off[r+1] - off[r]) ; 
    seg_AA              <- segment(b_AA_neigh, off[r] + 1, off[r+1] - off[r]) ; 
    seg_UNION           <- segment(b_UNION_neigh, off[r] + 1, off[r+1] - off[r]) ; 
    b_0_bar[r]          <- mean(seg_0) ;
    b_URB_bar[r]        <- mean(seg_URB) ;
    b_AA_bar[r]         <- mean(seg_AA) ;
    b_UNION_bar[r]      <- mean(seg_UNION) ;
    b_0_dev[r]          <- 1.0 / sqrt(Nneighs[r] * tau_b_0) ;
    b_URB_dev[r]        <- 1.0 / sqrt(Nneighs[r] * tau_b_URB) ;
    b_AA_dev[r]         <- 1.0 / sqrt(Nneighs[r] * tau_b_AA) ;
    b_UNION_dev[r]      <- 1.0 / sqrt(Nneighs[r] * tau_b_UNION) ;
    tau_b_0_lik[r]      <- Nneighs[r] * b_0[r] * (b_0[r] - b_0_bar[r]) ;
    tau_b_URB_lik[r]    <- Nneighs[r] * b_URB[r] * (b_URB[r] - b_URB_bar[r]) ;
    tau_b_AA_lik[r]     <- Nneighs[r] * b_AA[r] * (b_AA[r] - b_AA_bar[r]) ;
    tau_b_UNION_lik[r]  <- Nneighs[r] * b_UNION[r] * (b_UNION[r] - b_UNION_bar[r]) ;
  }
  
  d_b_0     <- 0.0001 + sum(tau_b_0_lik) / 2.0;
  d_b_URB   <- 0.0001 + sum(tau_b_URB_lik) / 2.0;
  d_b_AA    <- 0.0001 + sum(tau_b_AA_lik) / 2.0;
  d_b_UNION <- 0.0001 + sum(tau_b_UNION_lik) / 2.0;
  
}

model {
  
## priors ##
  tau_b_0       ~ gamma(shape, d_b_0) ;   
  tau_b_URB     ~ gamma(shape, d_b_URB) ; 
  tau_b_AA      ~ gamma(shape, d_b_AA) ; 
  tau_b_UNION   ~ gamma(shape, d_b_UNION) ; 
  
  b_DEM         ~ normal(0.0, 100.0) ; 
  b_LABORCOMM   ~ normal(0.0, 100.0) ;
  
## stochastic functions of hyperparameters ##
  for(r in 1:R){
    b_0[r]      ~ normal(b_0_bar[r], b_0_dev[r]) ;
    b_URB[r]    ~ normal(b_URB_bar[r], b_URB_dev[r]) ;
    b_AA[r]     ~ normal(b_AA_bar[r], b_AA_dev[r]) ;
    b_UNION[r]  ~ normal(b_UNION_bar[r], b_UNION_dev[r]) ;
  }
  
  
## likelihood ##
  for(i in 1:I) {
    real eta ; 
    eta <- Phi_approx( 
      b_0[ region[i] ] 
      + b_URB[ region[i] ] * URB[i] 
      + b_AA[ region[i] ] * AA[i] 
      + b_UNION[ region[i] ] * UNION[i] 
      + b_DEM * DEM[i] 
      + b_LABORCOMM * LABORCOMM[i] 
    ) ;
    VOTE[i] ~ bernoulli(eta) ;
  }
  
}