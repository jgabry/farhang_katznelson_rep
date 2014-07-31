data { # stuff to be passed to Stan from R
  int<lower=0>                      I ; # number of observations
  int<lower=0>                      N ; # length of map vector
  int<lower=0>                      R ; # number of regions
  int<lower=0>                      region[I] ; 
  real<lower=0, upper=100>          URB[I] ;
  real<lower=0, upper=100>          AA[I] ;
  real<lower=0>                     UNION[I] ;
  int<lower=0, upper=1>             DEM[I] ;
  int<lower=0, upper=1>             LABORCOMM[I] ;
  int<lower=0, upper=1>             VOTE[I] ; 
  int<lower=0, upper=R>             map[N] ;
  int<lower=0, upper=N>             off[R+1] ;
  real<lower=0>                     Nneighs[R] ;  
}

transformed data {
//  int   Nneighs[R] ;

//  for(r in 1:R){
//    Nneighs[r] <- off[r+1] - off[r] ;
//  }

}

parameters {
  real  b_DEM ;
  real  b_LABORCOMM ;
  real  b[R] ;
  real  b_URB[R] ;
  real  b_AA[R] ;
  real  b_UNION[R] ;
}

transformed parameters {

  real b_neigh[N] ;

  real b_bar[R] ;
  real b_URB_bar[R] ;
  real b_AA_bar[R] ;
  real b_UNION_bar[R] ;

  for(n in 1:N){
    b_neigh[n] <- b[ map[n] ] ;
  }

  for(r in 1:R){
    real seg[ off[r+1] - off[r] ] ;
    seg <- segment(b_neigh, off[r] + 1, off[r+1] - off[r]) ; 
    b_bar[r] <- mean(seg) ;
    b_URB_bar[r] <- mean(seg) ;
    b_AA_bar[r] <- mean(seg) ;
    b_UNION_bar[r] <- mean(seg) ;
  }
}

model {

### priors ###
  b_DEM       ~ normal(0, 100.0) ; 
  b_LABORCOMM ~ normal(0, 100.0) ; 

  for(r in 1:R){
    real sigma ; 
    sigma <- 1/sqrt(Nneighs[r]) ;

    b[r]        ~ normal(b_bar[r], sigma) ;
    b_URB[r]    ~ normal(b_URB_bar[r], sigma) ;
    b_AA[r]     ~ normal(b_AA_bar[r], sigma) ;
    b_UNION[r]  ~ normal(b_UNION_bar[r], sigma) ;
  }
  
### likelihood ###
  for(i in 1:I) {
    real eta ;
    eta <- ( 
            b[ region[i] ] 
          + b_URB[ region[i] ] * URB[i] 
          + b_AA[ region[i] ] * AA[i] 
          + b_UNION[ region[i] ] * UNION[i] 
          + b_DEM * DEM[i] 
          + b_LABORCOMM * LABORCOMM[i]
            ) ;

    VOTE[i] ~ bernoulli_logit(eta) ;
  }

}