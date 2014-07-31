data { # stuff to be passed to Stan from R
  int<lower=0>            I ;         # number of observations
  int<lower=0>            N ;         # length of map vector
  int<lower=0>            R ;         # number of unique regions
  int<lower=0, upper=R>   region[I] ; 
  int<lower=0, upper=1>   DEM[I] ;
  int<lower=0, upper=1>   VOTE[I] ; 
  int<lower=0, upper=R>   map[N] ;
  int<lower=0, upper=N>   off[R+1] ;
  real<lower=0>           Nneighs[R] ;  
}

parameters {
  real  b_DEM ;
  real  b[R] ;
}

transformed parameters {
  real b_neigh[N] ;
  real b_bar[R] ;

  for(n in 1:N){
    b_neigh[n] <- b[ map[n] ] ;
  }

  for(r in 1:R){
    real seg[off[r+1] - off[r]] ;
    seg <- segment(b_neigh, off[r] + 1, off[r+1] - off[r]) ; 
    b_bar[r] <- mean(seg) ;
  }
}

model {

  for(r in 1:R) {
    b[r] ~ normal(b_bar[r], 1.0/sqrt(Nneighs[r])) ;
  }
  
  b_DEM ~ normal(0, 100.0) ; 

  for(i in 1:I) {
    real eta ;
    eta <- ( b[region[i]]  + b_DEM * DEM[i] ) ;
    VOTE[i] ~ bernoulli_logit(eta) ;
  }

}