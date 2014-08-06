setwd("/Users/jgabry/Desktop/COLUMBIA/Stuff_for_Wawro/Rsync/The_Real_Deal")

library(rstan)
set_cppo('fast')

source('MRF_source.R')

D <- diag(rowSums(A))


### sparse matrix stuff ###
A_sparse <- which(A == 1, arr.ind=TRUE)
  # remove duplicates (because matrix is symmetric)
A_sparse <- A_sparse[A_sparse[,1] < A_sparse[,2],]
A_N <- dim(A_sparse)[1]
A1 <- A_sparse[,1]
A2 <- A_sparse[,2]

mvn_prec_data <- list(Z = Z,
                      X = X,
                      Y = Y,
                      M = M,
                      A = A, 
                      A_N = A_N, 
                      A1 = A1, 
                      A2 = A2,
                      d = diag(D),
                      R = R,
                      N = N,
                      C = 2, 
                      J = 4)

    
mvn_prec_inits <-function(){
      list(
        alpha0 = 0.0,
        beta = matrix(0.0, nrow = R, nc = 4),
        tau = rep(1.0, 4),
        rho = rep(0.99, 4),
        delta_noise = rnorm(2)
        )
    }
  
mvn_prec_fit <- stan(file = "SpatialStan_08_04_2014.stan", data = mvn_prec_data, chains = 0)
fit_test <- stan(fit = mvn_prec_fit, data = mvn_prec_data, chains = 1, iter = 30)


library(parallel)
options(mc.cores = 6L)
mvn_prec_fit_spatial <- mclapply(X = 1:2, 
                                mc.cores = 6, 
                                FUN = function(k) {
                                  stan(fit = mvn_prec_fit, 
                                       seed = 123478, 
                                       data = mvn_prec_data,
                                       chains = 1,
                                       iter = 2000,
                                       chain_id = k,
                                       refresh = -1) })

mvn_prec_fit_spatial_big <- mclapply(X = 1:6, mc.cores = 6, FUN = function(k) {
  stan(fit = mvn_prec_fit, seed = 13579, data = mvn_prec_data, 
       chains = 1, iter = 4000, chain_id = k, refresh = -1) 
  })
save(mvn_prec_fit_spatial_big, file = "mvn_prec_fit_spatial_big.RData")


