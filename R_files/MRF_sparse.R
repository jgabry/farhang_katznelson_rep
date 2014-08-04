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
                      N = N)


mvn_prec_fit <- stan(file = "my_kyle_edits.stan", data = mvn_prec_data, chains = 0)
fit_test <- stan(fit = mvn_prec_fit, data = mvn_prec_data, chains = 1, iter = 30)


library(parallel)
options(mc.cores = 6L)
mvn_prec_fit_1k_new <- mclapply(X = 1:6, 
                                mc.cores = 6, 
                                FUN = function(k) {
                                  stan(fit = mvn_prec_fit, 
                                       seed = 123456, 
                                       data = mvn_prec_data,
                                       chains = 1,
                                       iter = 1000,
                                       chain_id = k,
                                       refresh = -1) })




