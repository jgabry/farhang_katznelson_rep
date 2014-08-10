# Written by Ben Goodrich

library(rstan)
library(energy)

convergence <- function(object, R = 0, thin = 1, ...) {
  mat <- as.matrix(object)
  mat <- mat[1:nrow(mat) %% thin == 0,]
  chains <- if(is.list(object)) length(object) else ncol(object)
  n <- nrow(mat) / chains
  chain_id <- rep(1:chains, each = n)
  test <- disco(mat, factors = chain_id, distance = FALSE, R = R, ...)
  W <- test$within / test$Df.e
  B <- test$between /  (chains - 1)
  names(B) <- NULL
  total <- (n - 1) / n *  W + B / n
  loss <- sqrt(total / W)
  test$loss <- loss
  print(paste("Multivariate loss =", round(loss, 3)))
  return(test)
}


# good
# dogs <- stan_demo("dogs", chains = 8, thin = 10)
# convergence(dogs, R = 100)

# bad
# kidney <- stan_demo("kidney", chains = 8, thin = 10)
# convergence(kidney, R = 100)

# hepatitis <- stan_demo("hepatitis", chains = 8)
# convergence(hepatitis, R = 0)
