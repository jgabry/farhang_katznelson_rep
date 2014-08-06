make_incidence_matrix <- function(region.vector) {
# region.vector: a vector of length = number of observations in data with 
# region.vector[i] giving the region for observation i. the number of unique
# region values should be the same as the number of cells in the map matrix
  
  N <- length(region.vector)
  R <- length(unique(region.vector))
  
  M <- mat.or.vec(N, R)
  
  for(n in 1:N){
    r <- region.vector[n]
    M[n, r] <- 1
  }
  
  return(M)
}
