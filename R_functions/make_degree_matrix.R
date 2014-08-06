make_degree_matrix <- function(A, map.dims = NULL) {
  # A = adjacency matrix
  # map.dims = vector giving dimensions of map matrix (only needed if A is not given)  
  if (missing(A)) {
    A <- make_adjacency_matrix(map.dims)
  }
  
  Nneighs <- rowSums(A)
  D <- diag(Nneighs)
  
  return(D)
}