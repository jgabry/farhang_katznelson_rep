make_adjacency_matrix <- function(map.dims){
  # map.dims = vector giving dimensions of map matrix (e.g. dims = c(3, 4), for 3 rows, 4 cols)
  dim <- prod(map.dims)
  A <- mat.or.vec(nr = dim, nc = dim)
  neighs <- sapply(1:dim, function(i) find_neighbors(map.dims, i))
  for(i in 1:dim){
    A[i, neighs[[i]] ] <- 1
  }
  return(A)
}