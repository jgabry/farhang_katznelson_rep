#' Construct the adjacency matrix for any \code{NxK} map. 
#'
#'
#' @param map.dims A numeric vector with two elements specifying the dimensions 
#' of the map. For example, \code{map.dims = c(R,C)} for R rows and C columns. 
#' @return A symmetric \code{R*C} by \code{R*C} matrix \code{A}, with \code{a_ij = 1} iff
#' cells \code{i} and \code{j} are neighbors and \code{a_ij = 0} otherwise. 
#' @details The neighbors of a cell \code{c} are found by an internal call to \code{find_neighbors}.
#' @author Jonah Gabry <jsg2201@@columbia.edu>
#' @seealso \code{find_neighbors}
#' @export
#' @examples
#' A <- make_adjacency_matrix(map.dims = c(3, 4)) # A is then a symmetric 12x12 matrix. 


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