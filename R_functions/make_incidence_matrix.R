#' Construct an incidence matrix from a vector of region codes. 
#'
#'
#' @param region.vector A numeric vector with length \code{N} = number of observations in 
#' the data with \code{region.vector[i]} giving the region number for observation \code{i}. 
#' The number of unique region values (henceforth \code{K}) should be the same as the number 
#' of cells in the map. 
#' @return An \code{NxK} matrix \code{M}, with \code{m_ij = 1} iff observation \code{i} is
#' from region \code{j} and \code{m_ij = 0} otherwise. 
#' @author Jonah Gabry <jsg2201@@columbia.edu>
#' @export
#' @examples
#' # A simple region.vector
#' region <- c(1,1,2,3,1,3,2)
#' 
#' # The corresponding incidence matrix M will then have length(region) = 7 rows 
#' # and length(unique(region)) = 3 columns 
#' M <- make_incidence_matrix(region) 
#' dim(M)


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
