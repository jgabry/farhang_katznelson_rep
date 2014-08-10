#' Construct an array of design matrices to pass as data to \code{stan} (\pkg{rstan})
#'
#'
#' @param df The data frame containing the variables in \code{varlist}.
#' @param varlist A character vector naming variables in \code{df}.
#' @param M An incidence matrix with \code{nrow(M) = nrow(df)}. Only required 
#' if \code{region.vector} is not specified. See 'Details'. 
#' @param region.vector A numeric vector with length \code{N = \code{nrow(df)}, and where
#' \code{region.vector[i]} gives the region number for observation \code{i}. The number of 
#' unique region values (henceforth \code{K}) should be the same as the number of cells in 
#' the map. \code{region.vector} is only required if an incidence matrix \code{M} is not 
#' specified. See 'Details'.
#' @return An array \code{X} with (\code{dim = c(length(varlist), nrow(M), ncol(M))}). That is,
#' if \code{L = length(varlist)} then \code{X} is an array of \code{L} matrices of \code{dim(M)}.
#' The matrices in \code{X} are essentially incidence matrices with the 1s replaced by values
#' of the variables in \code{varlist}. This is accomplished via element-wise multiplication.    
#' @details The user must specify either \code{M} or \code{region.vector}. If \code{M} is not
#' specified then it will be constructed via an internal call to \code{make_incidence_matrix(region.vector)}.
#' @author Jonah Gabry <jsg2201@@columbia.edu>
#' @export
#' @examples
#' \dontrun{
#' # Suppose we have a data frame df which contains five variables V1, V2, ... V5. 
#' # Then the array \code{X} defined by
#' X <- make_design_matrix_array(df, varlist = c("V2, V3"), M)
#' # contains 2 matrices with dimensions dim(M) where corresponding to the element-wise 
#' # products M*V2 and M*V3.
#' }
#' 


make_design_matrix_array <- function(df, varlist, M, region.vector) {
# varlist = character vector naming variables in the data
# M = incidence matrix (optional) 
# region.vector = same as region.vector arg for make_incidence_matrix function.
# df = data frame in which to look for vars in varlist
  
# If M is not specified then region.vector MUST be specified
# If M is specified then region.vector is NOT needed
  
  no_M <- missing(M)
  no_region <- missing(region.vector)
  
  if(no_M & no_region) {
    stop("You must specify either M or region.vector")
  }
  if(!no_M & !no_region) {
    M_check <- make_incidence_matrix(region.vector)
    bad <- !identical(M_check, M)
    if(bad) {
      stop("The incidence matrix (M) that you specified is not 
           equivalent to the incidence matrix corresponding to 
           region.vector.")
    }
  }
  if (no_M & !no_region) {
    M <- make_incidence_matrix(region.vector)
  }
  
  nrow_bad <- nrow(M) != nrow(Data)
  if(nrow_bad){
    stop("df and M should have the same number of rows")
  }
    
  
  temp <- df[, varlist]
  ll <- length(varlist)
  R <- ncol(M)
  N <- nrow(M)
  
  X <- array(NA, dim = c(ll, N, R))
  
  for(l in 1:ll) {
    X[l,,] <- M*temp[, l] 
  }
  
  return(X)
}