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