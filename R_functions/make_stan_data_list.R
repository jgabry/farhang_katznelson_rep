make_stan_data_list <- function(df, map.dims, region.varname, region.vector, Xvarlist, Zvarlist, Yname) {
  # either specify region.varname or region.vector
    # region.varname <- character string naming the "region" variable in df
    # region.vector <- a vector to use that is outside of df
  
  
# checks
  lmd       <- length(map.dims)
  df_is_df  <- is.data.frame(df)
  Y_char    <- is.character(Yname)
  X_char    <- is.character(Xvarlist)
  Z_char    <- is.character(Zvarlist)
  no_rvarname <- missing(region.varname)
  no_rvector <- missing(region.vector)
  
  if (!X_char) {
    stop("Xvarlist must be of type 'character'")
  }
  if (!Z_char) {
    stop("Zvarlist must be of type 'character'")
  }
  if (!df_is_df) {
    stop(paste("df is type", class(df), "but it should be a data frame"))
  }
  if (lmd != 2) {
    stop(paste("length(map.dims) is", lmd, "but it should be 2."))
  }
  if (no_rvarname & no_rvector) {
    stop("Either region.varname or region.vector must be specified")
  }
  if (!Y_char) {
    Yname <- deparse(substitute(Yname)) 
  }
  if (no_rvarname & is.character(region.vector)) {
    region.vector <- get(region.vector)
  }
  if (no_rvector){
    if (!is.character(region.varname)) {
      region.varname <- deparse(substitute(region.varname))
    }
    region.vector <- df[, region.varname]
  }
  
  
  
  Y <- df[, Yname]
  A <- make_adjacency_matrix(map.dims)
  D <- diag(rowSums(A))
  M <- make_incidence_matrix(region.vector)
  X <- make_design_matrix_array(df, Xvarlist, M)
  Z <- as.matrix(df[, Zvarlist])

  R_is_ok <- nrow(A) == ncol(M) & nrow(A) == length(diag(D)) & nrow(A) == length(unique(region.vector))
  stopifnot(R_is_ok)
  
  data_list <- list(Y = Y,
                    Z = Z,
                    X = X,
                    M = M,
                    A = A, 
                    A_N = make_sparse(A)[["A_N"]], 
                    A1 = make_sparse(A)[["A1"]], 
                    A2 = make_sparse(A)[["A2"]],
                    d = diag(D),
                    R = nrow(A),
                    N = nrow(df),
                    J = dim(X)[1] + 1,
                    C = ncol(Z)
                    )

  return(data_list)
}