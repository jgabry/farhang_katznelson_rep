make_sparse <- function(A) {
# A = adjancency matrix
  
  ### sparse matrix stuff ###
  A_sparse <- which(A == 1, arr.ind=TRUE)
  # remove duplicates (because matrix is symmetric)
  A_sparse <- A_sparse[A_sparse[,1] < A_sparse[,2],]
  A_N <- dim(A_sparse)[1]
  A1 <- A_sparse[,1]
  A2 <- A_sparse[,2]
  
  out <- list(A_N = A_N, A1 = A1, A2 = A2)
  return(out)
}

