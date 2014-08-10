stans_posterior_medians <- function(stan) {
  M <- monitor(stan)
  M <- M[,"50%"]
  M
} 

retrieve_stans_posterior <- function(stan, fun) {
  
  M <- monitor(stan)
  easy_funs <- colnames(M)
  fun_match <- match(fun, possible_funs)
  
  if (!is.na(fun_match)) {
    M <- M[, fun_match]  
  }
  
  if (fun == "ci90") {
    M <- M[, c("2.5%", "97.5%")]
  }
  
  if (fun == "median") {
    M <- M[, "50%"]
  }
  
  return(M)
} 