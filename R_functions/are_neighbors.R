#' Quickly test if two cells in a map are neighbors. 
#'
#'
#' @param map.dims A numeric vector with two elements specifying the dimensions 
#' of the map. For example, \code{map.dims = c(R,C)} for R rows and C columns. 
#' @param cell1, cell2 The indices for the cells to test. 
#' @return A logical value, which is \code{TRUE} iff \code{cell1} and \code{cell2} are
#' neighbors. 
#' @details The neighbors of \code{cell1} are found by an internal call to \code{find_neighbors}.
#' If \code{cell2} is found among the neighbors of \code{cell1} then \code{are_neighbors} will
#' return a value of \code{TRUE}.
#' @author Jonah Gabry <jsg2201@@columbia.edu>
#' @seealso \code{find_neighbors}
#' @export
#' @examples
#' are_neighbors(map.dims = c(3,4), cell1 = 2, cell2 = 11)  #  FALSE
#' are_neighbors(map.dims = c(3,4), cell1 = 2, cell2 = 5)   #  TRUE


are_neighbors <- function(map.dims, cell1, cell2) {
  answer <- cell2 %in% find_neighbors(map.dims, cell = cell1)
  return(answer)
}