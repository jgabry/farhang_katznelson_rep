#' Find the indices of all neighbors for a cell in a map of any dimensions
#'
#'
#' @param map.dims A numeric vector with two elements specifying the dimensions 
#' of the map. For example, \code{map.dims = c(R,C)} for R rows and C columns. 
#' @param cell An integer giving the index number of the cell for which to find
#' the neighbors. For a \code{RxC} map (a map with \code{R*C} cells), the cell 
#' indices are the set \code{\{1,2, ..., R*C\}}, where indices \code{1} through \code{R} 
#' refer to the cells occupying column 1, indices \code{R+1} through \code{2R} refer to 
#' the cells occupying column 2, etc. See 'Examples' below. 
#' @param row, col In place of specifying a cell index number in the \code{cell} argument, 
#' the the same cell can be referred to by specifying its row and column in the map.
#' @return A numeric vector containing the cell indices corresponding to the neighbors 
#' of \code{cell}. Note that the returned vector is identical even if \code{row} and \code{col} are 
#' specified in place of \code{cell}. That is, the neighbors are referred to by their indices 
#' and not their row/column position in the map.
#' @details The \code{find_neighbors} function considers the neighbors of a cell \code{c} to
#' be the set of cells within one position (horizontally, vertically, or diagonally) of 
#' \code{c} in the map. 
#' @note The code for \code{find_neighbors} is based on a similar function \code{neigh_cell} (\pkg{spdynmod}).
#' @author Jonah Gabry <jsg2201@@columbia.edu>
#' @seealso \code{make_adjacency_matrix}
#' @export
#' @examples
#' ## Without loss of generality, the examples below use a 3 by 4 map ##
#' 
#' # Check the indices for the cells
#' map.indices <- matrix(1:12, nrow = 3, ncol = 4)
#' map.indices
#' 
#' # Find the neighbors for the cell in row 2, column 3 (i.e. index number 8)
# find_neighbors(map.dims = c(3,4), cell = 8) # using the cell argument
# find_neighbors(map.dims = c(3,4), row = 2, col = 3) # same output but using the row,col arguments.

find_neighbors <- function (map.dims, cell, row = NULL, col = NULL) {

  no_cell <- missing(cell)
  no_row <- missing(row)
  no_col <- missing(col)
  
  if (no_cell & (no_row | no_col)) {
    stop("If cell is not specified then you must specify both row and col")
  }

  nr <- map.dims[1]
  nc <- map.dims[2]
  
# get cell number if user inputs only row and col  
  if (no_cell) {
    cell <- nr*(col-1) + row
  }
  
  c <- ceiling(cell/nr)
  r <- (cell + nr - c * nr)
  vec_c <- c(c, rep(c + 1, 3), c, rep(c - 1, 3))
  vec_r <- c(rep(r + 1, 2), r, rep(r - 1, 3), r, r + 1)
  wc_neg <- which(vec_c < 1 | vec_c > nc)
  if (length(wc_neg) > 0) {
    vec_c <- vec_c[-wc_neg]
    vec_r <- vec_r[-wc_neg]
  }
  wr_neg <- which(vec_r < 1 | vec_r > nr)
  if (length(wr_neg) > 0) {
    vec_c <- vec_c[-wr_neg]
    vec_r <- vec_r[-wr_neg]
  }
  
  neighs <- vec_c*nr - nr + vec_r
  return(sort(neighs))
}
