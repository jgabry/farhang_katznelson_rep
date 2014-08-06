find_neighbors <- function (map.dims, cell, row = NULL, col = NULL) {
  # map.dims = vector giving dimensions of map matrix (e.g. dims = c(3, 4), for 3 rows, 4 cols)
  # either specify: cell = cell number, or row & col giving the row and column of the cell
  
  
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
#   vec_c <- c(c, c + 1, c + 1, c + 1, c, c - 1, c - 1, c - 1)
  vec_c <- c(c, rep(c + 1, 3), c, rep(c - 1, 3))
#   vec_r <- c(r + 1, r + 1, r, r - 1, r - 1, r - 1, r, r + 1)
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
#   cell.n <- function(r, c) return(c * nr - nr + r)
#   return(cell.n(r = vec_r, c = vec_c))
}

# Examples
# find_neighbors(dims = c(3,8), cell = 21)
# find_neighbors(dims = c(3,8), row = 3, col = 7)
