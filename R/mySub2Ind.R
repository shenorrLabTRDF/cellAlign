#' Title Convert a row-column vector to a linear index
#'
#' @param rowNum nuber of rows in the data matrix
#' @param rows a vector containing the rows of the converted items
#' @param cols a vector containing the columns of the converted items
#'
#' @return a vector of the linear indices
#' @export
#'
#' @examples sub2ind(rowNum, rows, cols)
sub2ind = function(rowNum, rows, cols){
  indices = sapply(1:length(rows), function(i){
    return((cols[i] - 1)*rowNum + rows[i])
  })
  return(indices)
}

#' Title Convert linear indices to a row-columns matrix
#'
#' @param rowNum number of rows in the desired matrix
#' @param ind a vector of linear indices
#'
#' @return a data.frame with the row and columns odf the linear indices
#' @export
#'
#' @examples ind2sub(rowNum, ind)
ind2sub = function(rowNum, ind){
  indices = do.call('rbind', lapply(ind, function(i){
    return(data.frame(row = ((i-1) %% rowNum) + 1, col = floor((ind-1) / rowNum) + 1))
  }))
  return(indices)
}
