#' Title Calculate dissimilarity matrix (local cost matrix)
#'
#' @param x single cell gene expression data: query trajectory (genes on rows, cells on columns)
#' @param y single cell gene expression data: reference trajectory (genes on rows, cells on columns)
#' @param dist.method distance metric (see: proxy)
#'
#' @return dissimilarity matrix (local cost matrix) between x and y
#' @export
#'
#' @examples calcDistMat(x,y)
calcDistMat <- function(x, y, dist.method){
  if(is.null(dim(x))){lm = proxy::dist(x,y,method=dist.method)}
  else{
    x = as.matrix(t(x))
    y = as.matrix(t(y))
    lm = proxy::dist(x,y,method=dist.method)
  }
  return(lm)
}
