#' Title cluster genes based on their pseudotime shifts
#'
#' @param x interpolated scaled expression of n genes (rows) in condition 1
#' @param y interpolated scaled expression of n genes (rows) in condition 2
#' @param k number of clusters
#'
#' @return clustering result
#' @export
#'
#' @examples pseudotimeClust(x, y, k)
pseudotimeClust <- function(x, y, k = 2){
  #check input validity (that the objects have the same number of rows)
  if(nrow(x) != nrow(y)){stop('expression matrices do not have the same dimensions')}
  if(nrow(x) == 1){stop('cannot cluster one gene')}

  #find the pseudotime shift per gene:
  ptShiftPerGene = do.call('rbind',lapply(1:nrow(x), function(i){
    gAlign = globalAlign(x[i,], y[i,])
    return(gAlign$ptShift)
  }))

  #kmeans clustering of the pseudotime shifts:
  kmeansRes = kmeans(ptShiftPerGene, k)
  return(kmeansRes)
}
