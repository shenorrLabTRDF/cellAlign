#' Title Plot the alignment over the dissimilarity matrix
#'
#' @param alignment output from the localAlign/globalAlign functions
#'
#' @return NA
#' @export
#' @import pheatmap
#'
#' @examples plotAlign(globalAlignment)
plotAlign <- function(alignment){
  costMat = alignment$localCostMatrix
  costMat = t(apply(costMat,1,function(x){return(as.numeric(x))}))
  linearInd = sub2ind(nrow(costMat), alignment$align[[1]]$index1, alignment$align[[1]]$index2)
  costMat[linearInd] = NA
  costMat = data.frame(costMat, row.names=1:nrow(costMat))
  colnames(costMat) = 1:ncol(costMat)
  #for global alignment, where there is a pseudotime shift vector:
  if(!is.null(alignment$ptShift)){
    annotCols = data.frame(ptShift = abs(alignment$ptShift), sign = factor(sign(alignment$ptShift)),row.names = colnames(costMat))
    pheatmap(costMat, cluster_cols = F, cluster_rows=F, border_color = NA,
             main = 'alignment plot',
             show_rownames = F, show_colnames = F, annotation_col = annotCols)
  }else{
    pheatmap(costMat, cluster_cols = F, cluster_rows=F, border_color = NA,
             main = 'alignment plot', show_rownames = F, show_colnames = F)
  }

  return(NA)
}
