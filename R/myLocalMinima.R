# function findLocalMinima ------------------------------------------------
#' This function calculates the local minima of the lm matrix;
#' These minima will be the potential starting points
#' The minima calculation is performed on the rows of the lm matrix:
#' i.e. for each row, it finds the local minima, and also checks the edges.

# findLocalMinima <- function(lm, costMatrix){
#   if(min(costMatrix) > 0){return(NULL)}
#   localMin = do.call('rbind',lapply(1:nrow(lm),function(row){
#     #internal (local) minima:
#     indMin = which(diff(sign(diff(lm[row,])))==2) + 1
#     #edges (partially global) minima:
#     if(lm[row,2] > lm[row,1]){indMin = c(indMin,1)}
#     if(lm[row,ncol(lm)-1] > lm[row,ncol(lm)]){indMin = c(indMin,ncol(lm))}
#     return(data.frame(row = rep(row, length(indMin)), col = indMin, value=lm[row,indMin], costVal = costMatrix[row,indMin]))
#   }))
#   return(localMin[localMin$value < 0,c('row', 'col','costVal')])
# }

findLocalMinima <- function(lm, similarityThresh){
  localMin = do.call('rbind',lapply(1:nrow(lm),function(row){
    #internal (local) minima:
    indMin = which(diff(sign(diff(lm[row,])))==2) + 1
    #edges (partially global) minima:
    if(lm[row,2] > lm[row,1]){indMin = c(indMin,1)}
    if(lm[row,ncol(lm)-1] > lm[row,ncol(lm)]){indMin = c(indMin,ncol(lm))}
    return(data.frame(row = rep(row, length(indMin)), col = indMin, value = lm[row,indMin]))
  }))
  localMin = localMin[localMin$value < similarityThresh,]
  return(localMin)
}

# function findGlobalMinima ------------------------------------------------
#' This function calculates the global minima of the cost matrix;
#' These minima will be the potential starting points
#' The minima calculation is performed on the rows of the lm matrix:
#' i.e. for each row, it finds the local minima, and also checks the edges.

findGlobalMinima <- function(lm, costMatrix){
  if(min(costMatrix) > 0){return(NULL)}
  minInd = ind2sub(nrow(costMatrix), which.min(costMatrix))
  return(data.frame(minInd, costVal = min(costMatrix)))
}
