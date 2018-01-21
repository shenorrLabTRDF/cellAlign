#' Title interpolation of the expression data along one trajectory
#'
#' @param expDataBatch raw single cell gene expression data (genes on rows, cells on columns)
#' @param trajCond a vector of pseudo-time scores of the data-points whose length equals to the number of samples
#' @param winSz window size of the interpolation
#' @param numPts number of desired interpolated points
#'
#' @return a list of 3: interpolatedVals, error and pseudo-time scores of the interpolated points;
#' interpolatedVals: a matrix of the interpolated points (genes*numPts)
#' error: a matrix of interpolated error for each interpolated data-point (genes*numPts)
#' traj: pseudo-time scores of the interpolated points
#'
#' @import pheatmap
#' @export
#'
interWeights <- function(expDataBatch, trajCond, winSz, numPts){
  #remove NAs from the trajectory (if exist)
  if(sum(is.na(trajCond)) > 0){
    expDataBatch = expDataBatch[,-which(is.na(trajCond))]
    trajCond = trajCond[-which(is.na(trajCond))]
  }
  #generate equally-spaced points along the trajectory:
  trajValNewPts = seq(from=min(trajCond), to=max(trajCond), length.out = numPts) #length - numPts
  ValNewPts = do.call('cbind',lapply(trajValNewPts, function(trajPt){
    dist2Others = trajCond - trajPt #length - samples number
    weightedData = exp(-(dist2Others^2)/(winSz^2))
    weightedData = weightedData/sum(weightedData)
    return(as.matrix(expDataBatch) %*% weightedData)
  }))
  #sapply on the traj - for each real data point, find the closest interpolated point:
  a = 1
  closestInt = sapply(trajCond, function(trajVal){return(which.min(abs(trajValNewPts - trajVal)))})
  #calculate the error of the smoothing function:
  errPerGene = do.call('rbind', lapply(1:nrow(expDataBatch), function(rowInd){
    return(abs(expDataBatch[rowInd,] - ValNewPts[rowInd,closestInt]))
  }))

  #interpolate the error at each interpolated point:
  errInterpolated = do.call('cbind',lapply(trajValNewPts, function(trajPt){
    dist2Others = trajCond - trajPt #length - samples number
    weightedData = exp(-(dist2Others^2)/(winSz^2))
    weightedData = weightedData/sum(weightedData)
    return(as.matrix(errPerGene) %*% weightedData)
  }))
  rownames(errInterpolated) = rownames(ValNewPts)

  return(list(interpolatedVals = ValNewPts, error = errInterpolated, traj = trajValNewPts))
}
