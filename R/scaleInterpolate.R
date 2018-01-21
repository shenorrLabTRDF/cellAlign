#' Title Scale interpolated data to [0,1]
#'
#' @param interDataConds interpolated data; output from the function: interWeights
#'
#' @return a list consisting from the interpolated-scaled values, interpolated-scaled error and pseudo-time scores of the interpolated points
#' @export
#'
#' @examples scaleInterpolate(interDataConds)
scaleInterpolate <- function(interDataConds){
  #scale the interpolated error and data; put them both in the same matrix:
  scaledDataError = do.call('rbind', lapply(1:nrow(interDataConds$interpolatedVals), function(rowInd){
    minVal = min(interDataConds$interpolatedVals[rowInd,])
    maxVal = max(interDataConds$interpolatedVals[rowInd,])
    scaledData = (interDataConds$interpolatedVals[rowInd,] - minVal)/(maxVal - minVal)
    scaledError = interDataConds$error[rowInd,]/(maxVal - minVal)
    return(rbind(scaledData, scaledError))
  }))

  #divide the total matrix into the scaled data and scaled error
  scaledData = scaledDataError[seq(1, nrow(scaledDataError), by=2),]
  scaledError = scaledDataError[seq(2, nrow(scaledDataError), by=2),]

  rownames(scaledData) = rownames(interDataConds$interpolatedVals)
  rownames(scaledError) = rownames(interDataConds$error)

  return(list(scaledData = scaledData, scaledError = scaledError, traj = interDataConds$traj))
}
