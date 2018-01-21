#' Title Performs local alignment between two single-cell gene-expression patterns under automatically set threshold
#'
#' @param x query data matrix: genes on rows and cells on columns OR local cost matrix
#' @param y reference data matrix: genes on rows and cells on columns, unused if x is
#' given as a cost matrix
#' @param dist.method distance metric, used by the function "dist" in the package 'proxy'
#' @param step.pattern a stepPattern object describing the local warping steps allowed with their cost
#' @param distance.only only compute distance (no backtrack, faster)
#' @param normDist logical; scale distance matrix?
#' @param ...
#'
#' @return local alignment (as in localAlign function); plot of the picked threshold
#' @export
#'
#' @examples localAlignAutoThreshold = function(x, y)
localAlignAutoThreshold <- function(x, y = NULL,
                            dist.method = "Euclidean",
                            step.pattern=symmetric2,
                            distance.only = F,
                            normDist = T, ...){

  #define the thresholds sequence:
  thresholds = seq(0.01,0.5, by = 0.01)

  #per threshold, find the optimal alignment and calculate distance to length ratio:
  lAlign = list()
  distLocalAlign = do.call('rbind',lapply(thresholds, function(thresh){
    message(sprintf('applying local alignment threshold %s',thresh))
    alignment = localAlign(x, y, threshPercent = thresh, verbose = F)
    lAlign[[as.character(thresh)]] <<- alignment
    return(data.frame(thresh = thresh, dist = alignment$distance, len = alignment$lenAlign))
  }))
  distLocalAlign$ratio = distLocalAlign$dist/distLocalAlign$len

  #pick the most stable ratio:
  optTresh = min(distLocalAlign$thresh[distLocalAlign$len == as.numeric(names(table(distLocalAlign$len))[which.max(table(distLocalAlign$len))])])

  #plot distance to length ratio as a function of threshold:
  p = ggplot(distLocalAlign, aes(x = thresh, y = ratio)) + geom_point() + ggtitle('dist 2 length ratio under increasing threhsolds') +
    geom_vline(xintercept = optTresh, linetype="dotted", color = 'red', size = 1.5)

  return(list(localAlign = lAlign[[as.character(optTresh)]], thresh = p))
}
