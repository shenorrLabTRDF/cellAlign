#' Title Performs local alignment between two single-cell gene-expression patterns
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
#' @return a list containing: local cost, cost, and step matrices, local alignment
#' @export
#' @import dtw
#'
#' @examples localAlign = function(x, y)
localAlign = function(x, y = NULL,
                      dist.method = "Euclidean",
                      step.pattern=symmetric2,
                      distance.only = F,
                      normDist = T,
                      threshPercent = 0.05,
                      userDefinedThresh = NA,
                      verbose = T,
                      ... ) {

  #check input validity:
  if(is.null(y) & !is.matrix(x)){stop('either two expression patterns or one dissimilarity matrix should be used')}
  if(is.null(y) & is.null(x)){stop('both arguments are null')}

  #calculate dissimilarity matrix:
  if(verbose){message('calculate dissimilarity matrix')}
  if(is.null(y)){lm = x
  }else{lm = calcDistMat(x, y, dist.method)}

  #normalized distance:
  if(normDist){lm = (lm-min(lm))/(max(lm) - min(lm))}

  #if only the distance matrix is required, return it:
  if(distance.only){return(lm)}

  #calculate similarity threshold based on the precentile given by the user:
  if(is.na(userDefinedThresh)){similarityThresh = threshPercent*(max(lm) - min(lm))
  }else{similarityThresh = userDefinedThresh}

  #build cost and step matrices:
  if(verbose){message('calculate cost and step matrices')}
  gcm = globalCostMatrix(lm, step.matrix = step.pattern)

  gcm$jmin = ncol(lm)

  #identify potential starting cells:
  if(verbose){message('find potential starting points')}
  localMin = findLocalMinima(lm, similarityThresh)
  localMin = localMin[localMin$row > 1 & localMin$col > 1,]

  #apply backtracking from each potential starting cells:
  if(verbose){message('find the optimal alignment')}
  alignIntervals = lapply(1:nrow(localMin), function(i){
    #build gcm matrix per potential starting point:
    gcmLoc = list(costMatrix = gcm$costMatrix[1:localMin$row[i],1:localMin$col[i]],
                  directionMatrix = gcm$directionMatrix[1:localMin$row[i],1:localMin$col[i]],
                  stepPattern = gcm$stepPattern,
                  jmin = localMin$col[i])

    #backtrack:
    backTr = backtrack2(gcmLoc)

    #prun the alignment elements that exhibit lm > similarityThresh:
    alignElements = sub2ind(nrow(lm), backTr$index1, backTr$index2)
    if(sum(lm[alignElements] > similarityThresh) > 0){
      removeInd = max(which(lm[alignElements] > similarityThresh))
      backTr = lapply(backTr, function(el){return(el[(removeInd + 1):length(el)])})}

    return(backTr)
  })

  #get the length and the distance of the found local alignments:
  lenAlign = sapply(alignIntervals, function(alignment){return(length(alignment$index1))})
  distAlign = mapply(function(alignment, len){
    dist = gcm$costMatrix[alignment$index1[length(alignment$index1)], alignment$index2[length(alignment$index1)]] -
      gcm$costMatrix[alignment$index1[1], alignment$index2[1]]
    return(dist/len)}, alignIntervals, lenAlign)

  #the alignment of the minimal length will be selected as local alignment:
  gcm$align = list()
  chosenAlign = order(lenAlign, -distAlign, decreasing = T)[1]
  gcm$align[[1]] = alignIntervals[[chosenAlign]]

  #other parameters:
  gcm$align.method = 'local'
  gcm$localCostMatrix = lm
  gcm$lenAlign = length(gcm$align[[1]]$index1)
  gcm$distance = gcm$costMatrix[gcm$align[[1]]$index1[gcm$lenAlign], gcm$align[[1]]$index2[gcm$lenAlign]] -
    gcm$costMatrix[gcm$align[[1]]$index1[1], gcm$align[[1]]$index2[1]]

  return(gcm)
}
