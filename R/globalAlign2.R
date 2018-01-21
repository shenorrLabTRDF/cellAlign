#' Title Performs global alignment between two single-cell gene-expression patterns
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
#' @return a list containing: local cost, cost, and step matrices, global alignment
#' @export
#' @import dtw
#'
#' @examples globalAlign(x,y)
#' @useDynLib cellAlign computeCM
globalAlign = function(x, y = NULL,
                       dist.method = "Euclidean",
                       step.pattern=symmetric2,
                       distance.only = F,
                       normDist = T,
                       verbose = T,
                       ... ) {

  #check input validity:
  if(is.null(y) & !is.matrix(x)){stop('either two expression patterns or one dissimilarity matrix should be used')}
  if(is.null(y) & is.null(x)){stop('both arguments equal to NA')}

  #calculate dissimilarity matrix:
  message('calculate dissimilarity matrix')
  if(is.null(y)){lm = x
  }else{lm = calcDistMat(x, y, dist.method)}

  #calculate the number of genes used for alignment:
  if(is.null(y)){numGenes = NA
  }else if(is.null(dim(x))){numGenes = 1
  }else{numGenes = nrow(x)}

  #normalized distance:
  if(normDist){lm = (lm-min(lm))/(max(lm) - min(lm))}

  #if only the distance matrix is required, return it:
  if(distance.only){return(lm)}

  #build cost and step matrices:
  message('calculate cost and step matrices')
  gcm = globalCostMatrix(lm, step.matrix = step.pattern)

  gcm$jmin = ncol(lm)

  #backtrack:
  message('backtracking')
  if(!distance.only) {
    gcm$align = list()
    gcm$align[[1]] = backtrack2(gcm)
  }

  #pseudotime shift:
  ptQuery = seq(0,1,length.out = nrow(lm))
  ptRef = seq(0,1,length.out = ncol(lm))
  gcm$ptShift = sapply(1:length(ptRef), function(i){
    return(mean(ptQuery[gcm$align[[1]]$index1[gcm$align[[1]]$index2 == i]]) - ptRef[i])
  })

  #build resulting object:
  #distance:
  gcm$distance = gcm$costMatrix[nrow(gcm$costMatrix),ncol(gcm$costMatrix)]

  #normalized distance:
  norm = attr(step.pattern,"norm")
  if(is.na(norm)) {
    gcm$normalizedDistance = NA
  } else if(norm == "N+M") {
    gcm$normalizedDistance = gcm$distance/(nrow(lm) + ncol(lm))
  } else if(norm == "N") {
    gcm$normalizedDistance = gcm$distance/nrow(lm)
  } else if(norm == "M") {
    gcm$normalizedDistance = gcm$distance/ncol(lm)
  }

  #other parameters:
  gcm$align.method = 'global'
  gcm$localCostMatrix = lm
  if(!is.na(numGenes)){gcm$normDistModule = gcm$distance/sqrt(numGenes)}

  return(gcm)
}
