#' Title Mapping of real cells into the trajectory for global alignment
#'
#' @param alignment global alignment
#' @param intTrajQuery pseudo-time scores for interpolated points along the query trajectory
#' @param realTrajQuery pseudo-time scores for real data points along the query trajectory
#' @param intTrajRef pseudo-time scores for interpolated points along the reference trajectory
#' @param realTrajRef pseudo-time scores for real data points along the reference trajectory
#'
#' @return a list consisting of meta-nodes mapping and assignment of real cells into these meta nodes
#' @export
#'
#' @examples mapRealDataGlobal(alignment, intTrajQuery, realTrajQuery, intTrajRef, realTrajRef)
mapRealDataGlobal <- function(alignment, intTrajQuery, realTrajQuery, intTrajRef, realTrajRef){
  #map meta-nodes alignments:
  metaNodesMap = metaNodesMapping(alignment, intTrajQuery, realTrajQuery, intTrajRef, realTrajRef)
  #assign cells into metaNodes:
  refAssign = metaNodesAssign(unique(metaNodesMap[,c('cellIDRef','ptRef')]), realTrajRef, option = 'global')
  queryAssign = metaNodesAssign(unique(metaNodesMap[,c('cellIDQuery','ptQuery')]), realTrajQuery, option = 'global')
  #calculate the pseudo-time of the metaNode:
  ptMetaNodeRef = calPtMetaNode(refAssign, realTrajRef)
  ptMetaNodeQuery = calPtMetaNode(queryAssign, realTrajQuery)
  #build objects:
  metaNodesPt = unique(do.call('rbind',lapply(1:nrow(metaNodesMap), function(i){
    return(data.frame(metaNodeQuery = metaNodesMap$cellIDQuery[i], metaNodeRef = metaNodesMap$cellIDRef[i],
                      ptQuery = ptMetaNodeQuery$pt[ptMetaNodeQuery$metaNode == metaNodesMap$cellIDQuery[i]],
                      ptRef = ptMetaNodeRef$pt[ptMetaNodeRef$metaNode == metaNodesMap$cellIDRef[i]],
                      ptQueryInt = metaNodesMap$ptQueryInt[i], ptRefInt = metaNodesMap$ptRefInt[i]))
  })))
  metaNodesPt$diff = abs((metaNodesPt$ptQueryInt - metaNodesPt$ptRefInt) - (metaNodesPt$ptQuery - metaNodesPt$ptRef))
  metaNodesPt = metaNodesPt[,!(colnames(metaNodesPt) %in% c('ptQueryInt','ptRefInt'))]
  return(list(metaNodesPt = metaNodesPt, refAssign = refAssign, queryAssign = queryAssign))
}

#' Title Mapping of real cells into the trajectory for local alignment
#'
#' @param alignment local alignment
#' @param intTrajQuery pseudo-time scores for interpolated points along the query trajectory
#' @param realTrajQuery pseudo-time scores for real data points along the query trajectory
#' @param intTrajRef pseudo-time scores for interpolated points along the reference trajectory
#' @param realTrajRef pseudo-time scores for real data points along the reference trajectory
#'
#' @return a list consisting of meta-nodes mapping and assignment of real cells into these meta nodes
#' @export
#'
#' @examples mapRealDataLocal(alignment, intTrajQuery, realTrajQuery, intTrajRef, realTrajRef)
mapRealDataLocal <- function(alignment, intTrajQuery, realTrajQuery, intTrajRef, realTrajRef){
  #map meta-nodes alignments:
  metaNodesMap = metaNodesMapping(alignment, intTrajQuery, realTrajQuery, intTrajRef, realTrajRef)
  #assign cells into metaNodes:
  refAssign = metaNodesAssign(unique(metaNodesMap[,c('cellIDRef','ptRef')]), realTrajRef, option = 'local')
  queryAssign = metaNodesAssign(unique(metaNodesMap[,c('cellIDQuery','ptQuery')]), realTrajQuery, option = 'local')
  #calculate the pseudo-time of the metaNode:
  ptMetaNodeRef = calPtMetaNode(refAssign, realTrajRef)
  ptMetaNodeQuery = calPtMetaNode(queryAssign, realTrajQuery)
  #build objects:
  metaNodesPt = unique(do.call('rbind',lapply(1:nrow(metaNodesMap), function(i){
    return(data.frame(metaNodeQuery = metaNodesMap$cellIDQuery[i], metaNodeRef = metaNodesMap$cellIDRef[i],
                      ptQuery = ptMetaNodeQuery$pt[ptMetaNodeQuery$metaNode == metaNodesMap$cellIDQuery[i]],
                      ptRef = ptMetaNodeRef$pt[ptMetaNodeRef$metaNode == metaNodesMap$cellIDRef[i]]))
  })))
  return(list(metaNodesPt = metaNodesPt, refAssign = refAssign, queryAssign = queryAssign))
}

#' Title Plot the alignment-based mapping between real cells from query and reference trajectories
#'
#' @param mapping mapping output from the mapRealDataGlobal or mapRealDataLocal
#'
#' @return plot
#' @export
#' @import ggplot2
#' @import reshape2
#'
#' @examples plotMapping(mapping)
plotMapping <- function(mapping){
  #plot resulting mapping:
  metaNodePt = mapping$metaNodesPt
  metaNodePt = metaNodePt[order(metaNodePt$ptQuery),]
  metaNodePt$align = 1:nrow(metaNodePt)
  metaNodePtLong = melt(metaNodePt[,c('ptQuery','ptRef','align')], id.vars = c('align'))
  metaNodePtLong = melt(metaNodePt, id.vars = c('align','metaNodeQuery','metaNodeRef'))

  ggplot(metaNodePtLong, aes(x = variable, y = value, group = align)) + geom_line(color = 'grey') + theme_bw() + geom_point() +
    coord_flip() + ggtitle('meta-nodes alignment')

}

#--------------------------------------------------------------------------
# Function: metaNodesMapping
# The function takes the alignment, the intepolated points from the query and reference
# trajectories and the real points trajectory scores.
# The function returns the mapping between meta-nodes using the alignments between interpolated points.
#--------------------------------------------------------------------------
metaNodesMapping <- function(alignment, intTrajQuery, realTrajQuery, intTrajRef, realTrajRef){
    alignMat = data.frame(queryInt = alignment$align[[1]]$index1, refInt = alignment$align[[1]]$index2)
    alignMat$ptQueryInt = intTrajQuery[alignMat$queryInt]
    alignMat$ptRefInt = intTrajRef[alignMat$refInt]
    alignMat$ptQuery = sapply(alignMat$ptQueryInt, function(ptQInt){
      return(realTrajQuery[which.min(abs(realTrajQuery - ptQInt))])})
    alignMat$ptRef = sapply(alignMat$ptRefInt, function(ptRInt){
      return(realTrajRef[which.min(abs(realTrajRef - ptRInt))])})
    alignMat$cellIDQuery = sapply(alignMat$ptQueryInt, function(ptQInt){
      return(names(realTrajQuery)[which.min(abs(realTrajQuery - ptQInt))])})
    alignMat$cellIDRef = sapply(alignMat$ptRefInt, function(ptRInt){
      return(names(realTrajRef)[which.min(abs(realTrajRef - ptRInt))])})
    alignMat = unique(alignMat)
    
    return(alignMat)
}

#--------------------------------------------------------------------------
# Function: metaNodesAssign
# The function assigns individual cells into meta nodes by NN method (from pseudo-time)
#--------------------------------------------------------------------------
metaNodesAssign <- function(metaNodeInfo, realTraj, option){
  colnames(metaNodeInfo) = c('cellID','pt')
  #Assign each cell to the closest metanode
  #for global alignment: map all the cells into meta nodes:
  if(option == 'global'){
    mapMeta = sapply(realTraj, function(x){
      return(metaNodeInfo$cellID[which.min(abs(metaNodeInfo$pt - x))])
    })
  }
  #for local option: map only the cells that are between the first and the last meta-nodes.
  if(option == 'local'){
    mapMeta = sapply(realTraj[realTraj >= min(metaNodeInfo$pt) & realTraj <= max(metaNodeInfo$pt)], function(x){
      return(metaNodeInfo$cellID[which.min(abs(metaNodeInfo$pt - x))])
    })
  }

  #per metanode, identify all the cells that were assigned to it
  metaAssignment = lapply(unique(metaNodeInfo$cellID), function(metaNode){
    return(names(mapMeta)[mapMeta == metaNode])
  })
  names(metaAssignment) = unique(metaNodeInfo$cellID)
  return(metaAssignment)
}

#--------------------------------------------------------------------------
# Function: calPtMetaNode
# The function calculates pseudo-time scores per metanode as the mean of the pseudo-time scores of the
# cells that were assigned to this metanode
#--------------------------------------------------------------------------
calPtMetaNode <- function(refAssign, realTraj){
  metaNodePt = sapply(names(refAssign), function(metaNode){
    return(mean(realTraj[metaNode]))
  })
  return(data.frame(metaNode = names(refAssign), pt = metaNodePt))
}
