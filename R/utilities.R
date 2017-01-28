# Tree-processing functions

NULL

#' Extract information for pruning a tree used as cache in poumm likelihood
#' calculation
#' 
#' @param tree a phylo object
#'   
#' @details This method should only be called if calculating poumm likelihood
#'   with impl='R5'.
#' @return a list of objects
pruneTree <- function(tree) {
  N <- length(tree$tip.label)                # number of tips
  M <- length(unique(as.vector(tree$edge)))  # number of all nodes 
  endingAt <- order(rbind(tree$edge, c(0, N+1))[, 2])
  
  edge <- tree$edge
  mode(edge) <- "integer"
  
  nonVisitedChildren <- rep(0, M)
  ee1 <- edge[, 1]
  while(length(ee1)) {
    matchp <- match((N+1):M, ee1)
    matchp <- matchp[!is.na(matchp)]
    nonVisitedChildren[ee1[matchp]] <- nonVisitedChildren[ee1[matchp]] + 1
    ee1 <- ee1[-matchp]
  }
  
  # start from the edges leading to tips
  nodesVector <- as.integer(c())
  nodesIndex <- as.integer(c(0))
  
  unVector <- as.integer(c())
  unIndex <- as.integer(c(0))
  
  nodes <- as.integer(1:N)
  
  while(nodes[1] != N+1) {
    nodesIndex <- c(nodesIndex, nodesIndex[length(nodesIndex)]+length(nodes))
    nodesVector <- c(nodesVector, nodes)
    
    es <- endingAt[nodes]
    nodes <- as.integer(c())
    edgeEnds <- edge[es, 2]
    
    #update parent pifs
    while(length(es)>0) {
      un <- match(unique(edge[es, 1]), edge[es, 1])
      unIndex <- c(unIndex, unIndex[length(unIndex)]+length(un))
      unVector <- c(unVector, un)
      #pif[edge[es[un], 1], ] <- pif[edge[es[un], 1], ] + pif[edge[es[un], 2], ]
      nonVisitedChildren[edge[es[un], 1]] <- nonVisitedChildren[edge[es[un], 1]] - 1
      nodes <- c(nodes, edge[es[un][nonVisitedChildren[edge[es[un], 1]] == 0], 1])
      es <- es[-un]
    }
  }
  list(M=M, edge = edge, endingAt=endingAt, 
       nodesVector=nodesVector, nodesIndex=nodesIndex, nLevels=length(nodesIndex)-1, 
       unVector=unVector, unIndex=unIndex)
}


#' Node indices of the direct descendants of n in the phylogeny tree.
#' 
#' @param tree an object of class phylo
#' @param n an index of a node (root, internal or tip) in tree
#' @return An integer vector.
chld <- function(tree, n) {
  as.integer(tree$edge[tree$edge[, 1]==n, 2])
}

#' Edge indices of the edges in tree starting from n
#' @param tree an object of class phylo
#' @param n an index of a node (root, internal or tip) in tree
#' @return An integer vector.
edgesFrom <- function(tree, n) {
  which(tree$edge[, 1]==n)
}

#' Calculate the time from the root to each node of the tree
#' 
#' @param tree An object of class phylo.
#' @param tipsOnly Logical indicating whether the returned results should be
#'   truncated only to the tips of the tree.
#' @return A vector of size the number of nodes in the tree (tips, root, 
#'   internal) containing the time from the root to the corresponding node in 
#'   the tree.
#' @importFrom stats reorder
#' @import ape
#' @export
nodeTimes <- function(tree, tipsOnly=FALSE) {
  rtree <- reorder(tree, 'postorder')
  es <- rtree$edge[dim(rtree$edge)[1]:1, ]
  nEdges <- dim(es)[1]
  ts <- rev(rtree$edge.length)
  nodeTimes <- rep(0, length(rtree$tip.label)+rtree$Nnode)
  for(e in 1:nEdges) 
    nodeTimes[es[e, 2]] <- nodeTimes[es[e, 1]]+ts[e]
  if(tipsOnly) {
    nodeTimes[1:length(tree$tip.label)]
  } else {
    nodeTimes
  }
}
