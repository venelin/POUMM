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
  
  nonVisitedChildren <- rep(0, M)
  ee1 <- edge[, 1]
  while(length(ee1)) {
    matchp <- match((N+1):M, ee1)
    matchp <- matchp[!is.na(matchp)]
    nonVisitedChildren[ee1[matchp]] <- nonVisitedChildren[ee1[matchp]] + 1
    ee1 <- ee1[-matchp]
  }
  
  # start from the edges leading to tips
  nodesVector <- c()
  nodesIndex <- c(0)
  
  unVector <- c()
  unIndex <- c(0)
  
  nodes <- 1:N
  
  while(nodes[1] != N+1) {
    nodesIndex <- c(nodesIndex, nodesIndex[length(nodesIndex)]+length(nodes))
    nodesVector <- c(nodesVector, nodes)
    
    es <- endingAt[nodes]
    nodes <- c()
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
  list(M=M, endingAt=endingAt, 
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

#' Make the ending edges of a tree longer/shorter by addLength
#' @param tree a phylo object
#' @param addLength numeric vector of length 1 or length(tree$tip.label)
#' @return the tree with modified edge lengths for the edges leading to tips
extendTipEdges <- function(tree, addLength) {
  tipEdgeIdx <- match(1:length(tree$tip.label), tree$edge[, 2])
  tree$edge.length[tipEdgeIdx] <- tree$edge.length[tipEdgeIdx]+addLength
  tree
}

#' Cut a tree to the largest ultrametric subtree not exceeding maxNTips tips
#' @param tree an object of class phylo
#' @param maxNTips integer indicating the desired number of tips in the resulting
#' ultrametric tree.
#' @return a phylo object
#' @description The tree gets cut at the closest tip from the root or as soon as
#' maxNTips lineages are counted in a breadth first traversal of the tree. 
#' @note no tips are dropped from the tree
cutUltrametric <- function(tree, maxNTips) {
  nTips <- length(tree$tip.label)
  if(maxNTips < 1 | maxNTips>nTips) {
    stop('N must be in the range [1, number of tips in the tree].')
  }
  
  ndTimes <- nodeTimes(tree)
  nNodes <- length(ndTimes)
  
  newTip <- c()
  newNode <- c(nTips+1)
  newEdge <- matrix(0, nrow=0, ncol=2)
  newEdgeLength <- c()
  
  # outgoing edges from the root
  es <- edgesFrom(tree, nTips+1)
  
  while(TRUE) {
    countLineages <- length(es)
    
    spikingEdge <- which.min(ndTimes[tree$edge[es, 2]])
    spikeParent <- tree$edge[es, 1][spikingEdge]
    spike <- tree$edge[es, 2][spikingEdge]
    
    #print(spike)
    if(countLineages >= maxNTips | length(edgesFrom(tree, spike)) == 0) {
      # cut tree at spike
      newTip <- tree$edge[es, 2]
      esLengths <- ndTimes[spike]-ndTimes[tree$edge[es, 1]]
      newEdge <- rbind(newEdge, tree$edge[es, ])
      newEdgeLength <- c(newEdgeLength, esLengths)
      newLabels <- c(rev(newEdge[,2])[1:length(newTip)], nTips+1, rev(newEdge[, 2])[-(1:length(newTip))])
      mapNewNodeIds <- rep(NA, max(newLabels))
      mapNewNodeIds[newLabels] <- 1:length(newLabels)
      edge <- apply(newEdge, c(1, 2), function(.) mapNewNodeIds[.])
      res <- list(edge=edge, Nnode=nrow(edge)-length(newTip)+1, edge.length=newEdgeLength, tip.label=as.character(newLabels[1:length(newTip)]), node.label=as.character(newLabels[-(1:length(newTip))]))
      class(res) <- 'phylo'
      return(res)
    } else {
      # add spike to new tree
      newNode <- c(newNode, spike)
      newEdge <- rbind(newEdge, c(spikeParent, spike))
      newEdgeLength <- c(newEdgeLength, tree$edge.length[es][spikingEdge])
      
      # remove edge leading to spike from es
      es <- es[-spikingEdge]
      # add edges descending from spike to es
      es <- c(es, edgesFrom(tree, spike))
    }
  }
}

#' get a sample list of possible transmission pair assignments
#' @param tree a phylo object with N tips.
#' @param g numerical vector of length the number of nodes in the tree: the
#'   known viral genetic contributions at the tips and internal nodes.
#' @param e numerical vector of length N: the known environmental deviations at 
#'   the tips.
#' @param nSamples : number of samples. Each sample corresponds to an assignment
#'   of the internal nodes of the tree with tip-labels and the correalation of 
#'   the formed set of source-recipient phenotypes.
#' @param sampTime : a character indicating how the phenotypes for the pairs 
#'   shall be formed. Currently two possible values are supported : 
#'   eachAfterGettingInfected - the source phenotype is measured at the moment
#'   of infection from its own source, the recipient phenotype is measured at
#'   the moment it got infected from the source; bothAtInfection - the source
#'   and the recipient phenotypes are taken at the moment of infection from
#'   source to recipient. atTips - the source and the recipient phenotypes are
#'   measured at their corresponding tips in the tree.
#' @return a list of nSamples matrices, each matrix having tree columns
#'   representing source phenotype, recipient phenotype and time of the
#'   infection from the root
sampleTransmissionPairs <- function(tree, g, e, nSamples, 
                                    sampTime=c('eachAfterGettingInfected', 'bothAtInfection', 'atTips')) {
  N <- length(tree$tip.label)
  gAfterInf <- g
  gAfterInf[tree$edge[, 2]] <- g[tree$edge[, 1]]
  gAfterInf[N+1] <- g[N+1]
  infectTimes <- nodeTimes(tree)[tree$edge[, 1]]
  
  lapply(1:nSamples, function(i) {
    id <- sampleNodeIds(tree)
    edge <- tree$edge
    
    es <- matrix(e[id[edge]], ncol=2)
    if(length(grep('eachAfter', sampTime[1]))>0) {
      gs <- matrix(gAfterInf[edge], ncol=2)
    } else if(length(grep('bothAt', sampTime[1]))>0) {
      gs <- matrix(g[edge[, 1]], nrow=nrow(edge), ncol=2)
    } else if(length(grep('atTips', sampTime[1]))>0) {
      gs <- matrix(g[id[edge]], ncol=2)
    }
    pp <- cbind(gSource=gs[,1], eSource=es[,1], gRecip=gs[,2], eRecip=es[,2], infectTimes)
    edge <- matrix(id[edge], ncol=2)
    pp <- pp[edge[, 1]!=edge[, 2], ]
    pp
  })  
}

#' Sample from the distribution of Pearson correlation indices
#' between possible source-recipient couples given an infectious tree.
#' @param tree a phylo object with N tips.
#' @param g numerical vector of length the number of nodes in the tree: the known viral genetic 
#'    contributions at the tips and internal nodes.
#' @param e numerical vector of length N: the known environmental deviations at
#'    the tips. 
#' @param nSamples : number of samples. Each sample corresponds to an assignment
#'    of the internal nodes of the tree with tip-labels and the correalation of
#'    the formed set of source-recipient phenotypes.
#' @param sampTime : a character indicating how the phenotypes for the pairs 
#'    shall be formed. Currently two possible values are supported : 
#'      justAfterInfected - phenotypes after infection;
#'      bothAtInfection - the source and the recipient phenotypes are taken at 
#'        the moment of infection.
#' @param reg logical indicating if regression slopes instead of Pearson correlation indices 
#' should be returned, FALSE by default.
#' @importFrom stats cor
#' @return a numerical vector of length nSamples
sampleCorrsPairs <- function(tree, g, e, nSamples=1000, 
                             sampTime=c('justAfterInfected', 'bothAtInfection'), 
                             reg=F) {
  N <- length(tree$tip.label)
  gAfterInf <- g
  gAfterInf[tree$edge[, 2]] <- g[tree$edge[, 1]]
  gAfterInf[N+1] <- g[N+1]
  sapply(1:nSamples, function(i) {
    id <- sampleNodeIds(tree)
    edge <- tree$edge
    es <- matrix(e[id[edge]], ncol=2)
    if(length(grep('justAfter', sampTime[1]))>0) {
      gs <- matrix(gAfterInf[edge], ncol=2)
    } else if(length(grep('bothAt', sampTime[1]))>0) {
      gs <- matrix(g[edge[, 1]], nrow=nrow(edge), ncol=2)
    }
    pp <- gs+es
    edge <- matrix(id[edge], ncol=2)
    pp <- pp[edge[, 1]!=edge[, 2], ]
    cor(pp[,1], pp[,2])
  })
}


#' Infection couples given a tree and ids of all tips and nodes
#'
#' @param tree a phylo object
#' @param id a vector containing identifiers of all tips and nodes
#' @return a matrix of the class of id of two columns and N-1 rows
infectionCouples <- function(tree, id) {
  edge <- tree$edge
  edge <- matrix(id[edge], ncol=2)
  edge[edge[, 1]!=edge[, 2],]
}

#' One (out of many possible) identifications of the nodes and root
#'   with tips in the tree.
#' @param tree a phylo object with N tips
#' @return an integer vector of length the number of nodes (including tips) in the tree 
#' and elements in 1:N, the identities of all nodes in the tree (including the tips at
#' positions 1:N)
#' 
#' @description In a transmission tree, every node is identified as 
#'   one of its descendants, the person who transmitted the disease 
#'   to the other descendants of the node. The function provides one 
#'   (out of many possible) such identification, assuming that each 
#'   of the direct descendants of the node is equally likely to be 
#'   the node itself.
sampleNodeIds <- function(tree) {
  N <- length(tree$tip.label)
  edge.table <- as.data.table(tree$edge)
  setnames(edge.table, c('V1', 'V2'))
  setkey(edge.table, V1)
  src <- c(1:N, edge.table[, list(src=sample(V2, 1)), by=V1][[2]])
  M <- length(unique(as.vector(tree$edge)))  # number of all nodes 
  id <- rep(0,M)
  id2 <- 1:M
  while(any((id2-id)!=0)) {
    id <- id2
    id2 <- id2[src[id2]]
  }
  id
}

#'Group tips according to distance from the root
#' @param tree a phylo object
#' @param nGroups integer, the desired number of groups (default 15)
#' @param rootTipDists numeric vector of root to tip distances in the order of 
#'   tree$tip.label. If not passed, this vector is calculated by the function
#'   nodeTimes().
#' 
#' 
groupByRootDist <- function(tree, nGroups=15, rootTipDists=NULL) {
  N <- length(tree$tip.label)
  if(is.null(rootTipDists))
    rootTipDists <- nodeTimes(tree)[1:N]
  
  rootTipDistGroups <- cut(rootTipDists, breaks=nGroups, include.lowest = TRUE)
  names(rootTipDistGroups) <- tree$tip.label
  
  groupMeans <- sapply(levels(rootTipDistGroups), function(str) {
    mean(eval(parse(text=paste0('c(',substr(str, 2, nchar(str)-1), ')'))))
  })
  
  list(rootTipDists=rootTipDists, rootTipDistGroups=rootTipDistGroups, groupMeans=groupMeans)  
  
}
