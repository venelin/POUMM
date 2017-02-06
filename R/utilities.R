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
  
  nodesVector <- as.integer(c())
  nodesIndex <- as.integer(c(0))
  
  unVector <- as.integer(c())
  unIndex <- as.integer(c(0))
  
  # start by pruning the tips
  nodes <- as.integer(1:N)
  
  while(nodes[1] != N+1) {
    nodesIndex <- c(nodesIndex, nodesIndex[length(nodesIndex)] + length(nodes))
    nodesVector <- c(nodesVector, nodes)
    
    # indices in the edge-matrix of the edges pointing to the to be pruned nodes
    es <- endingAt[nodes]
    
    nodes <- as.integer(c())
    edgeEnds <- edge[es, 2]
    
    
    # start with a un-vector that doesn't contain indices in es, so es[-un] returns es intact.
    unAll <- length(es)+1
    lenUnAll <- 0L
    
    #update parent pifs
    while(lenUnAll != length(es)) {
    #while(length(es) > 0) {
      # indices in current es of the edges not yet added to their parent node.
      un <- match(unique(edge[es[-unAll], 1]), edge[es, 1])
      # un <- match(unique(edge[es, 1]), edge[es, 1])
      unAll <- c(unAll, un)
      lenUnAll <- lenUnAll + length(un)
      
      # attach to the vector of such indices
      unVector <- c(unVector, un)
      # attach the index of the current last element of unVector to unIndex
      unIndex <- c(unIndex, unIndex[length(unIndex)]+length(un))
      
      #pif[edge[es[un], 1], ] <- pif[edge[es[un], 1], ] + pif[edge[es[un], 2], ]
    
      # For the parent nodes, decrement the amount of non-visited children
      nonVisitedChildren[edge[es[un], 1]] <- nonVisitedChildren[edge[es[un], 1]] - 1
      # those parent nodes with no more non-visited children will be pruned next
      nodes <- c(nodes, edge[es[un][nonVisitedChildren[edge[es[un], 1]] == 0], 1])
      
      es[un] <- NA
    }
  }
  list(M=M, edge = edge, endingAt=endingAt, 
       nodesVector=nodesVector, nodesIndex=nodesIndex, 
       nLevels=length(nodesIndex)-1, 
       unVector=unVector, unIndex=unIndex)
}


printPOUMMLikelihoodMainLoop <- function(tree) {
  pruneInfo <- POUMM:::pruneTree(tree)
  M <- pruneInfo$M
  endingAt <- pruneInfo$endingAt
  nodesVector <- pruneInfo$nodesVector
  nodesIndex <- pruneInfo$nodesIndex
  nLevels <- pruneInfo$nLevels
  unVector <- pruneInfo$unVector
  unIndex <- pruneInfo$unIndex
  
  unJ <- 1
  N <- length(tree$tip.label)
  edge <- rbind(pruneInfo$edge, c(0,1))
  
  for(i in 1:(nLevels + 1)) {
    if(i <= nLevels) {
      es <- endingAt[nodesVector[(nodesIndex[i] + 1):nodesIndex[i + 1]]]  
    } else {
      es <- M
    }
    
    if(es[1] != M) {
      cat("Processing edges ending at nodes:\n")
      print(edge[es,2])
    } else {
      cat("Processing root node.\n")
    }
    
    if(es[1] != M) {
      lenUnAll <- 0L
      
      while(lenUnAll != length(es)) {
      # while(length(es)>0) {
        un <- unVector[(unIndex[unJ]+1):unIndex[unJ+1]]
        unJ <- unJ+1
        cat("Updating parent-nodes:\n")
        
        print(cbind(edge[es[un], 1], " <- ", edge[es[un], 2]))
        
        lenUnAll <- lenUnAll + length(un)
        
        #es <- es[-un]
      }  
    }
    
  }
}
plotPOUMMLikelihoodMainLoop <- function(tree, x.lim=c(-2,16.4)) {
  pruneInfo <- POUMM:::pruneTree(tree)
  M <- pruneInfo$M
  endingAt <- pruneInfo$endingAt
  nodesVector <- pruneInfo$nodesVector
  nodesIndex <- pruneInfo$nodesIndex
  nLevels <- pruneInfo$nLevels
  unVector <- pruneInfo$unVector
  unIndex <- pruneInfo$unIndex
  
  unJ <- 1
  N <- length(tree$tip.label)
  
  edge <- rbind(pruneInfo$edge, c(0,1))
  
  edgeColor <- rep("black", nrow(edge))
  tipBg <- "white"
  tipColor <- "red"
  tipFrame <- "circle"
  
  nodeBg <- rep("black", length(tree$node.label))
  nodeCol <- rep("white", length(tree$node.label))
  
  names(nodeBg) <- names(nodeCol) <- tree$node.label
  
  for(i in 1:(nLevels + 1)) {
    if(i <= nLevels) {
      es <- endingAt[nodesVector[(nodesIndex[i] + 1):nodesIndex[i + 1]]]  
    } else {
      es <- M
    }
    
    edgeColor[es] <- "red"
    
    plot(tree, type="cladogram", show.tip.label=FALSE, edge.color = edgeColor, edge.width=2, x.lim=x.lim)
    
    if(tipBg == "darkgrey") {
      parColDefault <- par(col="darkgrey")
    } else {
      parColDefault <- par(col="red")
    }
    tiplabels(text = tree$tip.label, 
              frame = tipFrame, 
              bg = tipBg, 
              col = tipColor, cex = 1.4)
    par(col=parColDefault)
    
    # update these for the next iteration
    tipBg <- "darkgrey"
    tipColor <- "darkgrey"
    
    if(i > 1) {
      # not at tip-level
      nodeBg[as.character(edge[es, 2] - 1)] <- "white"
      nodeCol[as.character(edge[es, 2] - 1)] <- "red"
    }
    
    if(any(nodeCol == "red")) {
      parColDefault <- par(col="red")
      nodelabels(text = tree$node.label[nodeCol=="red"],
                 node = which(nodeCol == "red") + N,
                 frame = "circle",
                 bg = nodeBg[nodeCol == "red"],
                 col=nodeCol[nodeCol == "red"],
                 cex = 1)
      par(col=parColDefault)
    }
    if(any(nodeCol == "darkgrey")) {
      parColDefault <- par(col="darkgrey")
      nodelabels(text = tree$node.label[nodeCol=="darkgrey"],
                 node = which(nodeCol == "darkgrey") + N,
                 frame = "circle",
                 bg = nodeBg[nodeCol == "darkgrey"],
                 col = nodeCol[nodeCol == "darkgrey"],
                 cex = 1)
      par(col=parColDefault)
    }
    if(any(nodeCol == "white")) {
      parColDefault <- par(col="black")
      nodelabels(text = tree$node.label[nodeCol=="white"],
                 node = which(nodeCol == "white") + N,
                 frame = "circle",
                 bg = nodeBg[nodeCol == "white"],
                 col = nodeCol[nodeCol == "white"],
                 cex = 1)
      par(col=parColDefault)
    }
    
    if(i > 1) {
      # not at tip-level
      nodeBg[as.character(edge[es, 2] - 1)] <- "darkgrey"
      nodeCol[as.character(edge[es, 2] - 1)] <- "darkgrey"
    }
    
    edgeColor[es] <- "darkgrey"
    
    if(es[1] != M) {
      cat("Processing edges ending at nodes:\n")
      print(edge[es,2])
    } else {
      cat("Processing root node.\n")
    }
    
    if(es[1] != M) {
      queueNo <- 1
      letters <- c('a','b','c')
      
      #lenUnAll <- 0L
      
      #while(lenUnAll != length(es)) {
      while(length(es)>0) {
        un <- unVector[(unIndex[unJ]+1):unIndex[unJ+1]]
        unJ <- unJ+1
        cat("Updating parent-nodes:\n")
        
        print(cbind(edge[es[un], 1], " <- ", edge[es[un], 2]))
        parColDefault <- par(col="yellow")
        edgelabels(paste0(if(i!=3) "   " else "", letters[queueNo], if(i < nLevels) " " else ""), es[un], cex = 1, 
                   frame="none", col="red", 
                   adj=c(0.6,1.4))
        par(col=parColDefault)
        
        #lenUnAll <- lenUnAll + length(un)
        
        es <- es[-un]
        queueNo <- queueNo + 1
      }  
    }
    
  }
}

#' Writes verbose messages of the order of tree traversal during likelihood calculation
#' @param tree A phylo object.
#' @return Nothing
#' @export
simulatePOUMMLikelihoodMainLoop <- function(tree) {
  pruneInfo <- pruneTree(tree)
  M <- pruneInfo$M
  endingAt <- pruneInfo$endingAt
  nodesVector <- pruneInfo$nodesVector
  nodesIndex <- pruneInfo$nodesIndex
  nLevels <- pruneInfo$nLevels
  unVector <- pruneInfo$unVector
  unIndex <- pruneInfo$unIndex
  edge <- pruneInfo$edge
  unJ <- 1
  
  
  for(i in 1:nLevels) {
    es <- endingAt[nodesVector[(nodesIndex[i] + 1):nodesIndex[i + 1]]]
    
    cat("Processing edges ending at nodes:\n")
    print(edge[es,2])
    
    while(length(es)>0) {
      un <- unVector[(unIndex[unJ]+1):unIndex[unJ+1]]
      unJ <- unJ+1
      cat("Updating parent-nodes:\n")
      print(cbind(edge[es[un], 1], " <- ", edge[es[un], 2]))
    
      es <- es[-un]
    }
  }
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
