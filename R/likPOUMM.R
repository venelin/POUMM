#' Random generation of values along a phylogenetic tree following a branching
#' OU process
#' 
#' @param z0 Numeric value indicating the initial state at the root of the tree.
#' @param alpha,theta,sigma Numeric values, parameters of the OU process.
#' @param sigmae Numeric non-negative value (default 0). Specifies the standard
#'   deviation of random environmental contribution and or measurement standard 
#'   error to be added to the values (white noise). Note that if measurement standard
#'   error, se, is known and needs to be added to the environmental contribution,
#'   the right way to specify the parameter would be sqrt(sigmae^2+se^2).
#' @param tree An object of class phylo (package ape).
#' 
#' @return A numeric vector containing the generated values at all nodes
#'   (internal and tips) of the phylogenetic tree. 
#' @examples   
#' \dontrun{
#'   N <- 20 
#'   tree <- TreeSim::sim.bd.taxa(N, 1, lambda = 2, mu = 1, complete =
#'   FALSE)[[1]] 
#'   z <- rVNodesGivenTreePOUMM(tree, 8, 2, 2, 1)
#'   
#'   phytools::plotBranchbyTrait(tree, z, mode = "nodes", show.tip.label =
#'   FALSE, show.node.label = FALSE) 
#'   ape::nodelabels(round(z[-(1:N)], 2)) 
#'   ape::tiplabels(round(z[(1:N)], 2))
#' } 
#' @export
rVNodesGivenTreePOUMM <- function(tree, z0, alpha, theta, sigma, sigmae = 0) {
  g <- ape::rTraitCont(tree, 
                       function(x, l, .a, .t, .s) {
                         POUMM::rOU(n = 1, z0 = x, t = l, alpha = .a, theta = .t, sigma = .s)
                       }, root.value = z0, ancestor = TRUE, 
                       .a = alpha, .t = theta, .s = sigma)
  if(any(sigmae > 0)) {
    if(length(sigmae) < length(g))
      sigmae <- rep(sigmae[1], length(g))
    e <- rnorm(length(g), 0, sigmae) 
    g <- g+e
  } 
  
  g
}

#' Multivariate density of observed values along a tree given an OU process of 
#' evolution and root-value
#' 
#' @description Calculates the conditional probability density of observed 
#'   values at the tips and internal nodes of a tree, given that tree, the value
#'   at the root, z[N+1], where N is the number of tips in the tree, known 
#'   measurement error e for each value in z, and a POUMM model of evolution. 
#'   This function is mostly used to calculate the likelihood of simulated data
#'   under known model parameters.
#'   
#' @param z A numeric vector of size length(tree$tip.label)+tree$Nnode 
#'   representing the observed values at the tips, root and internal nodes.
#' @param tree An object of class phylo.
#' @param alpha,theta,sigma Numeric values, parameters of the OU model.
#' @param sigmae Numeric non-negative value or vector of length(z) elements 
#' (default 0). Specifies the
#'  standard deviation of random environmental contribution (and eventually 
#'  measurement error) to be added to the values. Note that if measurement standard
#'  error, se, is known and needs to be added to the environmental contribution,
#'  the right way to specify the parameter would be sqrt(sigmae^2+se^2), not
#'  sigmae+se.
#' @param e Numeric vector of size length(z) representing exactly known 
#'   error (sum of environmental contribution and measurement error). 
#'   Defaults to a vector of zeroes.
#' @param log Logical indicating whether a log-likelihood should be returned 
#'   instead of a likelihood. Default is TRUE.
#'   
#' @return A numeric value, the multivariate probability density of z under the 
#'   given parameters.
#' @export
dVNodesGivenTreePOUMM <- function(z, tree, alpha, theta, sigma, sigmae=0,
                                  e = rep(0, length(z)),  
                                  log=TRUE) {
  if (is.null(tree$edge.length)) {
    stop("tree has no branch length")
  }
  if (any(tree$edge.length < 0)) {
    stop("at least one branch length negative")
  }
  
  g <- z - e
  
  M <- dim(tree$edge)[1]            # number of edges in the tree
  anc <- tree$edge[, 1]
  des <- tree$edge[, 2]
  el <- tree$edge.length
  probs <- sapply(1:M, function(i) dOU(g[des[i]], g[anc[i]], el[i], 
                                       alpha, theta, sigma, log=TRUE))
  
  if(sigmae > 0) {
    de <- dnorm(e, 0, sigmae, log = TRUE)
  } else {
    de <- 0
  }
  
  ifelse(log, sum(probs, de), exp(sum(probs, de)))
}

# loading the IntegratorPOUMM C++ module
loadModule( "IntegratorPOUMM", TRUE )

#' Extract information for fast likelihood calculation using the breadth-first
#' pruning algorithm. 
#' 
#' @param tree a phylo object
#' @param z Numeric vector with length(tree$tip.label) values corresponding to
#' tree$tip.label.
#' @param se Non-negative numerical or N-vector indicating known standard 
#' measurement error.
#' 
#' @import Rcpp
#' @return a list of objects used for likelihood evaluation
#' 
#' @export
pruneTree <- function(tree, z, se = 0) {
  
  N <- length(tree$tip.label)                # number of tips
  M <- length(unique(as.vector(tree$edge)))  # number of all nodes 
  endingAt <- order(rbind(tree$edge, c(0, N+1))[, 2])
  
  edge <- tree$edge
  mode(tree$edge) <- "integer"
  
  if(length(se) != N) {
    se <- rep(se[1], N)
  }
  
  if(any(is.na(se) | se < 0)) {
    stop("All elements of se should be non-negative numbers.")
  }
  
  nonVisitedChildren <- rep(0, M)
  ee1 <- tree$edge[, 1]
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
    edgeEnds <- tree$edge[es, 2]
    
    
    # start with a un-vector that doesn't contain indices in es, so es[-un] returns es intact.
    unAll <- length(es)+1
    lenUnAll <- 0L
    
    #update parent pifs
    while(lenUnAll != length(es)) {
      #while(length(es) > 0) {
      # indices in current es of the edges not yet added to their parent node.
      un <- match(unique(tree$edge[es[-unAll], 1]), tree$edge[es, 1])
      # un <- match(unique(tree$edge[es, 1]), tree$edge[es, 1])
      unAll <- c(unAll, un)
      lenUnAll <- lenUnAll + length(un)
      
      # attach to the vector of such indices
      unVector <- c(unVector, un)
      # attach the index of the current last element of unVector to unIndex
      unIndex <- c(unIndex, unIndex[length(unIndex)]+length(un))
      
      #pif[tree$edge[es[un], 1], ] <- pif[tree$edge[es[un], 1], ] + pif[tree$edge[es[un], 2], ]
      
      # For the parent nodes, decrement the amount of non-visited children
      nonVisitedChildren[tree$edge[es[un], 1]] <- nonVisitedChildren[tree$edge[es[un], 1]] - 1
      # those parent nodes with no more non-visited children will be pruned next
      nodes <- c(nodes, tree$edge[es[un][nonVisitedChildren[tree$edge[es[un], 1]] == 0], 1])
      
      es[un] <- NA
    }
  }
  
  # create an Integrator object and initialize it with the pruning informaiton
  
  integrator <- Integrator$new()
  integrator$setPruningInfo(
    z, se, tree$edge, tree$edge.length,
    M, N, endingAt, nodesVector, nodesIndex, unVector, unIndex)
  
  list(M=M, N=N, 
       z = z, se = se, tree = tree, 
       endingAt=endingAt, 
       nodesVector=nodesVector, nodesIndex=nodesIndex, 
       nLevels=length(nodesIndex)-1, 
       unVector=unVector, unIndex=unIndex,
       integrator = integrator)
}

#' @name likPOUMM_C
#' 
#' @title Fast POUMM likelihood calculation based on the breadth-first 
#' pruning algorithm.
#'  
#' @description The function likPOUMMGivenTreeVTipsC is using elementwise 
#'   vector operations from the Armadillo library. The function
#'   likPOUMMGivenTreeVTipsC2 is using C++ for loops which benefit from
#'   better vectorization of native algebraic expressions. The function 
#'   likPOUMMGivenTreeVTipsC3 adds to likPOUMMGivenTreeVTipsC2 Open MP 
#'   parallelization (mostly available on multi-core Linux systems). 
#'
#' @param integrator An object of the Integrator C++ class (see example).
#'   trait values at the tip-nodes.
#' @param alpha the strength of the selection
#' @param theta long term mean value of the OU process
#' @param sigma the unit-time standard deviation of the random component in the
#'   OU process.
#' @param sigmae the standard deviation of the environmental deviation added to
#'   the genetic contribution at each tip, by default 0, meaning no
#'   environmental deviation.
#' @param g0 Numeric, NA or NaN, either a fixed genotypic value at the root of 
#'   tree or NA or NaN. A NA "Not Available" will cause to analytically
#'   calculate the value of g0 that would maximize the conditional likelihood of
#'   the data given g0. A NaN "Not a Number" will cause integration over g0
#'   taking values in (-Inf,+Inf) assuming that g0 is normally distributed with
#'   mean g0Prior$mean and variance g0Prior$var (see parameter g0Prior).
#' @param g0Prior Either NULL or a list with named numeric or character members
#'  "mean" and "var". Specifies a prior normal distribution for the parameter g0.
#'  If characters, the members mean and var are evaluated as R-expressions.
#' @param log Logical indicating whether log-likelihood should be returned
#'   instead of likelihood, default is TRUE.
#'
#' @details This function is the C++ equivalent of dVTipsGivenTreePOUMM (aliased also as
#' likPOUMMGivenTreeVTips). C++ implementation using contiguous vector operations and
#'  vectors a, b and c
#' @return A numeric with attributes "g0" and "g0LogPrior".
#' 
#' @examples 
#' \dontrun{
#' N <- 100
#' tr <- ape::rtree(N)
#' z <- POUMM::rVNodesGivenTreePOUMM(tr, 0, 2, 3, 1, 1)[1:N]
#' pruneInfo <- POUMM::pruneTree(tr, z)
#' microbenchmark::microbenchmark(
#'   likArma <- POUMM::likPOUMMGivenTreeVTipsC(pruneInfo$integrator, 2, 3, 1, 1),
#'   likR <- POUMM::likPOUMMGivenTreeVTips(z, tr, 2, 3, 1, 1, pruneInfo = pruneInfo))
#' 
#' # should be the same values
#' likArma
#' likR
#' }
#' 
NULL

#' @describeIn likPOUMM_C using RcppArmadillo and simd vectorization
#' @export
likPOUMMGivenTreeVTipsC <- function(
  integrator, alpha, theta, sigma, sigmae, g0 = NA, g0Prior = NULL, log = TRUE){
  
  if(anyNA(c(alpha, theta, sigma, sigmae))) { # case 9 and 10
    warning('Some parameters are NA or NaN')
    ifelse(log, -Inf, 0)
  } else if((alpha > 0 & sigma == 0 & sigmae == 0) |  # case 4
            all(c(alpha, sigma, sigmae) == 0) ) {        # case 8
    ifelse(log, -Inf, 0)
  } else {
    abc <- integrator$abc_arma(alpha, theta, sigma, sigmae)
    if(any(is.nan(abc)) | abc[1]==0 | is.infinite(abc[3])) {
      value <- NaN
      attr(value, "g0") <- g0
      attr(value, "g0LogPrior") <- NA
    } else {
      resList <- loglik_abc_g0_g0Prior(abc, alpha, theta, sigma, g0, g0Prior)
      loglik <- resList$loglik
      g0 <- resList$g0
      g0LogPrior <- resList$g0LogPrior
      
      value <- if(log) loglik else exp(loglik)
      attr(value, "g0") <- g0
      attr(value, "g0LogPrior") <- g0LogPrior
    }
    
    value
  }
}

#' @describeIn likPOUMM_C Using omp simd operations in C++ for-loops
#' @export
likPOUMMGivenTreeVTipsC2 <- function(
  integrator, alpha, theta, sigma, sigmae, g0 = NA, g0Prior = NULL, log = TRUE){
  
  if(anyNA(c(alpha, theta, sigma, sigmae))) { # case 9 and 10
    warning('Some parameters are NA or NaN')
    ifelse(log, -Inf, 0)
  } else if((alpha > 0 & sigma == 0 & sigmae == 0) |  # case 4
            all(c(alpha, sigma, sigmae) == 0) ) {        # case 8
    ifelse(log, -Inf, 0)
  } else {
    abc <- integrator$abc_omp_simd(alpha, theta, sigma, sigmae)
    if(any(is.nan(abc)) | abc[1]==0 | is.infinite(abc[3])) {
      value <- NaN
      attr(value, "g0") <- g0
      attr(value, "g0LogPrior") <- NA
    } else {
      resList <- loglik_abc_g0_g0Prior(abc, alpha, theta, sigma, g0, g0Prior)
      loglik <- resList$loglik
      g0 <- resList$g0
      g0LogPrior <- resList$g0LogPrior
      
      value <- if(log) loglik else exp(loglik)
      attr(value, "g0") <- g0
      attr(value, "g0LogPrior") <- g0LogPrior
    }
    
    value
  }
}

#' @describeIn likPOUMM_C Using OpenMP and SIMD parallelization if available on 
#' the system. This implementation differs from likPOUMMGivenTreeVTipsC3 by
#' the use of longer code-blocks in for-loops (better adapted for the Intel C++
#' compiler). 
#' @export
likPOUMMGivenTreeVTipsC4 <- function(
  integrator, alpha, theta, sigma, sigmae, g0 = NA, g0Prior = NULL, log = TRUE){
  
  if(anyNA(c(alpha, theta, sigma, sigmae))) { # case 9 and 10
    warning('Some parameters are NA or NaN')
    ifelse(log, -Inf, 0)
  } else if((alpha > 0 & sigma == 0 & sigmae == 0) |  # case 4
            all(c(alpha, sigma, sigmae) == 0) ) {        # case 8
    ifelse(log, -Inf, 0)
  } else {
    abc <- integrator$abc_omp_for_simd(alpha, theta, sigma, sigmae)
    if(any(is.nan(abc)) | abc[1]==0 | is.infinite(abc[3])) {
      value <- NaN
      attr(value, "g0") <- g0
      attr(value, "g0LogPrior") <- NA
    } else {
      resList <- loglik_abc_g0_g0Prior(abc, alpha, theta, sigma, g0, g0Prior)
      loglik <- resList$loglik
      g0 <- resList$g0
      g0LogPrior <- resList$g0LogPrior
      
      value <- if(log) loglik else exp(loglik)
      attr(value, "g0") <- g0
      attr(value, "g0LogPrior") <- g0LogPrior
    }
    
    value
  }
}


#' Processing of the root value and calculation of the maximum log-likelihood 
#' for the given coefficients abc, and parameters theta, g0 and g0Prior. This is
#' an internal function.
#' 
#' @param abc a vector of 3 numerical values denoting the corresponding 
#'   coefficients in the POUMM likelihood presented as exp(a g0^2 + b g0 + c).
#' @param alpha,theta,sigma parameters of the OU-process.
#' @param g0 initial value at the root of the tree (can be NA). See argument
#'   parMapping in ?specifyPOUMM.
#' @param g0Prior list object. See parameter g0Prior in ?specifyPOUMM.
#'   
loglik_abc_g0_g0Prior <- function(abc, alpha, theta, sigma, g0, g0Prior) {
  logpi <- log(pi)
  loge2 <- log(2)
  
  if(is.list(g0Prior)) {
    if(is.null(g0Prior$mean) | is.null(g0Prior$var)) {
      stop("g0Prior must have a 'mean' and a 'var' character members.")
    } else {
      g0Mean <- eval(parse(text=g0Prior$mean))
      g0Var <- eval(parse(text=g0Prior$var))
    }
  } else {
    g0Mean <- NA
    g0Var <- NA
  }
  
  # Processing of the root value
  if(is.finite(g0)) {
    # fixed g0 specified as parameter; 
    
    # pdf of the data conditioned on fixed g0
    g0_theta <- g0 - theta
    loglik <- (abc[1] * g0_theta + abc[2]) * g0_theta + abc[3]
    
    if(is.list(g0Prior)) {
      # Normal prior for g0 with mean g0Mean and variance g0Var
      g0LogPrior <- dnorm(as.double(g0), 
                          mean = as.double(g0Mean), 
                          sd = sqrt(as.double(g0Var)), log = log)
    } else {
      # No prior for g0 specified
      g0LogPrior <- NA
    }
  } else if(is.na(g0) & !is.nan(g0)) {
    # find g0 that would maximise one of the following:
    # (1) if a normal prior for g0 was specified, 
    #     pdf(z | alpha, theta, sigma, sigmae, g0, tree) x prior(g0)
    # (2) otherwise, pdf(z | alpha, theta, sigma, sigmae, g0, tree) 
    if(is.list(g0Prior)) {
      # (1) Normal prior for g0 with mean g0Mean and variance g0Var
      # max(loglik + g0Prior): derivative : 
      # 2 * (abc[1] - 1 / (2 * g0Var)) * g0 + (abc[2] + g0Mean / g0Var) == 0
      g0_theta <- 
        -0.5 * (abc[2] + (g0Mean - theta) / g0Var) / (abc[1] - 1 / (2 * g0Var))
      
      loglik <- (abc[1] * g0_theta + abc[2]) * g0_theta + abc[3]
      
      g0 <- g0_theta + theta
      
      g0LogPrior <- dnorm(as.double(g0), 
                          mean = as.double(g0Mean), 
                          sd = as.double(sqrt(g0Var)), log = log)
    } else {
      # (2) No prior for g0 specified
      # max(loglik): derivative : 2 * abc[1] * g0 + abc[2] == 0
      if(abc[1] != 0) {
        g0_theta <- -0.5 * abc[2] / abc[1]
      } else {
        g0_theta <- 0.5 * abc[2] / .Machine$double.xmin
      }
      
      
      loglik <- (abc[1] * g0_theta + abc[2]) * g0_theta + abc[3]
      
      g0 <- g0_theta + theta
      
      g0LogPrior <- NA
    }
  } else {
    # g0 is NaN or infinite. In this case, 
    # we treat g0 as not specified, so the only option is to integrate under 
    # its prior distribution. If g0Prior is not specified, assume that g0Prior 
    # is the stationary OU normal distribution with mean, theta, and variance, 
    # varOU(Inf, alpha, sigma)
    if(!is.list(g0Prior)) {
      g0Var <- varOU(Inf, alpha, sigma)
    }
    
    if(!is.na(g0Var) & g0Var != 0) {
      def <- c(-0.5 / g0Var, 0, - log(sqrt(g0Var)) - 0.5 * (loge2+logpi))
      abc <- abc + def
      loglik <- abc[3] + 0.5 * (logpi - log(-abc[1])) - abc[2]^2 / (4 * abc[1])
      g0 <- NaN
      g0LogPrior <- NaN
    } else {
      #g0Var = 0 or g0Var is NA or NaN: assume gr=g0Mean with probability 1
      loglik <- (abc[1] * (g0Mean[1]-theta) + abc[2]) * (g0Mean[1]-theta) + abc[3]
      
      g0 <- g0Mean
      
      g0LogPrior <- Inf
    }
  }
  list(loglik = loglik, g0 = g0, g0LogPrior = g0LogPrior)
}



#' Density of observed tip-values given a tree, assuming Ornstein-Uhlenbeck 
#' process for the genetic contributions along the tree and normally distributed
#' environmental deviations.
#' 
#' @description Calculates the (log-)probability density of trait values at the
#' tip-nodes given the tree and assuming that the trait value at the tips is the
#' sum of a genetic contribution, g, that evolved on the tree according to an OU
#' process with parameters alpha, theta, sigma and an environmental deviation, 
#' e, that is distributed normally and independently between the tips of the 
#' tree. Note: Without additional assumptions for the distribution of the value
#' at the root of the tree, the likelihood is not defined at alpha=0, although 
#' this corresponds to the limiting Brownian motion process with mean value 
#' theta and unit time variance sigma^2. Considering the observed data and tree as
#' a fixed parameter and the POUMM parameters as variables, this function is
#' interpreted as the POUMM likelihood.
#' 
#' @param tree an object of class phylo
#' @param z A numeric vector of size length(tree$tip.label) representing the
#'   trait values at the tip-nodes.
#' @param alpha the strength of the selection
#' @param theta long term mean value of the OU process
#' @param sigma the unit-time standard deviation of the random component in the
#'   OU process.
#' @param sigmae the standard deviation of the environmental deviation added to
#'   the genetic contribution at each tip, by default 0, meaning no
#'   environmental deviation.
#' @param g0 Numeric, NA or NaN, either a fixed genotypic value at the root of 
#'   tree or NA or NaN. A NA "Not Available" will cause to analytically
#'   calculate the value of g0 that would maximize the conditional likelihood of
#'   the data given g0. A NaN "Not a Number" will cause integration over g0
#'   taking values in (-Inf,+Inf) assuming that g0 is normally distributed with
#'   mean g0Prior$mean and variance g0Prior$var (see parameter g0Prior).
#' @param g0Prior Either NULL or a list with named numeric or character members
#'  "mean" and "var". Specifies a prior normal distribution for the parameter g0.
#'  If characters, the members mean and var are evaluated as R-expressions.
#' @param log Logical indicating whether log-likelihood should be returned
#'   instead of likelihood, default is TRUE.
#' @param pruneInfo list returned by pruneTree(tree) to be passed in explicit
#'   calls to dVGivenTreeOU.
#' @param usempfr integer indicating if and how mpfr should be used for small 
#'  parameter values (any(c(alpha, sigma, sigmae) < 0.01)). Using the mpfr 
#'  package can be forced by specifying an integer greater or equal to 2. 
#'  Setting usempfr=0 (default) causes high precision likelihood 
#'  calculation to be done on each encounter of parameters with at least 1 bigger
#'  log-likelihood value than any of the currently found
#'  maximum log-likelihood or the previously calculated log-likelihood value
#'  Requires the Rmpfr package. Note that using mpfr may increase the time for 
#'  one likelihood calculation more than 100-fold. Set usempfr to -1 or less
#'  to completely disable Rmpfr functionality. 
#' @param maxmpfr integer (not used)
#' @param precbits integer specifying precision bits for mpfr. Default is 512.
#' @param debug logical, if set to TRUE some debugging information is printed 
#'   during likelihood calculation
#' 
#' @return A numeric with attributes "g0" and "g0LogPrior".
#' 
#' @import data.table
#' @aliases likPOUMMGivenTreeVTips dVTipsGivenTreePOUMM
#' @export
likPOUMMGivenTreeVTips <- dVTipsGivenTreePOUMM <- function(
  z, tree, alpha, theta, sigma, sigmae = 0, g0 = NA, g0Prior = NULL, 
  log = TRUE, pruneInfo = pruneTree(tree, z), 
  usempfr = 0, maxmpfr = 2, precbits = 128, debug = FALSE) {
  
  availRmpfr <- requireNamespace("Rmpfr", quietly = TRUE) 
  
  alphaorig <- alpha
  thetaorig <- theta
  sigmaorig <- sigma
  sigmaeorig <- sigmae
  g0orig <- g0
  vorig <- z
  
  if(anyNA(c(alpha, theta, sigma, sigmae))) { # case 9 and 10
    warning('Some parameters are NA or NaN')
    ifelse(log, -Inf, 0)
  } else if((alpha > 0 & sigma == 0 & all(sigmae == 0)) |  # case 4
     all(c(alpha, sigma, sigmae) == 0) ) {        # case 8
    ifelse(log, -Inf, 0)
  } else {
    logpi <- log(pi)
    loge2 <- log(2)
    
    N <- length(tree$tip.label)                # number of tips
    
    if(is.null(pruneInfo)) {
      pruneInfo <- pruneTree(tree, z)
    }
    
    edge <- pruneInfo$tree$edge
    
    M <- pruneInfo$M
    endingAt <- pruneInfo$endingAt
    nodesVector <- pruneInfo$nodesVector
    nodesIndex <- pruneInfo$nodesIndex
    nLevels <- pruneInfo$nLevels
    unVector <- pruneInfo$unVector
    unIndex <- pruneInfo$unIndex
    
    t <- pruneInfo$tree$edge.length
    
    done <- FALSE
    
    # this loop governs the use of mpfr; mpfr is not used when backend = "C"
    while(!done & usempfr <= maxmpfr) { 
      if(availRmpfr & usempfr >= 2) {
        if(is.double(alphaorig)) {
          alpha <- Rmpfr::mpfr(alphaorig, precbits*2^(usempfr-2))
        }
        if(is.double(thetaorig)) {
          theta <- Rmpfr::mpfr(thetaorig, precbits*2^(usempfr-2))
        }
        if(is.double(sigmaorig)) {
          sigma <- Rmpfr::mpfr(sigmaorig, precbits*2^(usempfr-2))
        }
        if(is.double(sigmaeorig)) {
          sigmae <- Rmpfr::mpfr(sigmaeorig, precbits*2^(usempfr-2))
        }
        if(is.finite(g0)) {
          g0 <- Rmpfr::mpfr(g0orig, precbits*2^(usempfr-2))
        }
        if(is.double(z)) {
          z <- Rmpfr::mpfr(z, precbits*2^(usempfr-2))
        }
      }
      
      
      # use is.double(alpha) to check whether Rmpfr is being used
      
      sigma2 <- sigma^2
      sigmae2 <- sigmae^2
      logsigma <- log(sigma)
      logsigmae <- log(sigmae)
      
      talpha <- t*alpha
      etalpha <- exp(talpha)
      
      e2talpha <- etalpha*etalpha
      
      limitfe2talpha <- -0.5 / t
      if(alpha != 0) {
        fe2talpha <- alpha / (1 - e2talpha)
      } else {
        fe2talpha <- limitfe2talpha
      }
      
      if(availRmpfr & usempfr & usempfr < maxmpfr) {
        if(any(fe2talpha >= 0) | any(fe2talpha < limitfe2talpha)) {
          usempfr <- usempfr + 1
          next 
        }
      }
      
      pif <- rep(alpha*0, M*3)
      dim(pif) <- c(M, 3)
      
      gutalphasigma2 <- rep(alpha*0, M)
      
      unJ <- 1
      
      for(i in 1:nLevels) {
        es <- endingAt[nodesVector[(nodesIndex[i] + 1):nodesIndex[i + 1]]]
        edgeEnds <- edge[es, 2]
        
        if(edge[es[1], 2] <= N) {
          gutalphasigma2[es] <- e2talpha[es] + ((-0.5 / sigmae2) * sigma2) / fe2talpha[es]
          # all es pointing to tips
          if(all(sigmae==0)) {
            # no environmental deviation
            # (possibly reordered) shifted tip values
            g1 <- z[edgeEnds] - theta 
            
            pif[edgeEnds, 1] <- fe2talpha[es] / sigma2
            pif[edgeEnds, 2] <- -2 * etalpha[es] * g1 * pif[edgeEnds, 1]
            pif[edgeEnds, 3] <- talpha[es] + 0.5*log(-fe2talpha[es]) - 
              0.5*logpi - logsigma + (etalpha[es] * g1)^2 * pif[edgeEnds, 1]
          } else {
            # (possibly reordered) shifted tip values
            z1 <- z[edgeEnds] - theta 
            
            pif[edgeEnds, 3] <- -0.5 * log(gutalphasigma2[es]) -
              0.25 * sigma2 * (z1 / sigmae2)^2 / 
              (fe2talpha[es] - alpha + (-0.5 / sigmae2) * sigma2) +
              talpha[es] + (-0.5 * (loge2 + logpi  + z1 * (z1 / sigmae2)) - logsigmae)
            
            pif[edgeEnds, 2] <- (etalpha[es] * (z1 / sigmae2)) / gutalphasigma2[es]
            
            pif[edgeEnds, 1] <- (-0.5 / sigmae2) / gutalphasigma2[es]
          }
        } else {
          # edges pointing to internal nodes, for which all children nodes have been visited
          gutalphasigma2[es] <- e2talpha[es] + (pif[edgeEnds, 1] * sigma2) / fe2talpha[es]
          
          pif[edgeEnds, 3] <- -0.5 * log(gutalphasigma2[es]) -
            # use the formula log(a+c) = log(a) + log(1+c/a), where 
            # a=e2talpha, c = (pif[edgeEnds, 1] * sigma2) / fe2talpha[es], b = e.
            # 
            
            #pif[edgeEnds, 3] <- 
            #  -0.5 * (2 * talpha[es] + 
            #            log1p((pif[edgeEnds, 1] * sigma2) / fe2talpha[es] / e2talpha[es])) -
            0.25 * sigma2 * pif[edgeEnds, 2]^2 / 
            (fe2talpha[es] - alpha + pif[edgeEnds, 1] * sigma2) +
            talpha[es] + pif[edgeEnds, 3]
          
          pif[edgeEnds, 2] <- (etalpha[es] * pif[edgeEnds, 2]) / gutalphasigma2[es]
          
          pif[edgeEnds, 1] <- pif[edgeEnds, 1] / gutalphasigma2[es]
        } 
        
        #update parent pifs
        lenUnAll <- 0L
        
        while(lenUnAll != length(es)) {
          #while(length(es) > 0) {
          # cat('unIndex: length:', length(unIndex), ", value:", toString(unIndex), '\n')
          # cat('unJ:', unJ, "\n")
          # cat('lenUnAll:', lenUnAll, ', length(es):', length(es), "\n")
          # cat("=======")

          un <- unVector[(unIndex[unJ]+1):unIndex[unJ+1]]
          unJ <- unJ+1
          lenUnAll <- lenUnAll + length(un)
          
          pif[edge[es[un], 1], ] <- 
            pif[edge[es[un], 1], ] + pif[edge[es[un], 2], ]
          
          #es <- es[-un]
        }
      }
      #print(pif)
      # at the root: P[Vtips|tree, g_root, alpha, theta, sigma, sigmae] =
      #                  abc[1]*g_root^2 + abc[2]*g_root + abc[3]
      abc <- pif[N+1, ]
      
      if(availRmpfr & usempfr & (any(is.nan(abc)) | abc[1]==0 | is.infinite(abc[3]))) {
        usempfr <- usempfr + 1
        next 
      }
      done <- TRUE
    }
    
    resList <- loglik_abc_g0_g0Prior(abc, alpha, theta, sigma, g0, g0Prior)
    
    loglik <- resList$loglik
    g0 <- resList$g0
    g0LogPrior <- resList$g0LogPrior
    
    if(is.double(alphaorig) & is.double(thetaorig) & is.double(sigmaorig) & 
       is.double(sigmaeorig) & is.double(vorig)) {
      loglik <- as.double(loglik)
      g0 <- as.double(g0)
      g0LogPrior <- as.double(g0LogPrior)
    }
    
    value <- if(log) loglik else exp(loglik)
    attr(value, "g0") <- g0
    attr(value, "g0LogPrior") <- g0LogPrior
    
    if(debug) {
      debugdata <- data.table(
        alpha = list(alpha), theta = list(theta), sigma = list(sigma), 
        sigmae = list(sigmae),
        sigma2 = list(sigma2),
        sigmae2 = list(sigmae2), logsigma = list(logsigma), 
        logsigmae = list(logsigmae), 
        t = list(t),
        talpha = list(talpha), etalpha = list(etalpha), 
        e2talpha = list(e2talpha), 
        fe2talpha = list(fe2talpha),
        pif = list(pif), abc = list(abc), 
        g0 = list(g0), 
        g0LogPrior = list(g0LogPrior), loglik = list(loglik), 
        availRmpfr = list(availRmpfr), usempfr = list(usempfr), 
        precbits = list(precbits))
      attr(value, "debugdata") <- debugdata
    }
    
    value
  }
}

#' old likelihood version including theta in expressions
#' to be adapted for the case of multiple regimes.
# dVTipsGivenTreePOUMMwithThetaInExpressions <- function(
#   z, tree, alpha, theta, sigma, sigmae = 0, g0 = NA, g0Prior = NULL,
#   log = TRUE, pruneInfo = pruneTree(tree), 
#   usempfr = 0, maxmpfr = 2, precbits = 128, 
#   debug = FALSE) {
#   
#   availRmpfr <- requireNamespace("Rmpfr", quietly = TRUE) 
#   
#   alphaorig <- alpha
#   thetaorig <- theta
#   sigmaorig <- sigma
#   sigmaeorig <- sigmae
#   g0orig <- g0
#   vorig <- z
#   
#   if(any(tree$edge.length <= 0)) {
#     stop('All branch lengths in tree should be positive!')
#   }
#   if(anyNA(c(alpha, theta, sigma, sigmae))) { # case 9 and 10
#     warning('Some parameters are NA or NaN')
#     ifelse(log, -Inf, 0)
#   } 
#   
#   if((alpha > 0 & sigma == 0 & sigmae == 0) |  # case 4
#      all(c(alpha, sigma, sigmae) == 0) ) {        # case 8
#     ifelse(log, -Inf, 0)
#   } else {
#     logpi <- log(pi)
#     loge2 <- log(2)
#     
#     if(as.integer(usempfr) >= 2 |
#        (usempfr & ((as.double(alpha) != 0 & as.double(alpha) < 0.01) |
#                    (as.double(sigma) != 0 & as.double(sigma) < 0.01) |
#                    (as.double(sigmae) != 0 & as.double(sigmae) < 0.01)))) {
#       usempfr <- 2 # indicate for later code that we shall convert to mpfr!
#     }
#     
#     done <- FALSE
#     # this loop governs the use of mpfr
#     while(!done & usempfr <= maxmpfr) { 
#       if(availRmpfr & usempfr >= 2) {
#         if(is.double(alphaorig)) {
#           alpha <- Rmpfr::mpfr(alphaorig, precbits*2^(usempfr-2))
#         }
#         if(is.double(thetaorig)) {
#           theta <- Rmpfr::mpfr(thetaorig, precbits*2^(usempfr-2))
#         }
#         if(is.double(sigmaorig)) {
#           sigma <- Rmpfr::mpfr(sigmaorig, precbits*2^(usempfr-2))
#         }
#         if(is.double(sigmaeorig)) {
#           sigmae <- Rmpfr::mpfr(sigmaeorig, precbits*2^(usempfr-2))
#         }
#         if(is.finite(g0)) {
#           g0 <- Rmpfr::mpfr(g0orig, precbits*2^(usempfr-2))
#         }
#         if(is.double(z)) {
#           z <- Rmpfr::mpfr(z, precbits*2^(usempfr-2))
#         }
#       }
#       
#       if(is.list(g0Prior)) {
#         if(is.null(g0Prior$mean) | is.null(g0Prior$var)) {
#           stop("g0Prior must have a 'mean' and a 'var' character members.")
#         } else {
#           g0Mean <- eval(parse(text=g0Prior$mean))
#           g0Var <- eval(parse(text=g0Prior$var))
#         }
#       } else {
#         g0Mean <- NA
#         g0Var <- NA
#       }
#       
#       N <- length(tree$tip.label)                # number of tips
#       
#       if(is.null(pruneInfo)) {
#         pruneInfo <- pruneTree(tree)
#       }
#       
#       # use is.double(alpha) to check whether Rmpfr is being used
#       
#       edge <- tree$edge
#       
#       M <- pruneInfo$M
#       endingAt <- pruneInfo$endingAt
#       nodesVector <- pruneInfo$nodesVector
#       nodesIndex <- pruneInfo$nodesIndex
#       nLevels <- pruneInfo$nLevels
#       unVector <- pruneInfo$unVector
#       unIndex <- pruneInfo$unIndex
#       unJ <- 1
#       
#       t <- tree$edge.length
#       
#       alphasigma2 <- alpha/sigma/sigma
#       theta2 <- theta^2
#       sigma2 <- sigma^2
#       sigmae2 <- sigmae^2
#       logsigma <- log(sigma)
#       logsigmae <- log(sigmae)
#       
#       talpha <- t*alpha
#       etalpha <- exp(talpha)
#       if(alpha != 0) {
#         fetalpha <- alpha/(1 - etalpha)
#       } else {
#         fetalpha <- -1/t
#       }
#       
#       e2talpha <- etalpha*etalpha
#       limitfe2talpha <- -0.5 / t
#       if(alpha != 0) {
#         fe2talpha <- alpha / (1 - e2talpha)
#       } else {
#         fe2talpha <- limitfe2talpha
#       }
#       
#       if(availRmpfr & usempfr & usempfr < maxmpfr) {
#         if(any(fe2talpha >= 0) | any(fe2talpha < limitfe2talpha)) {
#           usempfr <- usempfr + 1
#           next 
#         }
#       }
#       
#       log_fe2talpha <- log(-fe2talpha)
#       
#       fe2talphasigma2 <- fe2talpha / sigma2
#       
#       # for sigmae=0
#       r0 <- talpha + 0.5*log_fe2talpha - 0.5*logpi - logsigma
#       
#       # matrix for paramsIntegForks one pif for every node
#       # the following code creates a matrix for the class of alpha,
#       # i.e. could be mpfr as well as double
#       pif <- rep(alpha*0, M*3)
#       dim(pif) <- c(M, 3)
#       
#       
#       for(i in 1:nLevels) {
#         es <- endingAt[nodesVector[(nodesIndex[i]+1):nodesIndex[i+1]]]
#         
#         edgeEnds <- edge[es, 2]
#         
#         if(edge[es[1], 2] <= N) {
#           # all es pointing to tips
#           if(sigmae==0) {
#             # no environmental deviation
#             g1 <- z[edgeEnds]  # tip values
#             etalphag1thetatheta <- etalpha[es] * (g1 - theta) + theta
#             
#             pif[edgeEnds, 1] <- fe2talphasigma2[es]
#             pif[edgeEnds, 2] <- -2 * etalphag1thetatheta * fe2talphasigma2[es]
#             pif[edgeEnds, 3] <- r0[es] + etalphag1thetatheta^2 * fe2talphasigma2[es]
#           } else {
#             z1 <- z[edgeEnds]  # tip values
#             
#             u <- -0.5/sigmae2
#             v <- z1/sigmae2
#             w <- -0.5 * loge2 - 0.5 * logpi  - z1^2 / (2 * sigmae2) - logsigmae 
#             
#             gutalphasigma2 <- (fe2talpha[es] - alpha + u * sigma2) / fe2talpha[es]
#             
#             pif[edgeEnds, 1] <- u / gutalphasigma2
#             pif[edgeEnds, 2] <- 
#               (etalpha[es] * (2 * theta * u + v) - 2 * theta * u) / gutalphasigma2
#             pif[edgeEnds, 3] <- -0.5 * log(gutalphasigma2) -
#               0.25 * sigma2 * v^2 / (fe2talpha[es] - alpha + sigma2 * u) +
#               alpha * theta * (theta * u - etalpha[es] * (theta * u + v)) /
#               ((1 + etalpha[es]) * (sigma2 * u - alpha) + fetalpha[es]) +
#               talpha[es] + w
#           }
#         } else {
#           # edges pointing to internal nodes, for which all children nodes have been visited
#           u <- pif[edgeEnds, 1]
#           v <- pif[edgeEnds, 2]
#           w <- pif[edgeEnds, 3]
#           
#           gutalphasigma2 <- (fe2talpha[es] - alpha + u * sigma2) / fe2talpha[es]
#           
#           pif[edgeEnds, 1] <- u / gutalphasigma2
#           pif[edgeEnds, 2] <- 
#             (etalpha[es] * (2 * theta * u + v) - 2 * theta * u) / gutalphasigma2
#           pif[edgeEnds, 3] <- -0.5 * log(gutalphasigma2) -
#             0.25 * sigma2 * v^2 / (fe2talpha[es] - alpha + sigma2 * u) +
#             alpha * theta * (theta * u - etalpha[es] * (theta * u + v)) /
#             ((1 + etalpha[es]) * (sigma2 * u - alpha) + fetalpha[es]) +
#             talpha[es] + w
#         } 
#         
#         #update parent pifs
#         while(length(es)>0) {
#           un <- unVector[(unIndex[unJ]+1):unIndex[unJ+1]]
#           unJ <- unJ+1
#           
#           pif[edge[es[un], 1], ] <- 
#             pif[edge[es[un], 1], ] + pif[edge[es[un], 2], ]
#           
#           es <- es[-un]
#         }
#       }
#       
#       # at the root: P[Vtips|tree, g_root, alpha, theta, sigma, sigmae] =
#       #                  abc[1]*g_root^2 + abc[2]*g_root + abc[3]
#       abc <- pif[N+1, ]
#       
#       if(availRmpfr & usempfr & (any(is.nan(abc)) | abc[1]==0 | is.infinite(abc[3]))) {
#         usempfr <- usempfr + 1
#         next 
#       }
#       done <- TRUE
#     }
#     
#     # Processing of the root value
#     if(is.finite(g0)) {
#       # fixed g0 specified as parameter; 
#       
#       # pdf of the data conditioned on fixed g0
#       loglik <- abc[1]*g0^2 + abc[2]*g0 + abc[3]
#       
#       if(is.list(g0Prior)) {
#         # Normal prior for g0 with mean g0Mean and variance g0Var
#         g0LogPrior <- dnorm(as.double(g0), mean = as.double(g0Mean), sd = sqrt(as.double(g0Var)), log = log)
#       } else {
#         # No prior for g0 specified
#         g0LogPrior <- NA
#       }
#     } else if(is.na(g0) & !is.nan(g0)) {
#       # find g0 that would maximise one of the following:
#       # (1) if a normal prior for g0 was specified, 
#       #     pdf(z | alpha, theta, sigma, sigmae, g0, tree) x prior(g0)
#       # (2) otherwise, pdf(z | alpha, theta, sigma, sigmae, g0, tree) 
#       if(is.list(g0Prior)) {
#         # (1) Normal prior for g0 with mean g0Mean and variance g0Var
#         # max(loglik + g0Prior): derivative : 
#         # 2 * (abc[1] - 1 / (2 * g0Var)) * g0 + (abc[2] + g0Mean / g0Var) == 0
#         g0 <- -0.5 * (abc[2] + g0Mean / g0Var) / (abc[1] - 1 / (2 * g0Var))
#         loglik <- abc[1] * g0^2 + abc[2] * g0 + abc[3]
#         g0LogPrior <- dnorm(as.double(g0), mean = as.double(g0Mean), sd = as.double(sqrt(g0Var)), log = log)
#       } else {
#         # (2) No prior for g0 specified
#         # max(loglik): derivative : 2 * abc[1] * g0 + abc[2] == 0
#         g0 <- -0.5 * abc[2] / abc[1]
#         loglik <- abc[1] * g0^2 + abc[2] * g0 + abc[3]
#         g0LogPrior <- NA
#       }
#     } else {
#       # g0 is NaN or infinite. In this case, 
#       # we treat g0 as not specified, so the only option is to integrate under 
#       # its prior distribution. If g0Prior is not specified, assume that g0Prior 
#       # is the stationary OU normal distribution with mean, theta, and variance, 
#       # varOU(Inf, alpha, sigma)
#       if(!is.list(g0Prior)) {
#         g0Mean <- theta
#         g0Var <- varOU(Inf, alpha, sigma)
#       }
#       
#       if(g0Var != 0) {
#         def <- c(-0.5 / g0Var, 
#                  g0Mean / g0Var, 
#                  -0.5 * g0Mean^2 / g0Var - log(sqrt(g0Var)) - 0.5 * (loge2+logpi))
#         abc <- abc + def
#         loglik <- abc[3] + 0.5 * (logpi - log(-abc[1])) - abc[2]^2 / (4 * abc[1])
#         g0 <- NaN
#         g0LogPrior <- NaN
#       } else {
#         #g0Var = 0, i.e. gr=g0Mean with probability 1
#         loglik <- abc[1]*g0Mean[1]^2 + abc[2]*g0Mean[1] + abc[3]
#         g0 <- g0Mean
#         g0LogPrior <- Inf
#       }
#     }
#     
#     if(debug) {
#       debugdata <- data.table(
#         alpha = list(alpha), theta = list(theta), sigma = list(sigma), 
#         sigmae = list(sigmae),
#         alphasigma2 = list(alphasigma2), theta2 = list(theta2), 
#         sigma2 = list(sigma2),
#         sigmae2 = list(sigmae2), logsigma = list(logsigma), 
#         logsigmae = list(logsigmae), 
#         t = list(t),
#         talpha = list(talpha), etalpha = list(etalpha), 
#         e2talpha = list(e2talpha), 
#         fe2talpha = list(fe2talpha), log_fe2talpha = list(log_fe2talpha), 
#         r0 = list(r0), pif = list(pif), abc = list(abc), 
#         g0 = list(g0), g0Mean = list(g0Mean), g0Var = list(g0Var),
#         g0LogPrior = list(g0LogPrior), loglik = list(loglik), 
#         availRmpfr = list(availRmpfr), usempfr = list(usempfr), 
#         precbits = list(precbits))
#       print(debudata)
#     }
#     
#     if(is.double(alphaorig) & is.double(thetaorig) & is.double(sigmaorig) & 
#        is.double(sigmaeorig) & is.double(vorig)) {
#       loglik <- as.double(loglik)
#       g0 <- as.double(g0)
#       g0LogPrior <- as.double(g0LogPrior)
#     }
#     
#     value <- if(log) loglik else exp(loglik)
#     attr(value, "g0") <- g0
#     attr(value, "g0LogPrior") <- g0LogPrior
#     value
#   }
# }
