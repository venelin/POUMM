#' Random generation of values along a phylogenetic tree following a branching
#' OU process
#' 
#' @param z0 Numeric value indicating the initial state at the root of the tree.
#' @param alpha,theta,sigma Numeric values, parameters of the OU process.
#' @param sigmae Numeric non-negative value (default 0). Specifies the standard
#'   deviation of random environmental contribution (white noise) to be added to
#'   the values. 
#' @param tree An object of class phylo (package ape).
#' 
#' @return A numeric vector containing the generated values at all nodes
#'   (internal and tips) of the phylogenetic tree. 
#' @examples    
#'   N <- 20 
#'   tree <- TreeSim::sim.bd.taxa(N, 1, lambda = 2, mu = 1, complete =
#'   FALSE)[[1]] 
#'   z <- rVNodesGivenTreePOUMM(tree, 8, 2, 2, 1)
#'   
#'   phytools::plotBranchbyTrait(tree, z, mode = "nodes", show.tip.label =
#'   FALSE, show.node.label = FALSE) 
#'   ape::nodelabels(round(z[-(1:N)], 2)) 
#'   ape::tiplabels(round(z[(1:N)], 2))
#' @export
rVNodesGivenTreePOUMM <- function(tree, z0, alpha, theta, sigma, sigmae = 0) {
  g <- ape::rTraitCont(tree, 
                       function(x, l, .a, .t, .s) {
                         POUMM::rOU(n = 1, z0 = x, t = l, alpha = .a, theta = .t, sigma = .s)
                       }, root.value = z0, ancestor = TRUE, 
                       .a = alpha, .t = theta, .s = sigma)
  if(sigmae > 0) {
    g <- g + rnorm(length(g), 0, sigmae) 
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
#' @param sigmae Numeric non-negative value (default 0). Specifies the standard 
#'   deviation of random environmental contribution (white noise) to be added to
#'   the values.
#' @param e Numeric vector of size length(z) representing exactly known 
#'   measurement errors. Defaults to a vector of zeroes.
#' @param log Logical indicating whether a log-likelihood should be returned 
#'   instead of a likelihood. Default is TRUE.
#'   
#' @return A numeric value, the multivariate probability density of z under the 
#'   given parameters.
#' @export
dVNodesGivenTreePOUMM <- function(z, tree, alpha, theta, sigma, sigmae=0,
                                  e = rep(0, length(z)), log=TRUE) {
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
#' @param usempfr,maxmpfr integer indicating if and how mpfr should be used for
#'   small parameter values (any(c(alpha, sigma, sigmae) < 0.01)). Using the
#'   mpfr package can be forced by specifying an integer greater or equal to 2.
#'   Setting usempfr=0 disables high precision likelihood calculation. Requires
#'   the Rmpfr package. Note that when using mpfr, the time for one likelihood
#'   calculation can increase more than 100-fold. Default (0).
#' @param precbits integer specifying precision bits for mpfr. Default is 512.
#' @param debug logical, if set to TRUE some debugging information is printed 
#'   during likelihood calculation
#'   
#' @import data.table
#' @aliases likPOUMMGivenTreeVTips dVTipsGivenTreePOUMM
#' @export
likPOUMMGivenTreeVTips <- dVTipsGivenTreePOUMM <- function(
  z, tree, alpha, theta, sigma, sigmae = 0, g0 = NA, g0Prior = NULL,
  log = TRUE, pruneInfo = pruneTree(tree), 
  usempfr = 0, maxmpfr = 2, precbits = 128, backend = "R",
  debug = FALSE) {
  
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
  } 
  
  if((alpha > 0 & sigma == 0 & sigmae == 0) |  # case 4
     all(c(alpha, sigma, sigmae) == 0) ) {        # case 8
    ifelse(log, -Inf, 0)
  } else {
    logpi <- log(pi)
    loge2 <- log(2)
    
    if(as.integer(usempfr) >= 2 |
       (usempfr & ((as.double(alpha) != 0 & as.double(alpha) < 0.01) |
                   (as.double(sigma) != 0 & as.double(sigma) < 0.01) |
                   (as.double(sigmae) != 0 & as.double(sigmae) < 0.01)))) {
      usempfr <- 2 # indicate for later code that we shall convert to mpfr!
    }
    
    
    N <- length(tree$tip.label)                # number of tips
    
    if(is.null(pruneInfo)) {
      pruneInfo <- pruneTree(tree)
    }
    
    edge <- tree$edge
    
    M <- pruneInfo$M
    endingAt <- pruneInfo$endingAt
    nodesVector <- pruneInfo$nodesVector
    nodesIndex <- pruneInfo$nodesIndex
    nLevels <- pruneInfo$nLevels
    unVector <- pruneInfo$unVector
    unIndex <- pruneInfo$unIndex
    unJ <- 1
    
    t <- tree$edge.length
    
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
      
      
      # use is.double(alpha) to check whether Rmpfr is being used
      
      sigma2 <- sigma^2
      sigmae2 <- sigmae^2
      logsigma <- log(sigma)
      logsigmae <- log(sigmae)
      
      talpha <- t*alpha
      etalpha <- exp(talpha)
      if(alpha != 0) {
        fetalpha <- alpha/(1 - etalpha)
      } else {
        fetalpha <- -1/t
      }
      
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
      
      if(backend == "C" & usempfr < 2) {
        # matrix for paramsIntegForks one pif for every node
        # the following code creates a matrix for the class of alpha,
        # i.e. could be mpfr as well as double
        pif <- mainLoopLikPOUMM_C(
          nLevels, M, N, z, 
          alpha, theta, sigma, sigmae, logsigmae,
          sigma2, logsigma, sigmae2, loge2, logpi,
          talpha, etalpha, e2talpha, fe2talpha,
          edge, endingAt, 
          nodesVector, nodesIndex, 
          unVector, unIndex)
      } else {
        pif <- rep(alpha*0, M*3)
        dim(pif) <- c(M, 3)
        
        gutalphasigma2 <- e2talpha + ((-0.5 / sigmae2) * sigma2) / fe2talpha
        
        
        for(i in 1:nLevels) {
          es <- endingAt[nodesVector[(nodesIndex[i] + 1):nodesIndex[i + 1]]]
          edgeEnds <- edge[es, 2]
          
          if(edge[es[1], 2] <= N) {
            # all es pointing to tips
            if(sigmae==0) {
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
                
                #pif[edgeEnds, 3] <- -0.5 * (2*talpha[es] + 
                #     log1p(((-0.5 / sigmae2) * sigma2) / fe2talpha[es] / e2talpha[es])) -
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
            un <- unVector[(unIndex[unJ]+1):unIndex[unJ+1]]
            unJ <- unJ+1
            lenUnAll <- lenUnAll + length(un)
            
            pif[edge[es[un], 1], ] <- 
              pif[edge[es[un], 1], ] + pif[edge[es[un], 2], ]
            
            #es <- es[-un]
          }
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
      
      if(g0Var != 0) {
        def <- c(-0.5 / g0Var, 0, - log(sqrt(g0Var)) - 0.5 * (loge2+logpi))
        abc <- abc + def
        loglik <- abc[3] + 0.5 * (logpi - log(-abc[1])) - abc[2]^2 / (4 * abc[1])
        g0 <- NaN
        g0LogPrior <- NaN
      } else {
        #g0Var = 0, i.e. gr=g0Mean with probability 1
        loglik <- (abc[1] * (g0Mean[1]-theta) + abc[2]) * (g0Mean[1]-theta) + abc[3]
        
        g0 <- g0Mean
        
        g0LogPrior <- Inf
      }
    }
    
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
        g0 = list(g0), g0Mean = list(g0Mean), g0Var = list(g0Var),
        g0LogPrior = list(g0LogPrior), loglik = list(loglik), 
        availRmpfr = list(availRmpfr), usempfr = list(usempfr), 
        precbits = list(precbits))
      attr(value, "debugdata") <- debugdata
    }
    
    value
  }
}

#' An "algebraic" implementation of the POUMM likelihood calculation based on
#' multivariate normal density function
#' This is an internal function used for tests only.
#' 
#' @importFrom mvtnorm dmvnorm
dVTipsGivenTreePOUMMg0Alg <- function(z, tree, alpha, theta, sigma, sigmae, g0) {
  tanc <- vcv(tree)
  tipTimes <- diag(tanc)
  N <- length(tree$tip.label)
  
  D0 <- sapply(tipTimes, function(t) exp(-alpha*t))
  D1 <- 1-D0
  Vt <- covVTipsGivenTreePOUMM(tree, alpha, sigma, 0, tanc)
  Ve <- diag(sigmae^2, N, N)
  
  muz <- D0 * g0  +  D1 * theta
  
  dmvnorm(z, muz, Vt+Ve, log=TRUE)
}