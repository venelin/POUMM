# Implementation of the POUMM method

NULL
                       
#' Variance covariance matrix of tips in an OUTree
#' 
#' @param tree a phylo object
#' @param alpha a non-negative numeric, default is 0
#' @param sigma a non-negative numeric, default is 1
#' @param sigmae a non-negative numeric, default is 0
#' @param corr logical indicating if a correlation matrix shall be returned
#' @param tanc numerical matrix with the time-distance from the root of the 
#' tree to the mrca of each tip-couple. If NULL it will be calculated.
#' 
#' @return a variance covariance or a correlation matrix of the tips in tree.
#' @references (Hansen 1997) Stabilizing selection and the comparative analysis
#'   of adaptation.
#' @export
covOUTree <- function(tree, alpha=0, sigma=1, sigmae=0, corr=F, tanc=NULL) {
  N <- length(tree$tip.label)
  
  # distances from the root to the mrca's of each tip-couple.
  if(is.null(tanc)) {
    tanc <- vcv(tree)
  }
  if(alpha > 0) {
    varanc <- 0.5*(1-exp(-2*alpha*tanc))/alpha
  } else {
    # limit case alpha==0
    varanc <- tanc
  }
  
  # time from the root to each tip
  ttips <- diag(tanc) # nodeTimes(tree)[1:N]
  # distance-matrix between tips
  tij <- sapply(ttips, function(.) .+ttips)-2*tanc
  
  covij <- exp(-alpha*tij)*varanc*sigma^2 + diag(rep(sigmae^2, N))
  if(corr) {
    sqvar <- sqrt(diag(covij))
    sccor <- sapply(sqvar, function(.) .*sqvar)
    covij/sccor
  } else {
    covij
  }  
}

#' Likelihood of observed node values given a tree
#' 
#' @description
#' Calculates the conditional probability density of the tree 
#' trait values at the tips and internal nodes of the tree 
#' given the tree and the trait-value at the root according to a model, such 
#' that the probability density function of a value 
#' at a descendent node gdesc given its ancestor-value ganc and edge-length l 
#' is given by pdf(gdesc, ganc, l, ...). 
#
#' @param g A numeric vector of size length(tree$tip.label)+tree$Nnode 
#'    representing the genetic contibution to the trait data at tips, root and 
#'    internal nodes. 
#' @param tree an object of class phylo
#' @param pdf A probability density function
#' @param log Logical indicating whether a log-likelihood should be returned 
#'    instead of a likelihood. Default is TRUE.
#' @param ... additional arguments passed to the pdf function
#' @return
#'    If !nodeprobs (default) log P(tree, g | pdf), otherwise the vector of the
#'    node-value's (log-)probabilities conditioned on their ancestor's values 
#'    otherwise (root not included)
#' @export
likGTree <- function(g, tree, pdf, log=T, ...) {
  
  if (is.null(tree$edge.length)) 
    stop("tree has no branch length")
  if (any(tree$edge.length < 0)) 
    stop("at least one branch length negative")
  
  N <- dim(tree$edge)[1]            # number of edges in the tree
  anc <- tree$edge[, 1]
  des <- tree$edge[, 2]
  el <- tree$edge.length
  probs <- sapply(1:N, function(i) pdf(g[des[i]], g[anc[i]], el[i], ..., log=T))
  
  ifelse(log, sum(probs), exp(sum(probs)))
}

#' Likelihood of observed genetic contributions at all nodes and environmental 
#' deviations at the tips given a tree.
#' 
#' @description
#' Calculates the probability density function of the genetic trait 
#' contributions at all nodes, g, and environmental 
#' trait deviations at the tip-nodes given the tree topology and assuming that 
#' the genetic contributions evolved on the tree according to an OU process 
#' with parameters alpha, theta and sigma and the environmental deviations 
#' represent iid samples from a normal distribution with mean mue and standard 
#' deviations sigmae.
#'
#' @param g A numeric vector of size length(tree$tip.label)+tree$Nnode 
#'    representing the genetic contibution to the trait data at tips, root and 
#'    internal nodes. Note that g[length(tree$tip.label)+1] is assumed to be 
#'    the root.
#' @param tree an object of class phylo
#' @param alpha the strength of the selective constraint
#' @param theta long term mean value of the OU process
#' @param sigma the unit-time standard deviation of the random component in the OU
#'     process.
#' @param e NULL or a numeric vector of size length(tree$tip.label) representing the 
#'    environmental deviations at the tip-nodes. If NULL the arguments mue and
#'    sigmae are ignored.
#' @param mue Mean of the environmental deviation distribution.
#' @param sigmae standard deviation of the environmental deviation distribution.
#' @param log Logical indicating whether a log-likelihood should be returned 
#'    instead of a likelihood. Default is TRUE.
#' @param ...  additional arguments passed to the pdf function
#' @return
#'    (log-)probability of g and e given tree, alpha, theta, sigma, mue, sigmae.
#' @export
likGETreeOU <- function(g, tree, alpha, theta, sigma, 
                        e=NULL, mue=0, sigmae=1, log=T) {
  lTG <- likGTree(g, tree, dCondOU, alpha=alpha, theta=theta, sigma=sigma, log=T)
  lRoot <- dStOU(g[length(tree$tip.label)+1], alpha, theta, sigma, log=T)
  lE <- 0
  if(!is.null(e)) 
    lE <- dnorm(e, mue, sigmae, log=T)
  
  ifelse(log, sum(lTG, lRoot, lE), exp(sum(lTG, lRoot, lE)))
}

#' Likelihood of observed tip-values given a tree, assuming Ornstein-Uhlenbeck
#' process for the genetic contributions along the tree and normally distributed
#' environmental deviations.
#' 
#' @description
#' Calculates the (log-)probability density of trait values at the tip-nodes 
#' given the tree and assuming that the trait value at the tips is the sum 
#' of a genetic contribution, g, that evolved on the tree according to an OU 
#' process with parameters alpha, theta, sigma and an environmental deviation,
#' e, that is distributed normally and independently between the tips of the 
#' tree.
#' Note: Without additional assumptions for the distribution of the value at 
#' the root of the tree, the likelihood is not defined at alpha=0, although 
#' this corresponds to the limiting Brownian motion process with mean value
#' theta and unit time variance sigma^2.
#'
#' @param tree an object of class phylo
#' @param v A numeric vector of size length(tree$tip.label) representing the trait
#'     values at the tip-nodes.
#' @param alpha the strength of the selection
#' @param theta long term mean value of the OU process
#' @param sigma the unit-time standard deviation of the random component in the OU
#'     process.
#' @param sigmae the standard deviation of the environmental deviation added to the
#'     genetic contribution at each tip, by default 0, meaning no environmental
#'     deviation.
#' @param distgr a character or a numeric describing the distribution or fixed value 
#'     to be assumed for the genetic contribution at the root. If a character the 
#'     acceptable values are : 
#'       'normal'- normal distribution with mean mugr and standard deviation sigmagr; 
#'       'maxlik' - fixed value that would maximise the conditional likelihood on gr.
#' @param mugr see distgr
#' @param sigmagr see distgr
#' @param log Logical indicating whether log-likelihood should be returned instead
#'    of likelihood, default is TRUE.
#' @param impl Character indicating which implementation of the method should 
#' be used. One of  'R1', 'R2', 'R4'. Defaults to 'R4'. 
#' @param pruneInfo list returned by prune(tree) to be passed when impl='R5'. 
#' @param usempfr,maxmpfr integer indicating if and how mpfr should be used for small parameter values 
#'    (any(c(alpha, sigma, sigmae) < 0.01)). Using the mpfr package can be forced by specifying an 
#'    integer greater or equal to 2. Setting usempfr=0 disables high precision likelihood calculation.
#'    Requires the Rmpfr package. Note that when using mpfr, the time for one likelihood calculation can increase
#'    substantially. 
#'    with the use of mpfr. 
#' @param precbits integer specifying precision bits for mpfr. Default is 512.
#' @param debug logical, if set to TRUE some debugging information is stored in a global 
#'    list called .likVTreeOUDebug (currently implemented only for impl='R3')
#' @export
likVTreeOU <- function(v, tree, alpha, theta, sigma, sigmae=0, lambda=1,
                       distgr=c('normal', 'maxlik'), 
                       mugr=theta, 
                       sigmagr=ifelse(alpha==0&sigma==0, 0, sigma/sqrt(2*alpha)),
                       log=T, impl=c('R4', 'R2', 'R1'), pruneInfo=NULL,
                       usempfr=1, maxmpfr=2, precbits=512, debug=F) {
  alphaorig <- alpha
  thetaorig <- theta
  sigmaorig <- sigma
  sigmaeorig <- sigmae
  mugrorig <- mugr
  sigmagrorig <- sigmagr
  vorig <- v
  
  if(any(tree$edge.length<=0)) {
    stop('All branch lengths in tree should be positive!')
  }
  if(any(is.na(c(alpha, theta, sigma, sigmae))) | 
    any(is.nan(c(alpha, theta, sigma, sigmae, lambda)))) {
    warning('Some parameters are NA or NaN')
    ifelse(log, -Inf, 0)
  } 
  if(any(is.infinite(c(alpha, sigma, sigmae))) | # case 9
       any(c(alpha, sigma, sigmae) < 0) |        # case 10
       (alpha > 0 & sigma == 0 & sigmae == 0) |  # case 4
       all(c(alpha, sigma, sigmae) == 0)        # case 8
       ) {
    ifelse(log, -Inf, 0)
  } else {
    logpi <- log(pi)
    loge2 <- log(2)
    
    availRmpfr <- require(Rmpfr) 
    
    if(as.integer(usempfr)>=2 |
         (usempfr & ((as.double(alpha) != 0 & as.double(alpha) < 0.01) |
                       (as.double(sigma) != 0 & as.double(sigma) < 0.01) |
                       (as.double(sigmae) != 0 & as.double(sigmae) < 0.01)))) {
      usempfr <- 2 # indicate for later code that we shall convert to mpfr!
    }
     
    done <- F
    # this loop is global for all implementations and only governs the use of mpfr
    while(!done & usempfr <= maxmpfr) { 
      if(availRmpfr & usempfr >= 2) {
        if(is.double(alphaorig)) alpha <- mpfr(alphaorig, precbits*2^(usempfr-2))
        if(is.double(thetaorig)) theta <- mpfr(thetaorig, precbits*2^(usempfr-2))
        if(is.double(sigmaorig)) sigma <- mpfr(sigmaorig, precbits*2^(usempfr-2))
        if(is.double(sigmaeorig)) sigmae <- mpfr(sigmaeorig, precbits*2^(usempfr-2))
        if(is.double(mugr)) mugr <- mpfr(mugr, precbits*2^(usempfr-2))
        if(is.double(sigmagr)) sigmagr <- mpfr(sigmagr, precbits*2^(usempfr-2))
        if(is.double(v)) v <- mpfr(v, precbits*2^(usempfr-2))
      }
      
      if(impl[1] == 'R4') {
        N <- length(tree$tip.label)                # number of tips
        M <- length(unique(as.vector(tree$edge)))  # number of all nodes 
        endingAt <- order(rbind(tree$edge, c(0, N+1))[, 2])
        
        if(lambda != 1) {
          # transform branch lengths
          times <- nodeTimes(tree)
          tree$edge.length <- tree$edge.length*lambda
          timeslambda <- nodeTimes(tree)
          tree$edge.length[endingAt[1:N]] <- 
            tree$edge.length[endingAt[1:N]] + times[1:N] - timeslambda[1:N]
          if(lambda==0) {
            #star topology
            tree$edge.length <- tree$edge.length[endingAt[1:N]]
            tree$edge <- tree$edge[endingAt[1:N], ]
            tree$edge[, 1] <- N+1
            M <- N+1
            endingAt <- order(rbind(tree$edge, c(0, N+1))[, 2])
          }
        }
        
        edge <- tree$edge
        t <- tree$edge.length
        
        alphasigma2 <- alpha/sigma/sigma
        theta2 <- theta^2
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
        limitfe2talpha <- -0.5/t
        if(alpha != 0) {
          fe2talpha <- alpha/(1 - e2talpha)
        } else {
          fe2talpha <- limitfe2talpha
        }
        
        if(availRmpfr & usempfr & usempfr < maxmpfr) {
          if(any(fe2talpha >= 0) | any(fe2talpha < limitfe2talpha)) {
            usempfr <- usempfr + 1
            next 
          }
        }
        
        log_fe2talpha <- log(-fe2talpha)

        fe2talphasigma2 <- fe2talpha/sigma2
        
        # for sigmae=0
        r0 <- talpha + 0.5*log_fe2talpha - 0.5*logpi - logsigma
        
        nonVisitedChildren <- rep(0, M)
        ee1 <- edge[, 1]
        while(length(ee1)) {
          matchp <- match((N+1):M, ee1)
          matchp <- matchp[!is.na(matchp)]
          nonVisitedChildren[ee1[matchp]] <- nonVisitedChildren[ee1[matchp]] + 1
          ee1 <- ee1[-matchp]
        }
        
        # matrix for paramsIntegForks one pif for every node
        # the following code creates a matrix for the class of alpha,
        # i.e. could be mpfr as well as double
        pif <- rep(alpha*0, M*3)
        dim(pif) <- c(M, 3)
        
        # start from the edges leading to tips
        nodes <- 1:N
        
        while(nodes[1] != N+1) {
          es <- endingAt[nodes]
          nodes <- c()
          edgeEnds <- edge[es, 2]
          
          if(edge[es[1], 2] <= N) {
            # all es pointing to tips
            if(sigmae==0) {
              # no environmental deviation
              g1 <- v[edgeEnds]  # tip values
              etalphag1thetatheta <- etalpha[es]*(g1-theta)+theta
              
              pif[edgeEnds, ]  <- 
                c(fe2talphasigma2[es],
                  -2 * etalphag1thetatheta * fe2talphasigma2[es],
                  r0[es] + etalphag1thetatheta^2 * fe2talphasigma2[es]
                )
            } else {
              v1 <- v[edgeEnds]  # tip values
              
              u <- -0.5/sigmae2
              w <- v1/sigmae2
              z <- -0.5*loge2 - 0.5*logpi  - v1^2/(2*sigmae2) - logsigmae 
              
              gutalphasigma2 <- (fe2talpha[es]-alpha+u*sigma2)/fe2talpha[es]
              
              pif[edgeEnds, ] <- c(u/gutalphasigma2,
                                      (etalpha[es]*(2*theta*u+w)-2*theta*u)/gutalphasigma2,
                                      -0.5*log(gutalphasigma2) -
                                        0.25*sigma2*w^2/(fe2talpha[es]-alpha + sigma2*u) + 
                                        alpha*theta*(theta*u-etalpha[es]*(theta*u+w))/((1+etalpha[es])*(sigma2*u-alpha)+fetalpha[es]) + 
                                        talpha[es]+z)
            }
          } else {
            # edges pointing to internal nodes, for which all children nodes have been visited
            u <- pif[edgeEnds, 1]
            w <- pif[edgeEnds, 2]
            z <- pif[edgeEnds, 3]
            
            gutalphasigma2 <- (fe2talpha[es]-alpha+u*sigma2)/fe2talpha[es]
            
            pif[edgeEnds, ] <- c(u/gutalphasigma2,
                                    (etalpha[es]*(2*theta*u+w)-2*theta*u)/gutalphasigma2,
                                    -0.5*log(gutalphasigma2) -
                                      0.25*sigma2*w^2/(fe2talpha[es]-alpha + sigma2*u) + 
                                      alpha*theta*(theta*u-etalpha[es]*(theta*u+w))/((1+etalpha[es])*(sigma2*u-alpha)+fetalpha[es]) + 
                                      talpha[es]+z)
          } 

          #update parent pifs
          while(length(es)>0) {
            un <- match(unique(edge[es, 1]), edge[es, 1])
            pif[edge[es[un], 1], ] <- pif[edge[es[un], 1], ] + pif[edge[es[un], 2], ]
            nonVisitedChildren[edge[es[un], 1]] <- nonVisitedChildren[edge[es[un], 1]] - 1
            nodes <- c(nodes, edge[es[un][nonVisitedChildren[edge[es[un], 1]] == 0], 1])
            es <- es[-un]
          }
        }
            
        # at the root: P[Vtips|tree, g_root, alpha, theta, sigma, sigmae] =
        #                  abc[1]*g_root^2 + abc[2]*g_root + abc[3]
        abc <- pif[N+1, ]

        if(availRmpfr & usempfr & (any(is.nan(abc)) | abc[1]==0 | is.infinite(abc[3]))) {
          usempfr <- usempfr + 1
          next 
        }
        done <- T
      } else if(impl[1]=='R5') {
        N <- length(tree$tip.label)                # number of tips
        
        if(is.null(pruneInfo)) {
          warning('You should supply pruneInfo with impl=="R5". see ?pruneInfo')
          pruneInfo <- pruneTree(tree)
        }
        
        if(lambda != 1) {
          stop('impl="R5" does not support lambda different than 1. Use impl="R4" instead.')
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
        
        alphasigma2 <- alpha/sigma/sigma
        theta2 <- theta^2
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
        limitfe2talpha <- -0.5/t
        if(alpha != 0) {
          fe2talpha <- alpha/(1 - e2talpha)
        } else {
          fe2talpha <- limitfe2talpha
        }
        
        if(availRmpfr & usempfr & usempfr < maxmpfr) {
          if(any(fe2talpha >= 0) | any(fe2talpha < limitfe2talpha)) {
            usempfr <- usempfr + 1
            next 
          }
        }
        
        log_fe2talpha <- log(-fe2talpha)
        
        fe2talphasigma2 <- fe2talpha/sigma2
        
        # for sigmae=0
        r0 <- talpha + 0.5*log_fe2talpha - 0.5*logpi - logsigma
        
        # matrix for paramsIntegForks one pif for every node
        # the following code creates a matrix for the class of alpha,
        # i.e. could be mpfr as well as double
        pif <- rep(alpha*0, M*3)
        dim(pif) <- c(M, 3)
        
        for(i in 1:nLevels) {
          nodes <- nodesVector[(nodesIndex[i]+1):nodesIndex[i+1]]
          es <- endingAt[nodes]
          nodes <- c()
          edgeEnds <- edge[es, 2]
          
          if(edge[es[1], 2] <= N) {
            # all es pointing to tips
            if(sigmae==0) {
              # no environmental deviation
              g1 <- v[edgeEnds]  # tip values
              etalphag1thetatheta <- etalpha[es]*(g1-theta)+theta
              
              pif[edgeEnds, ]  <- 
                c(fe2talphasigma2[es],
                  -2 * etalphag1thetatheta * fe2talphasigma2[es],
                  r0[es] + etalphag1thetatheta^2 * fe2talphasigma2[es]
                )
            } else {
              v1 <- v[edgeEnds]  # tip values
              
              u <- -0.5/sigmae2
              w <- v1/sigmae2
              z <- -0.5*loge2 - 0.5*logpi  - v1^2/(2*sigmae2) - logsigmae 
              
              gutalphasigma2 <- (fe2talpha[es]-alpha+u*sigma2)/fe2talpha[es]
              
              pif[edgeEnds, ] <- c(u/gutalphasigma2,
                                   (etalpha[es]*(2*theta*u+w)-2*theta*u)/gutalphasigma2,
                                   -0.5*log(gutalphasigma2) -
                                     0.25*sigma2*w^2/(fe2talpha[es]-alpha + sigma2*u) + 
                                     alpha*theta*(theta*u-etalpha[es]*(theta*u+w))/((1+etalpha[es])*(sigma2*u-alpha)+fetalpha[es]) + 
                                     talpha[es]+z)
            }
          } else {
            # edges pointing to internal nodes, for which all children nodes have been visited
            u <- pif[edgeEnds, 1]
            w <- pif[edgeEnds, 2]
            z <- pif[edgeEnds, 3]
            
            gutalphasigma2 <- (fe2talpha[es]-alpha+u*sigma2)/fe2talpha[es]
            
            pif[edgeEnds, ] <- c(u/gutalphasigma2,
                                 (etalpha[es]*(2*theta*u+w)-2*theta*u)/gutalphasigma2,
                                 -0.5*log(gutalphasigma2) -
                                   0.25*sigma2*w^2/(fe2talpha[es]-alpha + sigma2*u) + 
                                   alpha*theta*(theta*u-etalpha[es]*(theta*u+w))/((1+etalpha[es])*(sigma2*u-alpha)+fetalpha[es]) + 
                                   talpha[es]+z)
          } 
          
          #update parent pifs
          while(length(es)>0) {
            un <- unVector[(unIndex[unJ]+1):unIndex[unJ+1]]
            unJ <- unJ+1
            pif[edge[es[un], 1], ] <- pif[edge[es[un], 1], ] + pif[edge[es[un], 2], ]
            es <- es[-un]
          }
        }
        
        # at the root: P[Vtips|tree, g_root, alpha, theta, sigma, sigmae] =
        #                  abc[1]*g_root^2 + abc[2]*g_root + abc[3]
        abc <- pif[N+1, ]
        
        if(availRmpfr & usempfr & (any(is.nan(abc)) | abc[1]==0 | is.infinite(abc[3]))) {
          usempfr <- usempfr + 1
          next 
        }
        done <- T
      } else if(impl[1] == 'R1') {
        esubs <- edgesFrom(tree, length(tree$tip.label)+1)
        abc <- paramsForkOUR1(tree, v, esubs[1], alpha, theta, sigma, sigmae)
        for(es in esubs[-1]) 
          abc <- abc + paramsForkOUR1(tree, v, es, alpha, theta, sigma, sigmae)
        done <- T
      } else if(impl[1] == 'R2') {
        esubs <- edgesFrom(tree, length(tree$tip.label)+1)
        abc <- paramsForkOUR2(tree, v, esubs[1], alpha, theta, sigma, sigmae)
        for(es in esubs[-1]) 
          abc <- abc + paramsForkOUR2(tree, v, es, alpha, theta, sigma, sigmae)      
        done <- T
      } 
    }
    
    if(!is.character(distgr <- distgr[1])) {
      # fixed value
      loglik <- abc[1]*distgr[1]^2 + abc[2]*distgr[1] + abc[3]
    } else {        
      if(distgr[1] == 'normal') {
        mugr2 <- mugr^2
        if(sigmagr > 0) {
          sigmagr2 <- sigmagr^2
          def <- c(-0.5/sigmagr2, 
                   mugr/sigmagr2, 
                   -0.5*mugr2/sigmagr2-log(sigmagr)-0.5*(loge2+logpi))
          abc <- abc + def
          loglik <- abc[3] + 0.5*(log(pi)-log(-abc[1])) - abc[2]^2/(4*abc[1])
        } else {
          #sigmagr = 0, i.e. gr=mugr with probability 1
          loglik <- abc[1]*mugr[1]^2 + abc[2]*mugr[1] + abc[3]
        }
      } else if(distgr[1] == 'maxlik') {
        # maxlik: derivative : 2*abc[1]*gr+abc[2] == 0
        grmax <- -0.5*abc[2]/abc[1]
        loglik <- abc[1]*grmax^2 + abc[2]*grmax + abc[3]
      } 
    }  
    if(debug) {
      debugdata <- data.table(alpha=list(alpha), theta=list(theta), sigma=list(sigma), sigmae=list(sigmae),
                              lambda=list(lambda), 
                              alphasigma2=list(alphasigma2), theta2=list(theta2), sigma2=list(sigma2),
                              sigmae2=list(sigmae2), logsigma=list(logsigma), logsigmae=list(logsigmae), 
                              t=list(t),
                              talpha=list(talpha), etalpha=list(etalpha), e2talpha=list(e2talpha), 
                              fe2talpha=list(fe2talpha), log_fe2talpha=list(log_fe2talpha), 
                              r0=list(r0), pif=list(pif), abc=list(abc), 
                              distgr=list(distgr), mugr=list(mugr), sigmagr=list(sigmagr),
                              def=list(c(-0.5/sigmagr^2, 
                                         mugr/sigmagr^2, -0.5*mugr^2/sigmagr^2-log(sigmagr)-0.5*(loge2+logpi))),
                              loglik=list(loglik), grmax=list(-0.5*abc[2]/abc[1]), 
                              availRmpfr=list(availRmpfr), impl=list(impl[1]), usempfr=list(usempfr), precbits=list(precbits))
      if(!exists('.likVTreeOUDebug')) {
        .likVTreeOUDebug <<- debugdata
      } else {
        .likVTreeOUDebug <<- rbindlist(list(.likVTreeOUDebug, debugdata))
      }
    }

    if(is.double(alphaorig) & is.double(thetaorig) & is.double(sigmaorig) & 
         is.double(sigmaeorig) & is.double(vorig)) {
      loglik <- as.double(loglik)
    }
    
    value <- if(log) loglik else exp(loglik)
    if(exists('grmax'))
      attr(value, 'grmax') <- as.double(grmax)
    value
  }
}


#' Calculate alpha given heritability H2, sigma and sigmae
#' @export
calcAlpha <- function(H2, sigma, sigmae) {
  sigma^2*(1-H2)/(2*H2*sigmae^2)
}

#' Calculate Sigma given heritability at time, alpha and sigmae
#' @param H2 numeric, the heritability at time t
#' @param alpha numeric, selection strength of the OU process
#' @param sigmae numeric, environmental phenotypic deviation at the tips
#' @param t numeric, time since the beginning of the OU process, default Inf
#' 
#' 
#' @export
calcSigma <- function(H2, alpha, sigmae, t=Inf) {
  if(H2>1 | H2<0 | alpha < 0 | sigmae < 0 | t < 0) {
    stop("H2(t) should be in [0, 1], alpha, sigmae and t should be non-negative.")
  }
  if(is.infinite(t)) {
    if(alpha==0) {
      # BM
      sqrt(sigmae^2*H2/(t*(1-H2)))
    } else {
      # alpha>0, OU in equilibrium
      sqrt(2*alpha*H2*sigmae^2/(1-H2))
    }
  } else {
    # t is finite
    if(alpha==0) {
      # BM
      sqrt(sigmae^2*H2/(t*(1-H2)))
    } else {
      # alpha>0, OU not in equilibrium
      sqrt(2*alpha*H2*sigmae^2/((1-exp(-2*alpha*t))*(1-H2)))
    } 
  }
}

#' Calculate Sigmag (i.e. standard dev. of stationary distr.) given H2 and
#' sigmae
#' @export
calcSigmag <- function(H2, sigmae) {
  sqrt(H2*sigmae^2/(1-H2))
}

#' Calculate sigmae given alpha, sigma, and H2.
#' @details Uses the formula 
#' H2=varStOU(alpha, sigma)/(varStOU(alpha, sigma)+sigmae^2)
#' @export
calcSigmae <- function(alpha, sigma, H2) {
  sigmaStOU2 <- varStOU(alpha, sigma)
  sqrt(sigmaStOU2/H2 - sigmaStOU2)
}

#' Long-term (equilibrium) Broad-sense heritability H2 given alpha, sigma and sigmae
#' 
#' @note This function calculates the heritability by assuming
#' that all tips in the tree are sampled from the 
#' equilibrium distribution.
#' @export
calcH2 <- function(alpha, sigma, sigmae) {
  if(sigma==0) {
    if(sigmae==0) {
      NaN
    } else {
      # sigmae>0
      0
    }
  } else if(alpha==0) {
    1
  } else {
    sigma^2/(sigma^2+2*alpha*sigmae^2)  
  }
}

#' Broad-sense heritability estimated from the empirical variance of
#' the observed phenotypes and sigmae
#' @param v numerical vector of observed phenotypes
#' @param sigmae numerical standard deviation of the environmental deviation
#' @return numerical between 0 and 1 
#' @export
calcH2EmpSigmae <- function(v, sigmae, tree=NULL, tFrom=0, tTo=Inf) {
  if(!is.null(tree)) {
    tipTimes <- nodeTimes(tree)[1:length(tree$tip.label)]
    v <- v[which(tipTimes>=tFrom & tipTimes<=tTo)]
  }
  1-(sigmae^2/var(v))
}

#' Broad-sense heritability estimated at time t
#' @param alpha,sigma numeric, parameters of the OU process acting on the genetic contributions 
#' @param sigmae numeric, environmental standard deviation
#' @param t time from the beginning of the process at which heritability should be calculated, i.e.
#' epidemiologic time
#' @param tm average time for within host evolution from getting infected until getting measured or 
#' passing on the infection to another host
#' @export
calcH2AtTime <- function(alpha, sigma, sigmae, t, tm=0
                         #, as=c('vector', 'mean', 'median')
                         ) {
  lenoutput <- max(sapply(list(alpha, sigma, sigmae, t, tm), length))
  if(lenoutput>1) {
    alpha <- rep(alpha, lenoutput/length(alpha))
    sigma <- rep(sigma, lenoutput/length(sigma))
    sigmae <- rep(sigmae, lenoutput/length(sigmae))
    t <- rep(t, lenoutput/length(t))
    tm <- rep(tm, lenoutput/length(tm))
  }
  sigmag2 <- ifelse(alpha>0, 0.5*(1-exp(-2*alpha*t))/alpha*sigma^2, t*sigma^2)  
  sigmagm2 <- ifelse(alpha>0, 0.5*(1-exp(-2*alpha*tm))/alpha*sigma^2, tm*sigma^2)
  
  ifelse(t>tm, (sigmag2-sigmagm2)/(sigmag2+sigmae^2), rep(0, length(t)))
}

#' Genotypic variance at time
#' @param alpha,sigma numeric, parameters of the OU process acting on the genetic contributions 
#' @param sigmae numeric, environmental standard deviation
#' @param t time from the beginning of the process at which heritability should be calculated, i.e.
#' epidemiologic time
#' @export
calcSigmagAtTime <- function(alpha, sigma, sigmae, t) {
  lenoutput <- max(sapply(list(alpha, sigma, sigmae, t), length))
  if(lenoutput>1) {
    alpha <- rep(alpha, lenoutput/length(alpha))
    sigma <- rep(sigma, lenoutput/length(sigma))
    sigmae <- rep(sigmae, lenoutput/length(sigmae))
    t <- rep(t, lenoutput/length(t))
  }
  ifelse(alpha>0, 0.5*(1-exp(-2*alpha*t))/alpha*sigma^2, t*sigma^2)  
}

#' Expected value of the Ornstein-Uhlenbeck process conditioned on initial value
#' @param g0 numeric, initial value
#' @param alpha the strength of the selection
#' @param theta long term mean value of the OU process
#' @param sigma the unit-time standard deviation of the random component in the OU
#'     process.
#' @param t time at which the expected value should be estimated
#' @export
calcMuAtTime <- function(g0, alpha, theta, sigma, t=Inf) {
  if(t==Inf) {
    theta
  } else {
    theta + (g0 - theta)*exp(-alpha*t)
  }
}


#' Broad-sense heritability H2 for known genetic contributions, g, 
#' and known environmental deviations, e, at the tips of the tree.
#' @export
calcH2emp <- function(tree, g, e) {
  N <- length(tree$tip.label) 
  var(g[1:N])/var(g[1:N] + e[1:N])
}

# with standard double precision this code is numerically unstable for small
#
paramsForkOUR1 <- function(tree, x, e, alpha, theta, sigma, sigmae) {
  if(length(edgesFrom(tree, tree$edge[e, 2]))==0) {
    # not a true fork because the edge e is pointing to a tip. 
    # the probability is a function exp(ig3^2+jg3+k)
    t1 <- tree$edge.length[[e]]
    if(sigmae==0) {
      # no environmental deviation
      g1 <- x[[tree$edge[e, 2]]]  # tip value
      
      c(-(alpha/((-1 + exp(2*t1*alpha))*sigma^2)), 
           (2*(exp(t1*alpha)*(g1 - theta) + theta)*alpha)/
             ((-1 + exp(2*t1*alpha))*sigma^2),
           log(sqrt(alpha/(pi*(1 - exp(-2*t1*alpha))*sigma^2))) - 
             ((exp(t1*alpha)*(g1 - theta) + theta)^2*alpha)/
             ((-1 + exp(2*t1*alpha))*sigma^2)
      )
    } else {
      v1 <- x[[tree$edge[e, 2]]]  # tip value
      # integrate over g1
      # params for P(v1|g1,g3) = exp(mg1^2+ng1+og3^2+pg3+qg1g3+r)
      m <- -(alpha/((1 - exp(-2*alpha*t1))*sigma^2)) - 1/(2*sigmae^2)
      n <- (2*alpha*exp(alpha*t1)*theta)/((1 + exp(alpha*t1))*sigma^2) + 
        v1/sigmae^2
      o <- -alpha/((-1 + exp(2*alpha*t1))*sigma^2)
      p <- -(2*alpha*theta)/((1 + exp(alpha*t1))*sigma^2)
      q <- (2*alpha*exp(alpha*t1))/((-1 + exp(2*alpha*t1))*sigma^2)
      r <- -(alpha*(-1 + exp(alpha*t1))*theta^2)/((1 + exp(alpha*t1))*sigma^2) -
        v1^2/(2*sigmae^2) - log(pi) - log(sigma) - log(sigmae) -
        log(sqrt(2*(1 - exp(-2*alpha*t1))/alpha))
      
      # print(paste('TIP: t1=', t1))
      paramsIntegForkOUR(m, n, o, p, q, r)
    }
  } else {
    esubs <- edgesFrom(tree, tree$edge[e, 2])
    # print('PSUBS:')
    
    abc <- paramsForkOUR1(tree, x, esubs[1], alpha, theta, sigma, sigmae)
    for(es in esubs[-1]) 
      abc <- abc + paramsForkOUR1(tree, x, es, alpha, theta, sigma, sigmae)
    
    
    #abc <- c(paramsForkOUR1(tree, x, esubs[1], alpha, theta, sigma, sigmae),
    #           paramsForkOUR1(tree, x, esubs[2], alpha, theta, sigma, sigmae))
    
    t3 <- tree$edge.length[[e]]
    # params for fCondOU(g3|g4) = exp(mg3^2+ng3+og4^2+pg4+qg3g4+r)
    m <- (-(alpha/((1 - exp(-2*t3*alpha))*sigma^2)))
    n <- ((2*exp(t3*alpha)*alpha*theta)/((1 + exp(t3*alpha))*sigma^2))
    o <- (-(alpha/((-1 + exp(2*t3*alpha))*sigma^2)))
    p <- ((-2*alpha*theta)/((1 + exp(t3*alpha))*sigma^2))
    q <- ((2*exp(t3*alpha)*alpha)/((-1 + exp(2*t3*alpha))*sigma^2))
    r <- log(sqrt(alpha/(pi*(1 - exp(-2*t3*alpha))*sigma^2))) +
      ((1 - exp(t3*alpha))*alpha*theta^2)/((1 + exp(t3*alpha))*sigma^2)

    # print(paste('PFFC: t3=', t3))
    pffc <- paramsForkForkCondOUR(abc[1], abc[2], abc[3], 
                                   m, n, o, p, q, r)
    paramsIntegForkOUR(pffc[1], pffc[2], pffc[3], pffc[4], pffc[5], pffc[6])
  }
}

paramsForkOUR2 <- function(tree, x, e, alpha, theta, sigma, sigmae) {
  sigma2 <- sigma*sigma
  alphasigma2 <- alpha/sigma/sigma
  sigmae2 <- sigmae*sigmae
  
  if(length(edgesFrom(tree, tree$edge[e, 2]))==0) {
    # not a true fork because the edge e is pointing to a tip. 
    # the probability is a function exp(ig3^2+jg3+k)
    t1 <- tree$edge.length[[e]]
    et1alpha <- exp(t1*alpha)
    e2t1alpha <- exp(2*t1*alpha)
    
    if(sigmae==0) {
      # no environmental deviation
      g1 <- x[[tree$edge[e, 2]]]  # tip value
      
      c(alphasigma2 / (1-e2t1alpha), 
        2*(et1alpha*(g1-theta) + theta)/(e2t1alpha/alpha - 1/alpha)/sigma2,
        0.5*log(alphasigma2/(pi*(1 - 1/e2t1alpha))) - 
          (et1alpha*(g1-theta) + theta)^2/(e2t1alpha/alpha - 1/alpha)/sigma2
      )
    } else {
      v1 <- x[[tree$edge[e, 2]]]  # tip value
      # integrate over g1
      # params for P(v1|g1,g3) = exp(mg1^2+ng1+og3^2+pg3+qg1g3+r
      m <- -alphasigma2 - 1/(e2t1alpha/alpha-1/alpha)/sigma2 - 1/(2*sigmae2)
      
      p <- -2*alphasigma2*theta/(1 + et1alpha)
      n <- -p*et1alpha + v1/sigmae2
      
      o <- 1/(1/alpha - e2t1alpha/alpha)/sigma2
      q <- -2*o*et1alpha
      
      r <- -(alphasigma2*(et1alpha - 1)*theta^2)/(et1alpha + 1) -
        v1^2/(2*sigmae2) - log(pi) - log(sigma) - log(sigmae) -
        0.5*log(2*(1 - 1/e2t1alpha)/alpha)
      
     # print(paste('TIP: t1=', t1))
      paramsIntegForkOUR(m, n, o, p, q, r)
    }
  } else {
    esubs <- edgesFrom(tree, tree$edge[e, 2])
    #print('PSUBS:')
    
    abc <- paramsForkOUR2(tree, x, esubs[1], alpha, theta, sigma, sigmae)
    for(es in esubs[-1]) { 
      abc <- abc + paramsForkOUR2(tree, x, es, alpha, theta, sigma, sigmae)
    }
    
    t3 <- tree$edge.length[[e]]
    et3alpha <- exp(t3*alpha)
    e2t3alpha <- et3alpha^2
    
    # params for fCondOU(g3|g4) = exp(mg3^2+ng3+og4^2+pg4+qg3g4+r)
    m <- -alphasigma2 - 1/(e2t3alpha/alpha-1/alpha)/sigma2
    
    p <- (-2*alphasigma2*theta)/(1 + et3alpha)
    n <- -p*et3alpha
    
    o <- 1/(1/alpha - e2t3alpha/alpha)/sigma2
    q <- -2*o*et3alpha
    
    r <- 0.5*log(alphasigma2/(pi*(1 - 1/e2t3alpha))) +
      ((1 - et3alpha)*alpha*theta^2)/((1 + et3alpha)*sigma2)
    
    #print(paste('PFFC: t3=', t3))
    pffc <- paramsForkForkCondOUR(abc[1], abc[2], abc[3],
                                   m, n, o, p, q, r)
    paramsIntegForkOUR(pffc[1], pffc[2], pffc[3], pffc[4], pffc[5], pffc[6])
  }
}

paramsForkForkCondOUR <- function(i, j, k, 
                                   m, n, o, p, q, r) {
  # exp(i1*g3^2+j1*g3+k1)*exp(i2*g3^2+j2*g3+k2)*exp(mg3^2+ng3+og4^2+pg4+qg3g4+r)
  #  = exp(x*g3^2 + y*g3 + o*g4^2 + p*g4 + q*g3g4 + z)
  res <- c(i+m, j+n, o, p, q, k+r)
#   cat('paramsForkForkCondOUR: ')
#   if(class(m)=='mpfr') {
#     print(toNum(round(c(i, j, k, m, n, o, p, q, r, res))))
#   } else {
#     print(round(c(i, j, k, m, n, o, p, q, r, res), digits=4))
#   }
  res
}

paramsIntegForkOUR <- function(a, b, c, d, e, f) {
  # Integrate[exp(ag3^2+bg3+cg4^2+dg4+eg3g4+f), (g3, -Inf, Inf)]
  #    = x*g4^2 + y*g4 + z
  if(any(is.na(a) | a >= 0)) {
    warning(paste("paramsIntegForkOUR: Definite integral is not converging to a finite number, argument a should be negative, but was", a))
    c(NA, NA, NA)    
  } else {
    # a < 0
    res <- c(c-e*e/(4*a), d-b*e/(2*a), log(sqrt(-pi/a))-b*b/(4*a)+f) 
#         cat('paramsIntegForkOUR: ')
#         if(class(a)=='mpfr') {
#           print(toNum(round(c(a, b, c, d, e, f, res))))
#         } else {
#           print(round(c(a, b, c, d, e, f, res), digits=4))
#         }
    res
  }
}
