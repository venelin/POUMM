#' @import ape
NULL

density2 <- function (x) {
  #copied from coda::densplot
  xx <- as.matrix(x)
  y <- xx[, 1, drop = TRUE]
  
  bwf <- function(x) {
    x <- x[!is.na(as.vector(x))]
    return(1.06 * min(sd(x), IQR(x)/1.34) * length(x)^-0.2)
  }
  bw <- bwf(y)
  width <- 4 * bw
  
  scale <- "open"
  if (max(y) <= 1 && 1 - max(y) < 2 * bw) {
    if (min(y) >= 0 && min(y) < 2 * bw) {
      scale <- "proportion"
      y <- c(y, -y, 2 - y)
    }
  }
  else if (min(y) >= 0 && min(y) < 2 * bw) {
    scale <- "positive"
    y <- c(y, -y)
  }
  else scale <- "open"
  dens <- density(y, width = width)
  if (scale == "proportion") {
    dens$y <- 3 * dens$y[dens$x >= 0 & dens$x <= 
                           1]
    dens$x <- dens$x[dens$x >= 0 & dens$x <= 1]
  }
  else if (scale == "positive") {
    dens$y <- 2 * dens$y[dens$x >= 0]
    dens$x <- dens$x[dens$x >= 0]
  }
  
  dens
}

#' @export
analyseMCMCs <- function(chains, stat=function(x, ...) {x}, statName=NULL, 
                         start, end, thin, ...) {
  mcmcs <- lapply(chains, function(obj) {
    if(is.null(obj$mcmc)) 
      convertToMCMC(obj, thin)
    else
      obj
  })
  
  mcs <- as.mcmc.list(lapply(mcmcs, function(mc) {
    winmcmc <- window(mc$mcmc, start=start, end=end, thin=thin)
    data <- matrix(apply(winmcmc, 1, stat, ...), nrow=nrow(winmcmc), byrow=T)
    mcmc(data, start=start(winmcmc), end=end(winmcmc), thin=thin(winmcmc))
  }))
  names(mcs) <- names(chains)
  mcs.mat <- as.matrix(mcs)
  
  log.ps <- window(as.mcmc.list(lapply(mcmcs, function(mc) { mc$log.p })),
                   start=start, end=end, thin=thin)
  
  if(length(chains) > 1) {
    gelman <- gelman.diag(mcs, autoburnin=F)
    gelman <- ifelse(ncol(mcs.mat) > 1, gelman$mpsrf, gelman$psrf[2]) 
  } else {
    gelman <- NULL
  }
  ESS <- effectiveSize(mcs)
  names(ESS) <- NULL
  
  dens <- density2(mcs.mat)
  HPD <- HPDinterval(mcs)
  
  l <- list(ESS=ESS, Minx=min(dens$x), Maxx=max(dens$x), 
            MaxDens=max(dens$y), Mode=dens$x[which.max(dens$y)],
            HPD=HPD, G.R.=gelman, mcs=mcs, log.ps=log.ps)
  
  # without this we are running into trouble:
  names(l) <- c('ESS', 'Minx', 'Maxx', 'MaxDens', 'Mode', 'HPD', 'G.R.', 'mcs', 'log.ps')    
  l
}

#' MCMC-sampling from a posterior distribution of a PMM(OU) model given data and a prior distribution
#'
#' @param v Either a numeric vector containing the phenotypic values at the tips of tree or
#' a named list containing named elements v - a numeric vector and tree - a phylo object. 
#' @param tree a phylo object
#' @param distgr character vector, see parameter distgr of the likVTreeOU function
#' @param divideEdgesBy numeric, by which all branch lengths in the tree 
#' are to be divided prior to likelihood calculation. Only the branch lengths in the tree are changed 
#' and appropriate rescaling of the OU parameters alpha and sigma is not made. For example, 
#' if this parameter is set to 100, the resulting parameter alpha should be divided by 100 and the parameter 
#' sigma should be divided by 10. 
#' @param parInit
#' @param parPrior
#' @param parToATSSe
#' @param scale
#' @param n.mcmc
#' @param n.adapt
#' @param thin
#' @param acc.rate
#' @param n.chains
#' @param samplePrior
#' @param analysis
#' 
#' @details Currently, this function calls the MCMC function from the adaptMCMC package.
#' 
#' @return a list of coda objects
#' 
#' @export
mcmcOUTree <- function(v, tree, distgr=c('maxlik', 'normal'), divideEdgesBy=1, 
                       parToATSSe=function(par) {c(par[1], par[2], sqrt(par[3:4]))},
                       parInit=function(chainNo) {c(alpha=0, theta=0, sigma2=1, sigmae2=1)},
                       parPrior=function(par) {dexp(par[1], rate=.01, T)+dunif(par[2], 0, 100, T)+
                                                 dexp(par[3],  rate=.0001, T)+dexp(par[4], rate=.01, T)},
                       scale=matrix(c(100,  0.00, 0.00, 0.00,
                                      0.00, 10.0, 0.00, 0.00,
                                      0.00, 0.00, 100,  0.00,
                                      0.00, 0.00, 0.00, 100), nrow=4, ncol=4, byrow=T),
                       n.mcmc=1.2e6, n.adapt=2e5, thin=100, acc.rate=0.01, gamma=0.5, n.chains=1, samplePrior=F, ..., 
                       analysis=list(stat=function(x) {x}, statName=NULL, thin=100)) {
  if(!require(adaptMCMC)) {
    stop('Please, install the adaptMCMC package.')
  }
  if(!require(coda)) {
    stop('Please, install the coda package.')
  }
  convertToMCMC <- function(obj, thin=1) {
    codaobj <- convert.to.coda(obj)
    winmcmc <- window(codaobj, start=start(codaobj), end=end(codaobj), thin=thin)
    log.p <- window(mcmc(matrix(obj$log.p, ncol=1), start=start(codaobj), 
                         end=end(codaobj), thin=1), 
                    start=start(codaobj), end=end(codaobj), thin=thin)
    list(mcmc=winmcmc, log.p=log.p, acc.rate=obj$acceptance.rate, 
         adaption=obj$adaption, n.sample=obj$n.sample, cov.jump=obj$cov.jump,
         sampling.parameters=obj$sampling.parameters, 
         scale.start=obj$scale.start)
  }
  
  if(is.list(v)) {
    p <- v
    v <- p$v
    tree <- p$tree
  }
  if(divideEdgesBy!=1) {
    tree$edge.length <- tree$edge.length/divideEdgesBy
  }
  
  post <- function(par, ...) {
    pr <- if(any(par[c(1, 3, 4)]<0)) {
      -Inf
    } else {
      parPrior(par)
    }
    if(is.nan(pr)|is.na(pr)) {
      -Inf
    } else if(is.infinite(pr)) {
      pr
    } else {
      if(samplePrior) {
        lik <- 0
      } else {
        atsse <- parToATSSe(par)
        names(atsse) <- NULL
        lik <- likVTreeOU(v, tree, alpha=atsse[1], theta=atsse[2], sigma=atsse[3], sigmae=atsse[4], 
                          distgr=distgr[1], log=T, ...)
      }
      pr+lik  
    }
  }
  
  chains <- list()
  for(i in 1:n.chains) {
    init <- parInit(i)
    print(init)
    ch <- convertToMCMC(MCMC(post, n=n.mcmc, init=init, scale=scale, adapt=n.adapt, 
                             acc.rate=acc.rate, gamma=gamma, ...), thin=thin)
    ch$scale.start <- scale      
    chains[[i]] <- ch
  }
  
  an <- analyseMCMCs(chains, stat=analysis$stat, statName=analysis$statName, 
                     start=round(n.mcmc/5), end=n.mcmc, thin=ifelse(!is.null(analysis$thin), analysis$thin, thin))
  list(chains=chains, analysis=an)
}
                       
#' Find a maximum likelihood fit of the PMM(OU) model
#'
#' @param v Either a numeric vector containing the phenotypic values at the tips of tree or
#' a named list containing named elements v - a numeric vector and tree - a phylo object. 
#' @param tree a phylo object
#' @param distgr character vector, see parameter distgr of the likVTreeOU function
#' @param parFixed a named numeric vector, indicating fixed values for some of the parameters
#' alpha, theta, sigma and sigmae, by default, c().
#' @param parMin, parMax two named numeric vectors indicating the boundaries of the search region. 
#' Default values are parMin=c(alpha=0, theta=0, sigma=0, sigmae=0) and parMax=c(alpha=100, 
#' theta=10, sigma=20, sigmae=10).
#' @param divideEdgesBy numeric, by which all branch lengths in the tree 
#' are to be divided prior to likelihood calculation. Only the branch lengths in the tree are changed 
#' and appropriate rescaling of the OU parameters alpha and sigma is not made. For example, 
#' if this parameter is set to 100, the resulting parameter alpha should be divided by 100 and the parameter 
#' sigma should be divided by 10. 
#' @param control list of parameters passed on to optim
#' @param ... additional parameters passed to the likVTreeOU function
#' 
#' @return a list containing an element par and an element value as well as the parameters passed to the 
#' call of mlOUTree
#' 
#' 
#' @export
mlOUTree <- function(v, tree, distgr=c('maxlik', 'normal'), parFixed=c(), 
                     parMin=c(alpha=0, theta=0, sigma=0, sigmae=0), 
                     parMax=c(alpha=100, theta=10, sigma=20, sigmae=10),
                     divideEdgesBy=1, control=list(), ...) {
  if(is.list(v)) {
    p <- v
    v <- p$v
    tree <- p$tree
  }
  if(divideEdgesBy!=1) {
    tree$edge.length <- tree$edge.length/divideEdgesBy
  }
  
  if(is.null(names(parMin))) {
    names(parMin) <- c('alpha', 'theta', 'sigma', 'sigmae')
  }
  if(is.null(names(parMax))) {
    names(parMax) <- c('alpha', 'theta', 'sigma', 'sigmae')
  }
  
  parMin <- parMin[setdiff(names(parMin), names(parFixed))]
  parMax <- parMax[setdiff(names(parMax), names(parFixed))]
  parFixedNoAlpha <- parFixed[setdiff(names(parFixed), 'alpha')]
  
  memoriseMin <- function(f, ...) {
    memo <- new.env()
    assign('val', Inf, pos=memo)
    assign('par', NULL, pos=memo)
    assign('count', 0, pos=memo)
    
    function(..., report=F) {
      count <- get('count', memo)
      valMemo <- get('val', pos=memo)
      parMemo <- get('par', pos=memo)
      log <- F
      if(report) {
        list(count=count, value=valMemo, g0=attr(valMemo, 'grmax'), par=parMemo)
      } else {
        count <- count+1
        assign('count', count, pos=memo)
        val <- do.call(f, list(...))
        params <- unlist(as.list(unlist(list(...)))[c('alpha', 'theta', 'sigma', 'sigmae')])
        if(val < valMemo) {
          assign('par', params, pos=memo)
          assign('val', val, pos=memo)
          valMemo <- val
          parMemo <- params
          log <- T
        } 
#         if(count%%100==0|log) {
#           cat('Count calls=', count, '; parMemo=', toString(round(parMemo, 2)), 
#               '; valMemo=', valMemo, ';   ### par=', toString(round(params, 2)), 
#               '; value=', val, '\n', sep='')
#         }
        val
      }
    }
  }
  
  memoriseAll <- function(f, ...) {
    memo <- new.env()
    assign('results', NULL, pos=memo)
    
    function(..., report=F, valueOnly=T, hideLastResultOnNextCall=F) {
      if(hideLastResultOnNextCall) {
        # for the next call, don't pass the last result
        assign('hideLastResult', T, pos=memo)
      } else {
        results <- get('results', pos=memo)
        if(report) {
          list(results=results)
        } else {
          hideLastResult <- get('hideLastResult', pos=memo)
          if(hideLastResult) {
            res <- do.call(f, list(...))
            # enable passing last result on the next call
            assign('hideLastResult', F, pos=memo)
          } else {
            res <- do.call(f, list(..., lastResult=results[[length(results)]]))  
          }
          results <- c(results, list(res))
          assign('results', results, pos=memo)
          
          if(valueOnly) {
            res$value
          } else {
            res
          }
        }  
      }
    }
  }
  
  minusll <- memoriseMin(function(par, ...) {
    value <- -do.call(likVTreeOU, c(list(v, tree), as.list(parFixedNoAlpha), as.list(par), list(log=T, distgr=distgr), list(...)))
    if(is.na(value) | is.nan(value) | is.infinite(value)) {
      1e20
    } else {
      value
    }
  })
  
  optimFixedAlpha <- memoriseAll(function(alpha, parFixed, parMin, parMax, lastResult=NULL, useLastResult=T, ...) {
    
    if(useLastResult & is.list(lastResult) & all(lastResult$par>=parMin) & all(lastResult$par<=parMax)) {
      parInit <- lastResult$par
    } else {
      lenPar <- 4-1-length(parFixed)
      if(lenPar>0) {
        parInit <- sapply(1:(4-1-length(parFixed)), function(argi) {
          runif(1, min=parMin[argi], max=parMax[argi])
        })
        names(parInit) <- names(parMin)  
      } else {
        parInit <- NULL
      }
    }
    if(!is.null(parInit)) {
      cat('optimFixedAlpha called on alpha=', alpha, ', parFixed=', toString(parFixed), ', parInit=', parInit, ', ...=', toString(list(...)), '\n')
      res <- optim(par=parInit, fn=minusll, gr=NULL, alpha=alpha, ..., method='L-BFGS-B', lower=parMin, upper=parMax, control=control)
      res$alpha <- alpha
      res  
    } else {
      res <- list()
      res$value <- minusll(c(alpha=alpha))
      res$alpha <- alpha
      res$par <- c()
      res
    }
})
  
  if('alpha' %in% names(parFixed)) {
    alphaFixed=parFixed['alpha']
    # this is a bug : hiding last result only
    optimFixedAlpha(alpha=alphaFixed, parFixed=parFixedNoAlpha, parMin=parMin, parMax=parMax, valueOnly=F, useLastResult=F, hideLastResultOnNextCall=T, ...)
    # as a workaround we call it a second time
    optimFixedAlpha(alpha=alphaFixed, parFixed=parFixedNoAlpha, parMin=parMin, parMax=parMax, valueOnly=F, useLastResult=F, hideLastResultOnNextCall=F, ...)
  } else {
    lowerAlpha <- 2^c(floor(log2(parMin['alpha'])), seq(max(1, ceiling(log2(parMin['alpha']))), ceiling(log2(parMax['alpha']))+1, by=1))
    if(length(lowerAlpha)>2) {
      upperAlpha <- lowerAlpha[-(1:2)]
      lowerAlpha <- lowerAlpha[1:(length(lowerAlpha)-2)]
      upperAlpha[upperAlpha>parMax['alpha']] <- parMax['alpha']
      lowerAlpha[lowerAlpha<parMin['alpha']] <- parMin['alpha']
    } else {
      lowerAlpha <- parMin['alpha']
      upperAlpha <- parMax['alpha']
    }
    
    parMinNoAlpha <- parMin[setdiff(names(parMin), 'alpha')]
    parMaxNoAlpha <- parMax[setdiff(names(parMax), 'alpha')]
    res <- list()
    for(i in 1:length(lowerAlpha)) {
      optimFixedAlpha(hideLastResultOnNextCall=T)
      optimise(optimFixedAlpha, interval=c(lowerAlpha[i], upperAlpha[i]), parFixed=parFixedNoAlpha, 
               parMin=parMinNoAlpha, parMax=parMaxNoAlpha, valueOnly=T, ...)
    } 
  }
  minRep <- minusll(report=T)
  minRep$par <- c(minRep$par, parFixed)[c('alpha', 'theta', 'sigma', 'sigmae')]
  allRep <- optimFixedAlpha(report=T)
  c(minRep, allRep, list(distgr=distgr, parFixed=parFixed, 
                         parMin=parMin, 
                         parMax=parMax, divideEdgesBy=divideEdgesBy, control=control))
}
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
  ttips <- nodeTimes(tree)[1:N]
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
#' be used. One of  'R1', 'R2', 'RCPP1'. Defaults to 'R2'. 
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
                       log=T, impl=c('R4', 'R2', 'R1'), 
                       usempfr=1, maxmpfr=2, precbits=512, debug=F) {
  alphaorig <- alpha
  thetaorig <- theta
  sigmaorig <- sigma
  sigmaeorig <- sigmae
  mugrorig <- mugr
  sigmagrorig <- sigmagr
  vorig <- v
  
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
          
          if(edge[es[1], 2] <= N) {
            # all es pointing to tips
            if(sigmae==0) {
              # no environmental deviation
              g1 <- v[edge[es, 2]]  # tip values
              etalphag1thetatheta <- etalpha[es]*(g1-theta)+theta
              
              pif[edge[es, 2], ]  <- 
                c(fe2talphasigma2[es],
                  -2 * etalphag1thetatheta * fe2talphasigma2[es],
                  r0[es] + etalphag1thetatheta^2 * fe2talphasigma2[es]
                )
            } else {
              v1 <- v[edge[es, 2]]  # tip values
              
              u <- -0.5/sigmae2
              w <- v1/sigmae2
              z <- -0.5*loge2 - 0.5*logpi  - v1^2/(2*sigmae2) - logsigmae 
              
              gutalphasigma2 <- (fe2talpha[es]-alpha+u*sigma2)/fe2talpha[es]
              
              pif[edge[es, 2], ] <- c(u/gutalphasigma2,
                                      (etalpha[es]*(2*theta*u+w)-2*theta*u)/gutalphasigma2,
                                      -0.5*log(gutalphasigma2) -
                                        0.25*sigma2*w^2/(fe2talpha[es]-alpha + sigma2*u) + 
                                        alpha*theta*(theta*u-etalpha[es]*(theta*u+w))/((1+etalpha[es])*(sigma2*u-alpha)+fetalpha[es]) + 
                                        talpha[es]+z)
            }
          } else {
            # edges pointing to internal nodes, for which all children nodes have been visited
            u <- pif[edge[es, 2], 1]
            w <- pif[edge[es, 2], 2]
            z <- pif[edge[es, 2], 3]
            
            gutalphasigma2 <- (fe2talpha[es]-alpha+u*sigma2)/fe2talpha[es]
            
            pif[edge[es, 2], ] <- c(u/gutalphasigma2,
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

#' Construct a (log-)posterior function for parameters given a tree and phenotype values at the tips, assuming the PMM/OU model.
#' @param p a list containing an element called tree of class phylo and a numeric vector called v of length the number of tips in the 
#' tree
#' @param prior a function object calculating the prior (log-) density for given parameters 
#' @param parToATSSeL a function object converting the parameters par to a vector of at least 5 values, the first five of which should be
#' corresponding to alpha, theta, sigma, sigmae, and lambda, to be passed to likVTreeOU
#' @param addParam a function of p, returning a list of additional parameters to be passed to the prior and 
#' parToATSSeL functions, defaults to a funciton returning a list of a single element which is the input list p.
#' @param controlLikVTreeOU list, additional parameters to pass to likVTreeOU apart from tree, v, alpha, sigma, sigmae, lambda and log
#' @param log logical whether log-posterior should be calculated
#' @param reportFreq numeric between 0 and 1 indicating how often report messages should be printed to the console, default 
#' is 0 meaning that no output to the console should be generated
#' @return a function of one argument par
#' 
#' @note The function addParam is called only once during the posterior creation, and its result is passed to the 
#' prior and parToATSSeL function every time.
#' @export
makePostVTreeOU <- function(p, 
                            prior=function(par, ...) {stop('You should define a prior function!')}, 
                            parToATSSeL=function(par, ...) {par}, 
                            addParam=function(p) {list(p=p)},
                            controlLikVTreeOU=NULL, log=T, reportFreq=0, ...) {
  additional <- addParam(p)
  function(par) {
    pr <- do.call(prior, c(list(par=par), additional))
    
    if(is.infinite(pr)) {
      if(pr>0) {
        cat('Infinite prior for par=', toString(par), '\n', sep='')
      }
      pr
    } else {
      atssel <- do.call(parToATSSeL, c(list(par=par), additional)) 
      lik <- do.call(likVTreeOU, 
                     c(list(v=p$v, tree=p$tree, alpha=atssel[1], theta=atssel[2], sigma=atssel[3], sigmae=atssel[4], 
                            lambda=atssel[5], log=log), controlLikVTreeOU))
      
      post <- if(log) lik+pr else lik*pr
      
      if(runif(1)<reportFreq) {
        print(c(atssel, c(lik=lik, pr=pr, post=post)))
      }
      post
    }
  }
}

#' Node indices of the direct descendants of n in the phylogeny tree.
#' 
#' @param tree an object of class phylo
#' @param n an index of a node (root, internal or tip) in tree
#' @return
#'   An integer vector.
#'   
#' @export
chld <- function(tree, n) {
  as.integer(tree$edge[tree$edge[, 1]==n, 2])
}

#' Edge indices of the edges in tree starting from n
#' @param tree an object of class phylo
#' @param n an index of a node (root, internal or tip) in tree
#' @return
#'   An integer vector.
#' @export
edgesFrom <- function(tree, n) {
  which(tree$edge[, 1]==n)
}

#' Calculate the time from the root to each node of the tree
#' 
#' @param tree an object of class phylo
#' @return
#'   a vector of size the number of nodes in the tree (tips, root, internal)
#'   containing the time from the root to the corresponding node in the tree
#' @export
nodeTimes <- function(tree) {
  rtree <- reorder(tree, 'postorder')
  es <- rtree$edge[dim(rtree$edge)[1]:1, ]
  nEdges <- dim(es)[1]
  ts <- rev(rtree$edge.length)
  nodeTimes <- rep(0, length(rtree$tip.label)+rtree$Nnode)
  for(e in 1:nEdges) 
    nodeTimes[es[e, 2]] <- nodeTimes[es[e, 1]]+ts[e]
  nodeTimes
}

#' Make the ending edges of a tree longer/shorter by addLength
#' @param tree a phylo object
#' @param addLength numeric vector of length 1 or length(tree$tip.label)
#' @return the tree with modified edge lengths for the edges leading to tips
#' @export
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
#' @export
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
  
  while(T) {
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

#' Random generation of a trait along a tree following an OU process 
#' @param tree an object of class phylo
#' @param g0 initial trait value (at the root of tree)
#' @param alpha, theta, sigma parameters of the OU process.
#' @return a vector containing the trait values at all nodes of the tree
#' 
#' @export
generateTraitOU <- function(tree, g0, alpha, theta, sigma) {
  rTraitCont(tree, 
             function(x, l, .a, .t, .s) {
               OU(x, l, .a, .t, .s)
             }, root.value=g0, ancestor=T, .a=alpha, .t=theta, .s=sigma)
}



#' Random generation of a trait along a tree following a BM process 
#' @param tree an object of class phylo
#' @param g0 initial trait value (at the root of tree)
#' @param sigma parameter of the BM process.
#' @return a vector containing the trait values at all nodes of the tree
#' 
#' @export
generateTraitBM <- function(tree, g0, sigma) {
  rTraitCont(tree, 'BM', root.value=g0, ancestor=T, sigma=sigma)
}

#' Infection couples given a tree and ids of all tips and nodes
#' 
#' @param tree a phylo object
#' @param id a vector containing identifiers of all tips and nodes
#' @return a matrix of the class of id of two columns and N-1 rows
#' @export
infectionCouples <- function(tree, id) {
  edge <- tree$edge
  edge <- matrix(id[edge], ncol=2)
  edge[edge[, 1]!=edge[, 2],]
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
#'   
#' @export
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

#' get a sample list of possible transmission pair assignments
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
#'      eachAfterGettingInfected - the source phenotype is measured at the moment of infection from its own source,
#'        the recipient phenotype is measured at the moment it got infected from the source;
#'      bothAtInfection - the source and the recipient phenotypes are taken at 
#'        the moment of infection from source to recipient.
#'      atTips - the source and the recipient phenotypes are measured at their corresponding tips in the tree.
#' @return a list of nSamples matrices, each matrix having tree columns representing 
#' source phenotype, recipient phenotype and time of the infection from the root
#'   
#' @export
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
#' 
#' @return a numerical vector of length nSamples
#'    
#' @export
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
