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
#' 
#' 
#' @export
mcmcOUTree <- function(v, tree, distgr=c('maxlik', 'normal'), divideEdgesBy=1, 
                       parToATSSe=function(par) {c(par[1], par[2], sqrt(par[3:4]))},
                       parInit=function(chainNo) {c(alpha=0, theta=0, sigma2=1, sigmae2=1)},
                       parPrior=function(par) {
                         if(any(par[c(1, 3, 4)]<0)) {
                           -Inf
                         } else {
                           dexp(par[1], rate=.01, T)+dunif(par[2], 0, 100, T)+
                             dexp(par[3],  rate=.0001, T)+dexp(par[4], rate=.01, T)  
                         }
                       },
                       scale=matrix(c(100,  0.00, 0.00, 0.00,
                                      0.00, 10.0, 0.00, 0.00,
                                      0.00, 0.00, 100,  0.00,
                                      0.00, 0.00, 0.00, 100), nrow=4, ncol=4, byrow=T),
                       n.mcmc=1.2e6, n.adapt=2e5, thin=100, acc.rate=0.01, gamma=0.5, n.chains=1, samplePrior=F, ..., 
                       analysis=list(stat=function(x) {x}, statName=NULL, thin=100)) {
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
  
  pruneInfo=pruneTree(tree)
  
  post <- function(par, ...) {
    pr <- parPrior(par)
    
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
                          distgr=distgr[1], log=T, impl='R5', pruneInfo=pruneInfo, ...)
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
  
  an <- try(analyseMCMCs(chains, stat=analysis$stat, statName=analysis$statName, 
                         start=round(n.mcmc/5), end=n.mcmc, thin=ifelse(!is.null(analysis$thin), analysis$thin, thin)),
            silent=TRUE)
  list(chains=chains, analysis=an)
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

