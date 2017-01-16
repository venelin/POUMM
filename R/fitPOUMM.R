memoiseMax <- function(f, par, memo, verbose, ...) {
  countMemo <- mget('count', envir = memo, ifnotfound = list(0))$count
  valMemo <- mget('val', envir = memo, ifnotfound = list(-Inf))$val
  parMemo <- mget('par', envir = memo, ifnotfound = list(NULL))$par
  
  assign("count", countMemo + 1, pos = memo)

  val <- f(par, memo = memo, ...)

  if(valMemo < val) {
    assign('par', par, pos = memo)
    assign('val', val, pos = memo)
    valMemo <- val
    parMemo <- par
    if(verbose) {
      cat('Call ', countMemo, ': loglik on par=(', 
          toString(round(par, 6)), "): ", val, "\n", sep = "")
    }
  } 
  val
}

#' Find a maximum likelihood fit of the POUMM model
#'
#' @param loglik function(par, memo, parFixedNoAlpha) 
#' @param parFixed A named numeric vector, indicating fixed values for some of 
#'   the parameters alpha, theta, sigma and sigmae, by default, c().
#' @param parLower,parUpper Two named numeric vectors indicating the boundaries of 
#'   the search region. Default values are parLower = c(alpha = 0, theta = 0, 
#'   sigma = 0, sigmae = 0) and parUpper = c(alpha = 100, theta = 10, sigma = 20, 
#'   sigmae = 10).
#' @param tol Numeric accuracy parameter passed to optimize, defaults to 0.001. 
#'   (see Details)
#' @param control List of parameters passed on to optim, default 
#'   list(factr = 1e8, fnscale = -1), see ?optim. 
#' @param verbose A logical indicating whether to print informative messages on 
#'   the standard output.
#'   
#' @details If alpha is not fixed, then the optimization is done via calls to 
#'   optimize over alpha, calling optim on the remaining variable parameters for
#'   each alpha. The parameter tol is passed to optimize (default value is 
#'   0.001), which guarantees that no evaluation is done on points closer than 
#'   eps*|alpha*|+tol/3, where eps is .Machine$double.eps and alpha* is the 
#'   optimal alpha.
#' @return a list containing an element par and an element value as well as the 
#'   parameters passed
#'   
#' @importFrom stats optim optimise runif
maxLikPOUMMGivenTreeVTips <- function(
  loglik, parLower, parUpper, parInitML = NULL, g0 = NA, g0Prior = NULL,
  control=list(factr = 1e8, fnscale = -1), 
  verbose = FALSE, ...) {
  if(is.null(parLower) | is.null(parUpper)) {
    stop("parLower and/or parUpper are NULL.")
  }
  
  
  if(is.null(parInitML)) {
    listParInitML <- list(0.1 * (parLower + parUpper), 
                      0.5 * (parLower + parUpper),
                      0.9 * (parLower + parUpper))
  } else {
    listParInitML <- list(parInitML)
  }
  
  
  memoMaxLoglik <- new.env()
  fnForOptim <- function(par) {
    ll <- memoiseMax(
      loglik, par = c(par, g0 = g0), memo = memoMaxLoglik, verbose = verbose)
    g0LogPrior <- attr(ll, "g0LogPrior")
    
    if(is.finite(g0LogPrior)) {
      ll + g0LogPrior
    } else {
      ll
    }
  }
  
  for(iOptimTry in 1:length(listParInitML)) {
    parInitML <- listParInitML[[iOptimTry]]
    if(!(all(parInitML >= parLower) & all(parInitML <= parUpper))) {
      stop(paste0("All parameters in parInitML should be between parLower=", 
                  toString(round(parLower, 6)), "and parUpper=", 
                  toString(round(parUpper, 6)), ", but were ", 
                  toString(round(parInitML, 6)), "."))
    }
    
    if(is.null(control)) {
      control <- list()
    }
    
    # ensure that optim does maximization.
    control$fnscale <- -1
    
    if(length(parInitML) > 0) {
      res <- optim(fn = fnForOptim, 
                   par = parInitML, lower = parLower, upper = parUpper, 
                   method = 'L-BFGS-B', control = control)
    } else {
      stop("ML-fit: parInitML has zero length.")
    }
  }
  
  maxPar <- get("par", pos = memoMaxLoglik)
  maxValue <- get("val", pos = memoMaxLoglik)
  callCount <- get("count", pos = memoMaxLoglik)
  
  list(par = maxPar,
       value = maxValue, 
       g0 = attr(maxValue, which="g0", exact = TRUE), 
       g0LogPrior = attr(maxValue, which="g0LogPrior", exact = TRUE), 
       count = callCount, parLower=parLower, parUpper=parUpper, 
       control=control)
}

# Old likelihood optimizer, which used optimize on alpha with optim on the remaining parameters. 
ml.poumm <- function(
  loglik, parFixed=c(), 
  parLower=c(alpha=0, theta=0, sigma=0, sigmae=0), 
  parUpper=c(alpha=100, theta=10, sigma=20, sigmae=10),
  parInit = NULL,
  g0 = NA, g0Prior = NULL, 
  tol=0.001, control=list(factr = 1e8, fnscale = -1),  
  verbose=FALSE) {
  
  if(is.null(parLower) | is.null(parUpper)) {
    stop("parLower and/or parUpper are NULL.")
  }
  if(is.null(names(parLower))) {
    names(parLower) <- c('alpha', 'theta', 'sigma', 'sigmae')
  }
  if(is.null(names(parUpper))) {
    names(parUpper) <- c('alpha', 'theta', 'sigma', 'sigmae')
  }

  
  if(is.null(parInit)) {
    parInit <- (parLower + parUpper) / 2
  }
  if(is.null(names(parInit))) {
    names(parInit) <- c('alpha', 'theta', 'sigma', 'sigmae')
  }
  
  if(!(all(parInit >= parLower) & all(parInit <= parUpper))) {
    stop("All parameters in parInit should be between parLower and parUpper.")
  }
  
  if(is.null(control)) {
    control <- list()
  }
  
  # ensure that optim does maximization.
  control$fnscale <- -1
  
  parLower <- parLower[setdiff(names(parLower), names(parFixed))]
  parUpper <- parUpper[setdiff(names(parUpper), names(parFixed))]
  parInit <- parInit[setdiff(names(parInit), names(parFixed))]
  
  parFixedNoAlpha <- parFixed[setdiff(names(parFixed), 'alpha')]
  
  memoMaxLoglik <- new.env()
  fnForOptimOnFixedAlpha <- function(par, alpha) {
    names(alpha) <- "alpha"
    parFull <- 
      c(alpha, par, parFixedNoAlpha)[c("alpha", "theta", "sigma", "sigmae")]
    
    ll <- memoiseMax(
      loglik, par = c(parFull, g0 = g0), memo = memoMaxLoglik, verbose = verbose)
    g0LogPrior <- attr(ll, "g0LogPrior")
    
    if(is.finite(g0LogPrior)) {
      ll + g0LogPrior
    } else {
      ll
    }
  }
  
  # This is going to be called by optimise or by this code if alpha is fixed by 
  # the user; this calls optim for fnForOptimOnFixedAlpha defined above.
  fnForOptimiseOverAlpha <- function(alpha, parFixed, parLower, parUpper, 
                                     memoParInit) {
    parInit <- get("parInit", pos = memoParInit)
    
    if(length(parInit) > 0) {
      res <- optim(par = parInit, 
                   fn = fnForOptimOnFixedAlpha, alpha = alpha, 
                   method = 'L-BFGS-B', lower = parLower, upper = parUpper, 
                   control = control)
      if(res$value == get('val', pos = memoMaxLoglik)) {
        if(verbose) {
          cat("++++ Setting parInit to:", toString(res$par), "\n")
        } 
        assign("parInit", res$par, pos = memoParInit)
      } 
    } else {
      # all parameters except alpha are fixed, so we only evaluate loglik at
      # alpha.
      res <- list(
        value = fnForOptimOnFixedAlpha(par = c(), alpha = alpha)
      )
      res
    }
    res$value
  }
  
  memoParInit <- new.env()
  
  if('alpha' %in% names(parFixed)) {
    # alpha is among the fixed parameters, so no need to call optimise
    alphaFixed=parFixed['alpha']
    
    assign("parInit", parInit, pos = memoParInit)
    
    # the optimum will be stored in memoMaxLoglik
    fnForOptimiseOverAlpha(alpha = alphaFixed, parFixed = parFixedNoAlpha, 
                           parLower = parLower, parUpper = parUpper, 
                           memoParInit = memoParInit)
  } else {
    lowerAlpha <- parLower['alpha']
    upperAlpha <- parUpper['alpha']
    
    parLowerNoAlpha <- parLower[setdiff(names(parLower), 'alpha')]
    parUpperNoAlpha <- parUpper[setdiff(names(parUpper), 'alpha')]
    parInitNoAlpha <- parInit[setdiff(names(parInit), 'alpha')]
    
    assign("parInit", parInitNoAlpha, pos = memoParInit)
    
    for(i in 1:length(lowerAlpha)) {
      optimise(fnForOptimiseOverAlpha, 
               interval=c(lowerAlpha[i], upperAlpha[i]), 
               tol=tol, 
               maximum = TRUE,
               parFixed=parFixedNoAlpha, 
               parLower=parLowerNoAlpha, parUpper=parUpperNoAlpha, 
               memoParInit = memoParInit)
    }
  }
  maxPar <- get("par", pos = memoMaxLoglik)
  maxValue <- get("val", pos = memoMaxLoglik)
  callCount <- get("count", pos = memoMaxLoglik)
  
  list(par = c(maxPar, parFixed)[c('alpha', 'theta', 'sigma', 'sigmae')],
       value = maxValue, 
       g0 = attr(maxValue, which="g0", exact = TRUE), 
       g0LogPrior = attr(maxValue, which="g0LogPrior", exact = TRUE), 
       count = callCount, 
       parFixed=parFixed, parLower=parLower, parUpper=parUpper, 
       control=control)
}


#' @importFrom stats window start end
convertToMCMC <- function(obj, thinMCMC=1) {
  codaobj <- adaptMCMC::convert.to.coda(obj)
  
  winmcmc <- window(codaobj, start = start(codaobj), end = end(codaobj), 
                    thin = thinMCMC)
  
  log.p <- window(mcmc(obj$log.p, 
                       start = start(codaobj), 
                       end = end(codaobj), 
                       thin = 1), 
                  start = start(codaobj), end = end(codaobj), thin = thinMCMC)
  
  list(
    mcmc = winmcmc, log.p = log.p, 
    accRateMCMC = obj$acceptance.rate, 
    adaption = obj$adaption, n.sample = obj$n.sample, cov.jump = obj$cov.jump,
    sampling.parameters = obj$sampling.parameters, 
    scale.start = obj$scale.start)
}

analyseMCMCs <- function(chains, stat=NULL, statName="logpost", 
                         start, end, thinMCMC, as.dt=FALSE, ...) {
  mcmcs <- lapply(chains, function(obj) {
    if(is.null(obj$mcmc)) 
      convertToMCMC(obj, thinMCMC)
    else
      obj
  })
  
  log.ps <- window(coda::as.mcmc.list(lapply(mcmcs, function(mc) { mc$log.p })),
                   start=start, end=end, thin=thinMCMC)
  
  if(statName=='logpost') {
    mcs <- coda::as.mcmc.list(
      lapply(1:length(log.ps), function(i) log.ps[[i]][, 1, drop = FALSE]))
    names(mcs) <- names(chains)
  } else if(statName=='loglik') {
    mcs <- coda::as.mcmc.list(
      lapply(1:length(log.ps), function(i) log.ps[[i]][, 2, drop = FALSE]))
    names(mcs) <- names(chains)
  } else if(statName=='g0') {
    mcs <- coda::as.mcmc.list(
      lapply(1:length(log.ps), function(i) log.ps[[i]][, 3, drop = FALSE]))
    names(mcs) <- names(chains)
  } else {
    mcs <- coda::as.mcmc.list(lapply(mcmcs, function(mc) {
      winmcmc <- window(mc$mcmc, start = start, end = end, thin = thinMCMC)
      data <- matrix(stat(winmcmc, ...), nrow = nrow(winmcmc), byrow = TRUE)
      coda::mcmc(data, start = start(winmcmc), end = end(winmcmc), 
                 thin = thin(winmcmc))
    }))
    names(mcs) <- names(chains)
  }
  
  mcs.mat <- as.matrix(mcs)
  
  if(length(chains) > 1) {
    gelman <- coda::gelman.diag(mcs, autoburnin=FALSE)
    gelman <- ifelse(ncol(mcs.mat) > 1, gelman$mpsrf, gelman$psrf[2]) 
  } else {
    gelman <- NULL
  }
  ESS <- try(lapply(mcs, coda::effectiveSize), silent=TRUE)
  if(class(ESS)=='try-error') {
    warning(paste("For ", statName, 
                  "calling coda::effectiveSize resulted in an error:", 
                  toString(ESS)))
    ESS <- 0
  }
    
  names(ESS) <- NULL
  
  
  HPD <- lapply(mcs, function(.) {
    int <- try(coda::HPDinterval(.), silent = TRUE)
    if(class(int) == 'try-error') {
      int <- matrix(as.double(NA), nrow = ncol(as.matrix(.)), ncol = 2)
      colnames(int) <- c("lower", "upper")
      int
    } else {
      int
    }
  })
  
  HPD50 <- lapply(mcs, function(.) {
    int <- try(coda::HPDinterval(., 0.5), silent = TRUE)
    if(class(int) == 'try-error') {
      int <- matrix(as.double(NA), nrow = ncol(as.matrix(.)), ncol = 2)
      colnames(int) <- c("lower", "upper")
      int
    } else {
      int
    }
  })
  
  # Mode <- list()
  # for(i in 1:length(mcs)) {
  #   Mode[[i]] <- {
  #     par <- mcs[[i]][which.max(log.ps[[i]][, 1]), ]
  #     setattr(par, 'logpost', max(log.ps[[i]][, 1]))
  #     par
  #   }
  # }


  Mean <- sapply(1:length(mcs), function(i) {
    colMeans(mcs[[i]])
  })
  
  if(as.dt) {
    data.table(
      chain=1:length(mcs), stat=statName, start=start, end=end, thinMCMC=thinMCMC, 
      ESS=ESS, Mean=Mean, HPD=HPD, HPD50=HPD50, mcs=mcs)
  } else {
    l <- list(stat=statName, start=start, end=end, thinMCMC=thinMCMC, 
              ESS=ESS, Mean=Mean, HPD=HPD, HPD50=HPD50, G.R.=gelman, mcs=mcs, log.ps=log.ps)
    
    # without this we are running into trouble:
    names(l) <- c('stat', 'start', 'end', 'thinMCMC',
                  'ESS', Mean=Mean, 'HPD', 'HPD50', 'G.R.', 'mcs', 'log.ps')    
    l
  }
}

#' MCMC-sampling from a posterior distribution of a P(OU)MM model given tree, 
#' values at the tips and a prior distribution
#' 
#' @param z Either a numeric vector containing the phenotypic values at the tips
#'   of tree or a named list containing named elements z - a numeric vector and 
#'   tree - a phylo object.
#' @param tree a phylo object
#' @param distgr character vector, see parameter distgr of the
#'   likPOUMMGivenTreeVTips function
#' @param divideEdgesBy numeric, by which all branch lengths in the tree are to 
#'   be divided prior to likelihood calculation. Only the branch lengths in the 
#'   tree are changed and appropriate rescaling of the OU parameters alpha and 
#'   sigma is not made. For example, if this parameter is set to 100, the 
#'   resulting parameter alpha should be divided by 100 and the parameter sigma 
#'   should be divided by 10.
#' @param parInit a function(chainNumber) returning the starting point of the 
#'   MCMC as a vector.
#' @param parPriorMCMC a function(numeric-vector) returning the log-prior of the 
#'   supplied vector
#' @param parMapping a function(numeric-vector) transforming a sampled vector on
#'   the scale of the parameters alpha, theta, sigma, sigmae
#' @param parScaleMCMC numeric matrix indicating the initial jump-distribution matrix
#' @param nSamplesMCMC integer indicating how many iterations should the mcmc-chain 
#'   contain
#' @param nAdaptMCMC integer indicating how many initial iterations should be used 
#'   for adaptation of the jump-distribution matrix
#' @param thinMCMC integer indicating the thinning interval of the mcmc-chain
#' @param accRateMCMC (MCMC) numeric between 0 and 1 indicating the target 
#'   acceptance rate Passed on to adaptMCMC::MCMC.
#' @param gammaMCMC (MCMC) controls the speed of adaption. Should be between 0.5 and
#'   1. A lower gammaMCMC leads to faster adaption. Passed on to adaptMCMC::MCMC.
#' @param nChainsMCMC integer indicating the number of chains to run. Defaults to 
#'   1.
#' @param samplePriorMCMC logical indicating if only the prior distribution should 
#'   be sampled. This can be useful to compare with mcmc-runs for an overlap 
#'   between prior and posterior distributions.
#' @param ... Additional arguments passed to likPOUMMGivenTreeVTips. If ...
#'   includes debug = TRUE, some debug messages will be written also outside of
#'   the call to likPOUMMGivenTreeVTips.
#' @param zName,treeName characters used when the parameter z is a list; 
#'   indicate the names in the list of the values-vector and the tree. Default: 
#'   'z' and 'tree'.
#'   
#' @details Currently, this function calls the MCMC function from the adaptMCMC 
#'   package.
#'   
#' @return a list of coda objects
#' @importFrom  foreach foreach %dopar% %do%
#'   
#' @aliases mcmcPOUMMGivenPriorTreeVTips mcmc.poumm
mcmcPOUMMGivenPriorTreeVTips <- mcmc.poumm <- function(
  loglik, fitML = NULL,
  parMapping, parInitMCMC, parPriorMCMC, parScaleMCMC, nSamplesMCMC, nAdaptMCMC, thinMCMC, accRateMCMC, 
  gammaMCMC, nChainsMCMC, samplePriorMCMC, ..., 
  verbose = FALSE, parallelMCMC = FALSE, zName = 'z', treeName = 'tree') {
  
  samplePriorMCMC <- c(samplePriorMCMC, rep(FALSE, nChainsMCMC - 1))
  
  if(!is.null(list(...)$debug)) {
    debug <- list(...)$debug
  } else {
    debug <- FALSE
  }
  
  post <- function(par, memoMaxLoglik, chainNo) {
    pr <- parPriorMCMC(par)
    if(debug) {
      cat("par:", toString(round(par, 6)), 'prior:', pr, ";")
    }
    if(is.nan(pr) | is.na(pr)) {
      warning(paste0("NA or NaN prior for ", toString(par)))
      c(-Inf, NA, NA, NA)
    } else if(is.infinite(pr)) {
      c(pr, NA, NA, NA)
    } else {
      if(samplePriorMCMC[chainNo]) {
        ll <- 0
        attr(ll, "g0") <- par[5]
        attr(ll, "g0LogPrior") <- NA
      } else {
        ll <- memoiseMax(
          loglik, par = par, memo = memoMaxLoglik, verbose = verbose)
        
        if(debug) {
          cat("par: ", toString(round(par, 6)), ";", "ll:", ll)
        }
      }
      if(debug) {
        cat("\n")
      }
      c(pr + ll, ll, attr(ll, "g0"), attr(ll, "g0LogPrior"))
    }
  }
  
  doMCMC <- function(i) {
    memoMaxLoglik <- new.env()
    
    if(verbose) {
      cat('Chain No:', i, ', nSamplesMCMC = ', nSamplesMCMC, ', initial state: \n')
    }
    
    init <- parInitMCMC(i, fitML)
    
    if(verbose) {
      print(init)  
    }
    
    ch <- convertToMCMC(
      MCMC(
        post, n = nSamplesMCMC, init = init, scale = parScaleMCMC, adapt = nAdaptMCMC, 
        acc.rate = accRateMCMC, gamma = gammaMCMC, 
        memoMaxLoglik = memoMaxLoglik, chainNo = i), 
      thinMCMC = thinMCMC)
    
    print(paste("Finished chain no:", i))
    
    ch$parScaleMCMC.start <- parScaleMCMC    
    
    ch$parMaxLoglik <- mget("par", envir = memoMaxLoglik, 
                            ifnotfound = list(NULL))$par
    ch$valueMaxLoglik <- mget("val", envir = memoMaxLoglik,
                              ifnotfound = list(-Inf))$val
    
    ch
  }
  
  if(parallelMCMC) {
    chains <- foreach::foreach(
      i = 1:nChainsMCMC, .packages = c('coda', 'adaptMCMC', 'poumm')) %dopar% {
        chain <- try(doMCMC(i), silent = TRUE)
      }
  } else {
    chains <- foreach::foreach(
      i = 1:nChainsMCMC, .packages = c('coda', 'adaptMCMC', 'poumm')) %do% {
        chain <- try(doMCMC(i), silent = TRUE)
      }
  }
  
  for(i in 1:length(chains)) {
    if(class(chains[[i]]) == "try-error") {
      warning(paste0("Error in MCMC chain no ", i, ":", toString(chains[[i]])))
    }
  }
  # maximum log-likelihood from each chain. This is used to correct ML fit in
  # case it got stuck in a local optimum.
  maxLogliksFromChains <- sapply(chains, function(.) .$valueMaxLoglik)
  
  list(chains = chains, 
       post = post, 
       parMaxLoglik = chains[[which.max(maxLogliksFromChains)]]$parMaxLoglik,
       valueMaxLoglik = max(maxLogliksFromChains)
       )
}

#' (Adaptive) Metropolis Sampler
#' 
#' @description Copied and modified from Andreas Scheidegger's adaptMCMC
#'   package. The only difference is that the function p can return a vector
#'   instead of a single numeric. The first element of this vector is the 
#'   log-posterior, while the other elements can be any additional information
#'   such as the log-likelihood, or the root-value for that log-likelihood in 
#'   the case of POUMM.
#' @return Same list as the one returned from adaptMCMC::MCMC but the member 
#' log.p is a matrix instead of a vector.
MCMC <- function (p, n, init, scale = rep(1, length(init)), adapt = !is.null(acc.rate), 
          acc.rate = NULL, gamma = 0.5, list = TRUE, n.start = 0,  ...) 
{
  if (adapt & !is.numeric(acc.rate)) 
    stop("Argument \"acc.rate\" is missing!")
  if (gamma < 0.5 | gamma > 1) 
    stop("Argument \"gamma\" must be between 0.5 and 1!")
  if (is.numeric(adapt)) 
    n.adapt <- adapt
  if (adapt == TRUE) 
    n.adapt <- Inf
  if (adapt == FALSE) 
    n.adapt <- 0
  d <- length(init)
  X <- matrix(NA, ncol = d, nrow = n)
  colnames(X) <- names(init)
  X[1, ] <- init
  
  p.val.init <- p(X[1, ], ...)
  
  p.val <- matrix(NA, nrow = n, ncol=length(p.val.init))
  p.val[1, ] <- p.val.init
  
  if (length(scale) == d) {
    M <- diag(scale)
  }
  else {
    M <- scale
  }
  if (ncol(M) != length(init)) {
    stop(paste("Length or dimension of 'init' and 'scale' do not match!", 
               "init:", toString(init), "scale:", toString(M)))
  }
  if (length(init) == 1) {
    stop("One-dimensional sampling is not possible!")
  }
  S <- t(chol(M))
  cat("  generate", n, "samples \n")
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  update.step <- max(5, floor(n/100))
  k <- 0
  for (i in 2:n) {
    if (i%%update.step == 0) {
      setTxtProgressBar(pb, i)
    }
    U <- rnorm(d)
    X.prop <- c(X[i - 1, ] + S %*% U)
    names(X.prop) <- names(init)
    p.val.prop <- p(X.prop, ...)
    
    alpha <- min(1, exp(p.val.prop[1] - p.val[i - 1, 1]))
    if (!is.finite(alpha)) 
      alpha <- 0
    if (runif(1) < alpha) {
      X[i, ] <- X.prop
      p.val[i, ] <- p.val.prop
      k <- k + 1
    }
    else {
      X[i, ] <- X[i - 1, ]
      p.val[i, ] <- p.val[i - 1, ]
    }
    ii <- i + n.start
    if (ii < n.adapt) {
      adapt.rate <- min(5, d * ii^(-gamma))
      M <- S %*% (diag(d) + adapt.rate * (alpha - acc.rate) * 
                    U %*% t(U)/sum(U^2)) %*% t(S)
      eig <- eigen(M, only.values = TRUE)$values
      tol <- ncol(M) * max(abs(eig)) * .Machine$double.eps
      if (!isSymmetric(M) | is.complex(eig) | !all(Re(eig) > 
                                                   tol)) {
        M <- as.matrix(Matrix::nearPD(M)$mat)
      }
      S <- t(chol(M))
    }
  }
  close(pb)
  acceptance.rate <- round(k/(n - 1), 3)
  if (list) {
    return(list(samples = X, log.p = p.val, cov.jump = M, 
                n.sample = n, acceptance.rate = acceptance.rate, 
                adaption = adapt, sampling.parameters = list(sample.density = p, 
                                                             acc.rate = acc.rate, gamma = gamma)))
  }
  else {
    cat("Acceptance rate:", acceptance.rate, "\n")
    return(X)
  }
}
