# Implementation of the POUMM likelihood and heritability estimators


#' @rdname PhylogeneticH2
#' @title Phylogenetic Heritability
#'   
#' @description The phylogenetic heritability, \eqn{H^2}, is defined as the 
#'   ratio of the genetic variance over the total phenotypic variance expected 
#'   at a given evolutionary time t (measured from the root of the tree). Thus,
#'   the phylogenetic heritability connects the parameters alpha, sigma and
#'   sigmae of the POUMM model through a set of equations. The functions
#'   described here provide an R-implementation of these equations.
#'   
#' @param t Numeric value denoting evolutionary time (i.e. distance from the 
#'   root of a phylogenetic tree).
#' @param alpha,sigma Numeric values or n-vectors, parameters of the OU process;
#'   alpha and sigma must be non-negative. A zero alpha is interpreted as the 
#'   Brownian motion process in the limit alpha -> 0.
#' @param sigmae Numeric, environmental phenotypic deviation at the tips.
#' @param H2 Phylogenetic heritability at time t.
#' @param z Numerical vector of observed phenotypes.
#' @param tree A phylo object.
#' @param tFrom,tTo Numerical minimal and maximal root-tip distance to limit the
#'   calculation.
#'   
#' @return All functions return numerical values or NA, in case of invalid 
#'   parameters
#' @examples 
#' # At POUMM stationary state (equilibrium, t=Inf)
#' H2 <- H2(alpha = 0.75, sigma = 1, sigmae = 1, t = Inf)     # 0.4
#' alpha <- alpha(H2 = H2, sigma = 1, sigmae = 1, t = Inf)    # 0.75
#' sigma <- sigma(H2 = H2, alpha = 0.75, sigmae = 1, t = Inf) # 1
#' sigmae <- sigmae(H2 = H2, alpha = 0.75, sigma = 1, t = Inf) # 1
#' 
#' # At finite time t = 0.2
#' H2 <- H2(alpha = 0.75, sigma = 1, sigmae = 1, t = 0.2)     # 0.1473309
#' alpha <- alpha(H2 = H2, sigma = 1, sigmae = 1, t = 0.2)    # 0.75
#' sigma <- sigma(H2 = H2, alpha = 0.75, sigmae = 1, t = 0.2) # 1
#' sigmae <- sigmae(H2  =  H2, alpha = 0.75, sigma = 1, t = 0.2) # 1
#' 
#' # Comparing with the empirical H2e from a simulation
#' N <- 20 
#' tree <- TreeSim::sim.bd.taxa(N, 1, lambda = 2, mu = 1, complete =
#' FALSE)[[1]] 
#' tMean <- mean(nodeTimes(tree, tipsOnly = TRUE))
#' z <- rVNodesGivenTreePOUMM(tree, 8, 2, 2, 1)
#'   
#' phytools::plotBranchbyTrait(tree, z, mode = "nodes", show.tip.label =
#'   FALSE, show.node.label = FALSE) 
#' ape::nodelabels(round(z[-(1:N)], 2)) 
#' ape::tiplabels(round(z[(1:N)], 2))
#'    
#' @name PhylogeneticH2
#' @seealso OU
NULL

#' @describeIn PhylogeneticH2 Calculate alpha given time t, H2, sigma and sigmae
#' 
#' @details The function alpha invokes the gsl function lambert_W0.
#' 
#' @importFrom gsl lambert_W0
#'
#' @export
alpha <- alpha.poumm <- function(H2, sigma, sigmae, t = Inf) {
  if(!(H2>=0 & H2 <=1 & sigma >= 0 & sigmae >=0 & t>=0)) {
    warning(paste0("Function poumm::alpha was called on invalid parameters. ", 
                   "Check that H2>=0 & H2 <=1 & sigma >= 0 & sigmae >=0 & t>=0.",
                   "Parameter values were: ", 
                   toString(c(H2=H2, sigma=sigma, sigmae=sigmae, t=t))))
    NA
  }
  if(is.infinite(t)) {
    if(H2 == 1) {
      # Assume correctly defined PMM, i.e. Brownian motion with sigma > 0 and 
      # environmental deviation with sigmae >= 0
      0
    } else {
      # H2 is the phylogenetic heritability at equilibrium
      sigma^2 * (1 - H2) / (2 * H2 * sigmae^2)    
    }
  } else {
    # finite t
    y <- sigma^2 / sigmae^2 * (1 / H2 - 1)
    if(is.na(y) | is.infinite(y)) {
      as.double(NA)
    } else {
      (t * y + lambert_W0((-exp(-t * y)) * t * y)) / (2 * t)
    }
  }
}

#' @describeIn PhylogeneticH2 Calculate sigma given time t, H2 at time t, alpha
#'   and sigmae
#'   
#' @export
sigma <- sigmaOU.poumm <- function(H2, alpha, sigmae, t=Inf) {
  if(H2>1 | H2<0 | alpha < 0 | sigmae < 0 | t < 0) {
    stop("H2(t) should be in [0, 1], alpha, sigmae and t should be non-negative.")
  }
  res <- if(is.infinite(t)) {
    if(alpha == 0) {
      # BM
      sqrt(sigmae^2 * H2 / (t * (1 - H2)))
    } else {
      # alpha>0, OU in equilibrium
      sqrt(2 * alpha * H2 * sigmae^2 / (1 - H2))
    }
  } else {
    # t is finite
    if(alpha == 0) {
      # BM
      sqrt(sigmae^2 * H2 / (t * (1 - H2)))
    } else {
      # alpha>0, OU not in equilibrium
      sqrt(2 * alpha * H2 * sigmae^2 / ((1 - exp(-2 * alpha * t)) * (1 - H2)))
    } 
  }
  names(res) <- NULL
  res
}

#' @describeIn PhylogeneticH2 Calculate sigmae given alpha, sigma, and H2 at
#'   time t
#' @details The function sigmae uses the formula H2 = varOU(t, alpha, sigma) /
#'   (varOU(t, alpha, sigma) + sigmae^2)
#'   
#' @export
sigmae <- sigmae.poumm <- function(H2, alpha, sigma, t = Inf) {
  sigmaG2 <- varOU(t, alpha, sigma)
  sqrt(sigmaG2 / H2 - sigmaG2)
}

#' @describeIn PhylogeneticH2 "Empirical" phylogenetic heritability estimated
#'   from the empirical variance of the observed phenotypes and sigmae
#'   
#' @importFrom stats var
#' @export
H2e <- H2e.poumm <- function(z, sigmae, tree=NULL, tFrom=0, tTo=Inf) {
  if(!is.null(tree)) {
    tipTimes <- nodeTimes(tree)[1:length(tree$tip.label)]
    z <- z[which(tipTimes >= tFrom & tipTimes <= tTo)]
  }
  1 - (sigmae^2 / var(z))
}

#' Phylogenetic heritability estimated at time t
#' @param alpha,sigma numeric, parameters of the OU process acting on the
#'   genetic contributions
#' @param sigmae numeric, environmental standard deviation
#' @param t time from the beginning of the process at which heritability should
#'   be calculated, i.e. epidemiologic time
#' @param tm average time for within host evolution from getting infected until
#'   getting measured or passing on the infection to another host
#' @export
H2 <- H2.poumm <- function(alpha, sigma, sigmae, t = Inf, tm = 0) {
  lenoutput <- max(sapply(list(alpha, sigma, sigmae, t, tm), length))
  if(lenoutput>1) {
    alpha <- rep(alpha, lenoutput / length(alpha))
    sigma <- rep(sigma, lenoutput / length(sigma))
    sigmae <- rep(sigmae, lenoutput / length(sigmae))
    t <- rep(t, lenoutput / length(t))
    tm <- rep(tm, lenoutput / length(tm))
  }
  sigmag2 <- ifelse(alpha > 0, 
                    0.5*(1 - exp(-2 * alpha * t)) / alpha * sigma^2, 
                    t * sigma^2)  
  sigmagm2 <- ifelse(alpha > 0, 
                     0.5*(1 - exp(-2 * alpha * tm)) / alpha * sigma^2, 
                     tm * sigma^2)
  
  values <- ifelse(t > tm, 
                  (sigmag2 - sigmagm2) / (sigmag2 + sigmae^2), 
                  rep(0, length(t)))
  
  values <- ifelse(alpha == 0 & sigma > 0 & is.infinite(t) & !is.infinite(sigmae),
                   1, values)
  values
}

#' @title The Phylogenetic (Ornstein-Uhlenbeck) Mixed Model
#' @name POUMM_PMM
#'   
#' @description This is the high-level entry point to the POUMM method. The 
#'   POUMM function fits the PMM method to data and returns an object of class 
#'   "POUMM", the PMM function is a shortcut to calling POUMM with alpha fixed 
#'   to 0 and returns an object of class "PMM" which is a subclass of "POUMM". 
#'   The functions described here perform an ML-fit and an MCMC-fit to a 
#'   phylogenetic tree with observed values at its tips.
#'   
#' @param z Either a numeric vector containing the phenotypic values at the tips
#'   of tree or a named list containing named elements z - a numeric vector and 
#'   tree - a phylo object (it is possible to specify different element names 
#'   using the arguments zName and treeName).
#' @param tree A phylo object or NULL in case z is a list.
#' @param zName,treeName Character strings used when the parameter z is a list; 
#'   indicate the names in the list of the values-vector and the tree. Default: 
#'   'z' and 'tree'.
#' @param removeTips A character vector denoting names (tip.labels) of (outlier)
#'   tips to be removed from the tree before the model fit.
#' @param divideEdgesBy Numeric, by which all branch lengths in the tree are to 
#'   be divided prior to model fit (default: 1, meaining that the branch lengths
#'   are left unchanged). This can be used to change the time unit on the branch
#'   lengths, i.e. divide by 24 to switch from hours to days. Note that this
#'   parameter only affects the branch lengths in the tree and appropriate
#'   rescaling of the OU parameters alpha and sigma is not made. Thus,
#'   interpreting the fit parameters on the original time-scale would require
#'   appropriate rescaling. For example, if divideEdgesBy is set to 100, after
#'   the model fit, the resulting parameter alpha should be divided by 100 and
#'   the parameter sigma should be divided by 10. Note alose, that in some 
#'   cases, rescaling the branch lengths can improve the precision of floating 
#'   point operations during the model fit.
#' @param ... additional arguments passed to the dVGivenTreeOU function 
#'   (?dVGivenTreeOU for details).
#' @param verbose A logical indicating whether to print informative messages on 
#'   the standard output.
#'   
#' @param distgr a character or a numeric describing the distribution or fixed 
#'   value to be assumed for the genetic contribution at the root of tree. If a 
#'   character the acceptable values are : 'normal'- normal distribution with 
#'   mean mugr and standard deviation sigmagr; 'maxlik' - fixed value that would
#'   maximise the conditional likelihood on g0.
#'   
#' @param doML logical: should a maximum likelihood fit be performed. This is 
#'   TRUE by default. The following arguments until doMCMC specify how the 
#'   ML-fit should be done. Note that the ML-fit provides only a point-estimate.
#'   To obtain confidence (HPD) intervals, the doMCMC argument should be set to 
#'   TRUE.
#' @param parFixed (ML) a named numeric vector, indicating fixed values for some
#'   of the parameters alpha, theta, sigma and sigmae, by default, c().
#' @param parLower,parUpper (ML) two named numeric vectors indicating the 
#'   boundaries of the search region for the ML-fit. Default values are 
#'   parLower=c(alpha=0, theta=0, sigma=0, sigmae=0) and parUpper=c(alpha=100, 
#'   theta=10, sigma=20, sigmae=10).
#' @param optimizeInSteps (ML) logical indicating whether to split the interval 
#'   (parLower['alpha'], parUpper['alpha']) into sub-intervals and call optimize
#'   separately in each of these intervals. Using this may be of help in case of
#'   multiple local optimums, but results in slower search. Default: FALSE. This
#'   parameter does not have any effect if alpha is a fixed parameter (see 
#'   parFixed).
#' @param tol  (ML) numeric accuracy parameter passed to optimize, defaults to 
#'   0.001. (see Details)
#' @param control (ML) list of parameters passed on to optim, default 
#'   list(factr=1e9), see ?optim.

#' @param doMCMC logical: should a MCMC fit be performed. This is TRUE by 
#'   default. Note, however, that, to obtain meaningful estimates MCMC may need 
#'   to run for hours, especially on bigger trees. An MCMC fit provides a sample
#'   from the posterior distribution of the parameters given a prior 
#'   distribution and the data. Unlike the ML-fit, it allows to estimate 
#'   confidence intervals for the estimated parameters.
#' @param parInit (MCMC) a function(chainNumber) returning the starting point of
#'   the MCMC as a vector.
#' @param parPrior (MCMC) a function(numeric-matrix). Each row in the matrix is 
#'   treated as a vector of parameters. The function should return a numeric 
#'   vector with the log-priors for each row matrix.
#' @param parMapping (MCMC) a function(numeric-matrix) Each row in the matrix is
#'   treated as a vector of parameters. The function should transform each such 
#'   row-vector into a row-vector of the parameters alpha, theta, sigma, sigmae.
#'   The function should return the corresponding matrix. 
#'   \preformatted{ 
#'   function(par) { 
#'    cbind(par[, 1:2, drop = FALSE], 
#'             sqrt(par[, 3:4, drop = FALSE]), 
#'             NA) 
#'   } 
#'   alpha = 3 
#'   theta = 4 
#'   }
#' @param parScale (MCMC) numeric matrix indicating the initial jump-distribution 
#'   matrix
#' @param nSamples (MCMC) integer indicating how many iterations should the 
#'   mcmc-chain contain
#' @param nAdapt (MCMC) integer indicating how many initial iterations should 
#'   be used for adaptation of the jump-distribution matrix
#' @param thin (MCMC) integer indicating the thinning interval of the mcmc-chain
#' @param accRate (MCMC) numeric between 0 and 1 indicating the target 
#'   acceptance rate Passed on to adaptMCMC::MCMC.
#' @param gamma (MCMC) controls the speed of adaption. Should be between 0.5 and
#'   1. A lower gamma leads to faster adaption. Passed on to adaptMCMC::MCMC.
#' @param nChains (MCMC) integer indicating the number of chains to run. 
#'   Defaults to 2 chains, from which the first one is a sample from the prior 
#'   distribution (see samplePrior).
#' @param samplePrior (MCMC) logical indicating if sampling from the prior 
#'   should be done for the first chain (see nChains). This is useful to 
#'   compare mcmc's for an overlap between prior and posterior distributions. 
#'   Default is TRUE.
#'   
#' @details The arguments z, tree, zName, treeName, distgr, divideEdgesBy and 
#'   ... are used by both, the ML- and the MCMC-fit; The arguments after doML 
#'   until doMCMC concern the ML-fit only; the arguments after doMCMC concern 
#'   the MCMC-fit only. The PMM function fits the PMM model to the tree and data
#'   by fixing the parameter alpha to 0.
#'   
#' @return An object of class 'POUMM'.
#'   
#' @importFrom stats var sd rnorm dnorm dexp rexp dunif runif
NULL

#' @describeIn POUMM_PMM Maximum likelihood and Bayesian fit of the POUMM
#' 
#' @export
POUMM <- function(
  z, tree, zName = 'z', treeName = 'tree', 
  removeTips = NULL, divideEdgesBy = 1,
  parDigits = 6, usempfr = 0, 
  ..., 
  doML = TRUE, specML = NULL, doMCMC = TRUE, specMCMC = NULL, 
  verbose = FALSE) {
  
  if(is.list(z)) {
    p <- z
    z <- p[[zName]]
    if(is.null(z)) {
      z <- p[['v']]
    }
    tree <- p[[treeName]]
    
    if(is.null(z) | is.null(tree)) {
      stop('If a list is supplied as argument z, this list should contain a ',
           'vector of trait values named "z" or zName and a phylo-object named ',
           '"tree" or treeName')
    }
    if(any(is.na(z)) | any(is.infinite(z))) {
      stop('Check z for infinite or NA values!')
    }
    if(any(tree$edge.length <= 0) | any(is.infinite(tree$edge.length)) | 
       any(is.na(tree$edge.length))) {
      stop('Check the tree for non-finite or non-positive edge-lengths!')
    }
  }
  
  if(!doMCMC & !doML) {
    stop('Meaningless call to POUMM. ',
         'At least one of doMCMC or doML should be TRUE.')
  } 
  
  specMLDefault <- list(
    parFixed = c(), 
    parLower = c(alpha = 0, theta = 0, sigma = 0, sigmae = 0), 
    parUpper = c(alpha = 100, theta = 10, sigma = 20, sigmae = 10),
    g0 = NA, g0Prior = NULL, 
    tol = 0.001, control = list(factr = 1e8)
  )
  
  if(is.list(specML)) {
    for(name in names(specMLDefault)) {
      if(is.null(specML[[name]])) {
        specML[[name]] <- specMLDefault[[name]]
      } 
    } 
  } else {
    specML <- specMLDefault
  }
  
  specMCMCDefault <- list(
    parMapping = function(par) {
      cbind(par[, 1:2, drop = FALSE], 
            sqrt(par[, 3:4, drop = FALSE]),
            NA) 
    },
    parInit = function(chainNo, fitML = NULL) {
      if(!is.null(fitML)) {
        parML <- c(fitML$par[c('alpha', 'theta')], 
                   fitML$par[c('sigma', 'sigmae')]^2)
        names(parML) <- c('alpha', 'theta', 'sigma2', 'sigmae2')
      } else {
        parML <- NULL
      }
      
      init <- rbind(
        c(alpha = 0, theta = 0, sigma2 = 1, sigmae2 = 1),
        parML,
        c(alpha = 0, theta = 0, sigma2 = 1, sigmae2 = 0),
        c(alpha = 10, theta = 0, sigma2 = 1, sigmae2 = 0),
        c(alpha = 20, theta = 0, sigma2 = 1, sigmae2 = 1),
        c(alpha = 100, theta = 0, sigma2 = 1, sigmae2 = 1))
      
      init[(chainNo - 1) %% nrow(init) + 1, ]
    },
    
    parPrior = function(par) {
      if(any(par[, c(1, 3, 4)] < 0)) {
        -Inf
      } else {
        dexp(par[,1], rate = .01, TRUE) + dunif(par[, 2], 0, 100, TRUE)+
          dexp(par[,3],  rate = .0001, TRUE) + dexp(par[, 4], rate = .01, TRUE)  
      }
    },
    parScale = diag(c(100,  10, 100, 100)),
    nSamples = 5e4, nAdapt = 1e4, thin = 100, accRate = 0.01, gamma = 0.5, 
    nChains = 3, samplePrior = TRUE
  )
  
  if(is.list(specMCMC)) {
    for(name in names(specMCMCDefault)) {
      if(is.null(specMCMC[[name]])) {
        specMCMC[[name]] <- specMCMCDefault[[name]]
      } 
    } 
  } else {
    specMCMC <- specMCMCDefault
  }
  
  result <- list(z = z, tree = tree, removeTips = removeTips, 
                 divideEdgesBy = divideEdgesBy, 
                 N = length(tree$tip.label), 
                 tMax = max(nodeTimes(tree, tipsOnly = TRUE)),
                 tMean = mean(nodeTimes(tree, tipsOnly = TRUE)), 
                 specML = specML, specMCMC = specMCMC, 
                 ...)
  
  if(divideEdgesBy != 1) {
    tree$edge.length <- tree$edge.length/divideEdgesBy
  }
  
  if(!is.null(removeTips) & length(removeTips) > 0) {
    if(!is.character(removeTips)) {
      stop("Remove tips should be a character vector denoting tips to be ",
           "removed from the tree before the model fit.")
    } else {
      pos <- match(removeTips, tree$tip.label)
      pos <- pos[!is.na(pos)]
      if(verbose) {
        cat('Removing tips', toString(tree$tip.label[pos]), '.')
      }
      if(length(pos)>0) {
        names(z) <- tree$tip.label
        tree <- ape::drop.tip(tree, pos)
        z <- z[tree$tip.label]
        result$N <- length(tree$tip.label)
      }
    }
  }
  
  pruneInfo <- pruneTree(tree)
  
  # define a loglik function to be called during likelihood maximization as well
  # as MCMC-sampling. 
  # the argument par is a named numeric vector with elements "alpha", "theta",
  # "sigma" and, optionally, "sigmae" and "g0". If not specified, the default
  # values for the latter two are sigmae = 0 and g0 = NA. If g0 = NA the loglik
  # funciton finds the value of g0 that maximizes the likelihood, i.e. 
  # p(z | tree, alpha, theta, sigma, sigmae, g0), or if g0Prior specifiies a 
  # normal distribution with mean g0Prior$mean and variance g0Prior$var, 
  # p(z | tree, alpha, theta, sigma, sigmae, g0) x N(g0 | g0Prior$mean,
  # g0Prior$var). 
  loglik <- function(par, memo = NULL) {
    if(parDigits >= 0) {
      par <- round(par, parDigits)
    }
    
    val <- 
      do.call(likPOUMMGivenTreeVTips, c(list(z, tree), as.list(par), list(
        g0Prior = specML$g0Prior, log = TRUE, pruneInfo = pruneInfo, 
        usempfr = usempfr, ...)))
    
    if(is.na(val) | is.infinite(val)) {
      val <- -1e100
      attr(val, "g0") <- specML$g0
      attr(val, "g0LogPrior") <- NA
    } 
    
    if(!is.null(memo)) {
      valMemo <- mget('val', memo, ifnotfound = list(-Inf))$val  
    } else {
      valMemo <- NA
    }
    
    if(!is.na(valMemo) & valMemo + 1 < val) {
      if(usempfr == 0) {
        if(verbose) {
          print(par)
          cat('Rmpfr check on: ', val, '>', valMemo)
        }
        
        val <- 
          do.call(likPOUMMGivenTreeVTips, c(list(z, tree), as.list(par), list(
            g0Prior = specML$g0Prior, log = TRUE, pruneInfo = pruneInfo, 
            usempfr = usempfr + 1, ...)))
        
        if(verbose) {
          cat('. New: ', val, '.\n')
        }  
      }
    }
    attr(val, 'par') <- par
    val
  }
  
  result$loglik <- loglik
  
  if(doML) {
    if(verbose) {
      print('Performing ML-fit...')
    }
    fitML <- do.call(maxLikPOUMMGivenTreeVTips, 
                     c(list(loglik = loglik, verbose = verbose), specML))
    
    result[['fitML']] <- fitML
  } 
  
  if(doMCMC) {
    if(verbose) {
      print('Performing MCMC-fit...')
    }
    
    fitMCMC <- do.call(
      mcmcPOUMMGivenPriorTreeVTips, 
      c(list(loglik = loglik, fitML = fitML, verbose = verbose), specMCMC))
    
    if(doML) {
      # if the max likelihood from the MCMC is bigger than the max likelihood 
      # from the ML-fit, it is likely that the ML fit got stuck in a local 
      # optimum, so we correct it.
      if(fitML$value < fitMCMC$valueMaxLoglik) {
        parInit <- as.vector(fitMCMC$parMaxLoglik)
        
        if(all(parInit[1:4] >= specML$parLower) & 
           all(parInit[1:4] <= specML$parUpper)) {
          if(verbose) {
            print("The MCMC-fit found a better likelihood than the ML-fit. Performing ML-fit starting from the MCMC optimum.")
          }
          
          # change the search region for alpha since the search algorithm of
          # optimise doesn't allow to specify a fixed initial value for alpha.
          # This will cause optimise to start from parInit[1], i.e. the optimal
          # alpha found by the MCMC fit.
          
          phi <- 0.61803 # golden section ratio
          a <- specML$parLower["alpha"]
          specML$parUpper["alpha"] <- (parInit[1] - a) / (1 - phi) + a
          specML$parInit <- parInit[1:4]
          
          fitML2 <- do.call(maxLikPOUMMGivenTreeVTips, 
                  c(list(loglik = loglik, verbose = verbose), specML))
          
          result[['fitML']] <- fitML2
          
        } else {
          warning("The MCMC-fit found a better likelihood outside of the search-region of the ML-fit.")
        }
      }
    }
    
    result[['fitMCMC']] <- fitMCMC
  }
  class(result) <- c('POUMM', class(result))
  result
}

#' @describeIn POUMM_PMM Maximum likelihood and Bayesian fit of the PMM
#'   
#' @details The PMM function is a shortcut to calling the POUMM function after 
#'   fixing the parameters alpha and theta to 0. Note that, with alpha == 0, the
#'   value of the parameter theta is irrelevant for the POUMM fit.
#' @export
PMM <- function(
  z, tree, zName = 'z', treeName = 'tree', 
  removeTips = NULL,
  distgr = c('maxlik', 'normal'), divideEdgesBy = 1, parDigits = 6, usempfr = 0,
  ...,
  
  doML = TRUE, 
  parFixed = c(), 
  parLower = c(sigma = 0, sigmae = 0), 
  parUpper = c(sigma = 20, sigmae = 10),
  g0 = NA, g0Prior = NULL, 
  tol = 0.001, control = list(factr = 1e8),  
  verbose = FALSE, 
  
  doMCMC = TRUE, 
  parMapping = function(par) {
    cbind(0, 0, sqrt(par[,1:2,drop = FALSE]), NA)
  },
  parInit = function(chainNo) {
    rbind(c(sigma2 = 1, sigmae2 = 1),
          c(sigma2 = 0, sigmae2 = 1),
          c(sigma2 = 1, sigmae2 = 0),
          c(sigma2 = 10, sigmae2 = 1),
          c(sigma2 = 1, sigmae2 = 10))[chainNo%%5+1, ]},
  parPrior = function(par) {
    if(any(par<0)) { 
      -Inf 
    } else {
      dexp(par[,1],  rate = 10^-4, TRUE)+dexp(par[,2], rate = 10^-4, TRUE)
    }
  },
  parScale = rbind(c(0.02,  0.00),
                c(0.00,  0.02)), 
  nSamples = 2.2e6, nAdapt = 2e5, thin = 100, accRate = 0.01, 
  gamma = 0.5, nChains = 2, samplePrior = TRUE) {
  
  parFixed <-  c(alpha = 0, theta = 0, parFixed)
  parLower <- c(alpha = 0, theta = 0, parLower)
  parUppler <- c(alpha = 100, theta = 10, parUpper)
  
  res <- POUMM(
    z = z, tree = tree, zName = zName, treeName = treeName, 
    removeTips = removeTips, divideEdgesBy = divideEdgesBy, 
    parDigits = parDigits, usempfr = usempfr, ..., 
    
    doML = doML, 
    parFixed = parFixed, 
    parLower = parLower, 
    parUpper = parUpper,
    g0 = g0, g0Prior = g0Prior,
    verbose = verbose,
    
    doMCMC = doMCMC,
    parMapping = parMapping,
    parInit = parInit,
    parPrior = parPrior,
    parScale = parScale,
    nSamples = nSamples, nAdapt = nAdapt, thin = thin, accRate = accRate, 
    gamma = gamma, nChains = nChains, samplePrior = samplePrior)
  class(res) <- c('PMM', class(res))
  res
}

#' Summarize the results of a P(OU)MM-fit
#' 
#' @param object a POUMM object returned by POUMM or PMM functions.
#' @param ... Not used, but declared for consistency with the generic method summary.
#' @param startMCMC,endMCMC integers indicating the range of the MCMC chains
#' to be used for the analysis (excluding the initial warm-up phase)
#' @param thinMCMC thinning interval of the MCMC chain to avoid strong 
#' autocorrelation between sampled elements;
#' @param statsMCMC a named list of functions of the form function(par) { number },
#' which are called for each sample of each mcmc chain in object. Defaults to 
#' a call of statisticsMCMC(object) returning a list of statistics functions relevant for 
#' the type of object (either a POUMM or a PMM). See also statisticsMCMC.
#' @param mode a character indicating the desired format of the returned summary 
#' as follows:
#' 'short' - a data.table with the ML and MCMC estimates of heritability, 
#' model parameters, root-value and other statistics. 
#' 'long' - same information as in 'short' but including also the samples, which
#' can be convenient for 
#' 
#' @import data.table
#' @import coda
#' 
#' @export
summary.POUMM <- function(object, ...,
                          startMCMC = NA, endMCMC = NA, thinMCMC = 1000, 
                          statsMCMC = statisticsMCMC(object),
                          mode = c('short', 'long', 'expert')) {
  mode <- tolower(mode)
  
  tipTimes <- nodeTimes(object$tree, tipsOnly = TRUE)
  tMax <- max(tipTimes)
  tMean <- mean(tipTimes)
  
  if(!is.null(object$fitML)) {
    anlist <- list(
      H2e.ML = data.table(stat = 'H2e', lowerML = 0, upperML = 1, 
                          estML = H2e(z = object$z, sigmae = object$fitML$par['sigmae'])),
      
      H2tMean.ML = data.table(stat = 'H2tMean', lowerML = 0, upperML = 1, 
                              estML = H2(alpha = object$fitML$par[1], 
                                          sigma = object$fitML$par[3], 
                                          sigmae = object$fitML$par[4],
                                          t = object$tMean /
                                            object$divideEdgesBy,
                                          tm = 0)),
      
      H2tMax.ML = data.table(stat = 'H2tMax', lowerML = 0, upperML = 1, 
                             estML = H2(alpha = object$fitML$par[1], 
                                         sigma = object$fitML$par[3], 
                                         sigmae = object$fitML$par[4],
                                         t = object$tMax /
                                           object$divideEdgesBy,
                                         tm = 0)),
      
      H2tInf.ML = data.table(stat = 'H2tInf', lowerML = 0, upperML = 1, 
                             estML = H2(alpha = object$fitML$par[1], 
                                         sigma = object$fitML$par[3], 
                                         sigmae = object$fitML$par[4],
                                         t = Inf)),
      
      loglik.ML = data.table(stat = 'loglik', lowerML = -Inf, upperML = Inf, 
                             estML = object$fitML$value),
      alpha.ML = data.table(stat = 'alpha', 
                            lowerML = object$fitML$parLower[1], 
                            upperML = object$fitML$parUpper[1],
                            estML = object$fitML$par[1]),
      theta.ML = data.table(stat = 'theta', 
                            lowerML = object$fitML$parLower[2], 
                            upperML = object$fitML$parUpper[2],
                            estML = object$fitML$par[2]),
      sigma2.ML = data.table(stat = 'sigma2', 
                             lowerML = object$fitML$parLower[3]^2, 
                             upperML = object$fitML$parUpper[3]^2,
                             estML = object$fitML$par[3]^2),
      sigmae2.ML = data.table(stat = 'sigmae2', 
                              lowerML = object$fitML$parLower[4]^2, 
                              upperML = object$fitML$parUpper[4]^2,
                              estML = object$fitML$par[4]^2),
      sigmaG2tMean.ML = data.table(stat  =  'sigmaG2tMean',
                                   lowerML = varOU(alpha = object$fitML$parUpper[1], 
                                                  sigma = object$fitML$parLower[3], 
                                                  t = object$tMean /
                                                    object$divideEdgesBy),
                                   upperML = varOU(alpha = object$fitML$parLower[1], 
                                                  sigma = object$fitML$parUpper[3], 
                                                  t = object$tMean /
                                                    object$divideEdgesBy),
                                   estML = varOU(alpha = object$fitML$par[1], 
                                                  sigma = object$fitML$par[3], 
                                                  t = object$tMean /
                                                    object$divideEdgesBy)),
      sigmaG2tMax.ML=data.table(stat='sigmaG2tMax',
                                lowerML = varOU(alpha=object$fitML$parUpper[1], 
                                               sigma=object$fitML$parLower[3], 
                                               t=object$tMax /
                                                 object$divideEdgesBy),
                                upperML = varOU(alpha = object$fitML$parLower[1], 
                                               sigma = object$fitML$parUpper[3], 
                                               t = object$tMax /
                                                 object$divideEdgesBy),
                                estML = varOU(alpha = object$fitML$par[1], 
                                               sigma = object$fitML$par[3], 
                                               t = object$tMax /
                                                 object$divideEdgesBy)),
      sigmaG2tInf.ML=data.table(stat = 'sigmaG2tInf',
                                lowerML = varOU(alpha = object$fitML$parUpper[1], 
                                               sigma = object$fitML$parLower[3],
                                               t = Inf),
                                upperML = varOU(alpha = object$fitML$parLower[1], 
                                               sigma = object$fitML$parUpper[3],
                                               t = Inf),
                                estML = varOU(alpha = object$fitML$par[1], 
                                               sigma = object$fitML$par[3],
                                               t = Inf)),
      g0.ML=data.table(stat='g0', lowerML=-Inf, upperML=Inf,
                       estML=object$fitML$g0)
    )
  } else {
    anlist <- list(
      H2e.ML=data.table(stat='H2e', lowerML=0, upperML=1, 
                        estML=NA),
      
      H2tMean.ML=data.table(stat='H2tMax', lowerML=0, upperML=1, 
                            estML=NA),
      
      H2tMax.ML=data.table(stat='H2tMax', lowerML=0, upperML=1, 
                           estML=NA),
      
      H2tInf.ML=data.table(stat='H2tInf', lowerML=0, upperML=1, 
                           estML=NA),
      loglik.ML=data.table(stat='loglik', lowerML=NA, upperML=NA, 
                           estML=NA),
      alpha.ML=data.table(stat='alpha', lowerML=NA, upperML=NA,
                          estML=NA),
      theta.ML=data.table(stat='theta', lowerML=NA, upperML=NA,
                          estML=NA),
      sigma2.ML=data.table(stat='sigma2', lowerML=NA, upperML=NA,
                           estML=NA),
      sigmae2.ML=data.table(stat='sigmae2', lowerML=NA, upperML=NA,
                            estML=NA),
      
      sigmaG2tMean.ML=data.table(stat='sigmaG2tMean', lowerML=NA, upperML=NA, 
                                 estML=NA),
      sigmaG2tMax.ML=data.table(stat='sigmaG2tMax', lowerML=NA, upperML=NA, 
                                estML=NA),
      sigmaG2tInf.ML=data.table(stat='sigmaG2tInf', lowerML=NA, upperML=NA, 
                                estML=NA),
      g0.ML=data.table(stat='g0', lowerML=NA, upperML=NA,
                       estML=NA)
    )
  }
  
  an.ML <- rbindlist(anlist)
  an.ML[, N:=object$N]
  setcolorder(an.ML, c('stat', 'N', 'lowerML', 'upperML', 'estML'))
  
  if(!is.null(object$fitMCMC)) {
    if(is.na(startMCMC)) {
      startMCMC <- object$fitMCMC$nSamples / 10
    } 
    if(is.na(endMCMC)) {
      endMCMC <- object$fitMCMC$nSamples
    } 
    
    anlist <- lapply(1:length(statsMCMC), function(i) {
      analyseMCMCs(object$fitMCMC$chains, 
                   stat=statsMCMC[[i]], statName=names(statsMCMC)[i],
                   start=startMCMC, end=endMCMC, thin=thinMCMC, 
                   as.dt=TRUE)
    })
    
    anlist <- c(anlist, list(
      analyseMCMCs(object$fitMCMC$chains, 
                   stat=NULL, statName='logpost',
                   start=startMCMC, end=endMCMC, thin=thinMCMC, 
                   as.dt=TRUE),
      analyseMCMCs(object$fitMCMC$chains, 
                   stat=NULL, statName='loglik', logprior=object$fitMCMC$parPrior,
                   start=startMCMC, end=endMCMC, thin=thinMCMC, 
                   as.dt=TRUE),
      analyseMCMCs(object$fitMCMC$chains, 
                   stat=NULL, statName='g0', logprior=object$fitMCMC$parPrior,
                   start=startMCMC, end=endMCMC, thin=thinMCMC, 
                   as.dt=TRUE)
    ))
    
    an.MCMC <- rbindlist(anlist)
    an.MCMC[, samplePrior:=object$fitMCMC$samplePrior]
    
    if(mode[1]!='expert') {
      an.MCMC <- an.MCMC[, 
                         list(Mode = Mode[[which.max(sapply(Mode, function(m) 
                           attr(m, 'logpost')))]],
                           Mean = mean(unlist(Mean)),
                           HPD = list(colMeans(do.call(rbind, HPD))),
                           HPD50 = list(colMeans(do.call(rbind, HPD50))),
                           start = start(mcs), 
                           end = end(mcs), 
                           thin = thin(mcs),
                           ESS = sum(unlist(ESS)), 
                           G.R. = if(length(mcs)>1) { 
                             gelman.diag(mcs, autoburnin=FALSE)$psrf[1] 
                           } else {
                             as.double(NA)
                           },
                           nChains = length(mcs),
                           mcmc = mcmc.list(mcmc(do.call(rbind, mcs), thin=thin(mcs)))),
                         by=list(stat, samplePrior)]
    }
    
    if(mode[1] == 'short') {
      an.MCMC <- an.MCMC[samplePrior == FALSE, 
                         list(stat, Mode, Mean, HPD, HPD50, ESS, G.R., nChains)]
    } else if(mode[1] == 'long') {
      an.MCMC <- an.MCMC[samplePrior == FALSE, 
                         list(stat, Mode, Mean, HPD, HPD50, start, end, thin, ESS, G.R., nChains, mcmc)]
    } else if(mode[1] == 'expert') {
      an.MCMC <- an.MCMC[, list(stat, samplePrior, Mode, Mean, HPD, HPD50, start, end, thin, ESS, mcmc = mcs)]
    } else {
      warning(paste('mode should be one of "short", "long" or "expert", but was', mode[1], '.'))
    }
    
  } else {
    an.MCMC <- NULL
  }
  
  if(mode[1]%in%c('short', 'long')) {
    if(!is.null(an.ML) & !is.null(an.MCMC)) {
      res <- merge(an.ML, an.MCMC, by='stat', all=TRUE, sort=FALSE)
      res[sapply(HPD, is.null), HPD:=list(list(c(NA,NA)))]
      res[sapply(HPD50, is.null), HPD50:=list(list(c(NA,NA)))]
    } else if(!is.null(an.ML)) {
      res <- an.ML
    } else if(!is.null(an.MCMC)) {
      res <- an.MCMC
    }
  } else {
    res <- list(specML = object$specML, specMCMC = object$specMCMC, 
                ML = an.ML, MCMC = an.MCMC)
  }
  
  class(res) <- c('POUMM.summary', class(res))
  res
}


#' Extract statistics from an MCMC sample of a POUMM/PMM fit
#' @param object An object of class "POUMM" or "PMM".
#'
#' @details This is a generic method.
#' @export
statisticsMCMC <- function(object) {
  UseMethod('statisticsMCMC')
}

#' @describeIn statisticsMCMC Relevant statistics from the sampled parameters of a
#'   POUMM fit
#' @export
statisticsMCMC.POUMM <- function(object) {
  list(
    H2e = function(par) H2e(z = object$z,
                            sigmae = object$fitMCMC$parMapping(par)[,4]),
    H2tInf = function(par) H2(alpha = object$fitMCMC$parMapping(par)[,1],
                              sigma = object$fitMCMC$parMapping(par)[,3],
                              sigmae = object$fitMCMC$parMapping(par)[,4],
                              t = Inf),
    H2tMax = function(par) H2(alpha = object$fitMCMC$parMapping(par)[,1],
                              sigma = object$fitMCMC$parMapping(par)[,3],
                              sigmae = object$fitMCMC$parMapping(par)[,4],
                              t = object$tMax /
                                object$divideEdgesBy),
    H2tMean = function(par) H2(alpha = object$fitMCMC$parMapping(par)[,1],
                               sigma = object$fitMCMC$parMapping(par)[,3],
                               sigmae = object$fitMCMC$parMapping(par)[,4],
                               t = object$tMean /
                                 object$divideEdgesBy),
    alpha = function(par) object$fitMCMC$parMapping(par)[,1],
    theta = function(par) object$fitMCMC$parMapping(par)[,2],
    sigma2 = function(par) object$fitMCMC$parMapping(par)[,3]^2,
    sigmae2 = function(par) object$fitMCMC$parMapping(par)[,4]^2,
    sigmaG2tMean = function(par) varOU(alpha = object$fitMCMC$parMapping(par)[,1],
                                       sigma = object$fitMCMC$parMapping(par)[,3],
                                       t = object$tMean /
                                         object$divideEdgesBy),
    sigmaG2tMax = function(par) varOU(alpha = object$fitMCMC$parMapping(par)[,1],
                                      sigma = object$fitMCMC$parMapping(par)[,3],
                                      t = object$tMax /
                                        object$divideEdgesBy),
    sigmaG2tInf = function(par) varOU(alpha = object$fitMCMC$parMapping(par)[,1],
                                      sigma = object$fitMCMC$parMapping(par)[,3],
                                      t = Inf))
  
}

#' @describeIn statisticsMCMC Relevant statistics from the sampled parameters of a
#'   PMM fit
#' @export
statisticsMCMC.PMM <- function(object) {
  list(
    H2e = function(par) H2e(z = object$z,
                            sigmae = object$fitMCMC$parMapping(par)[,4]),
    H2tInf = function(par) H2(alpha = object$fitMCMC$parMapping(par)[,1],
                              sigma = object$fitMCMC$parMapping(par)[,3],
                              sigmae = object$fitMCMC$parMapping(par)[,4],
                              t = Inf),
    H2tMax = function(par) H2(alpha = object$fitMCMC$parMapping(par)[,1],
                              sigma = object$fitMCMC$parMapping(par)[,3],
                              sigmae = object$fitMCMC$parMapping(par)[,4],
                              t = object$tMax /
                                object$divideEdgesBy),
    H2tMean = function(par) H2(alpha = object$fitMCMC$parMapping(par)[,1],
                               sigma = object$fitMCMC$parMapping(par)[,3],
                               sigmae = object$fitMCMC$parMapping(par)[,4],
                               t = object$tMean /
                                 object$divideEdgesBy),
    alpha = function(par) rep(0,nrow(par)),
    theta = function(par) rep(0,nrow(par)),
    sigma2 = function(par) object$fitMCMC$parMapping(par)[,3]^2,
    sigmae2 = function(par) object$fitMCMC$parMapping(par)[,4]^2,
    sigmaG2tMean = function(par) varOU(alpha = object$fitMCMC$parMapping(par)[,1],
                                       sigma = object$fitMCMC$parMapping(par)[,3],
                                       t = object$tMean /
                                         object$divideEdgesBy),
    sigmaG2tMax = function(par) varOU(alpha = object$fitMCMC$parMapping(par)[,1],
                                      sigma = object$fitMCMC$parMapping(par)[,3],
                                      t = object$tMax /
                                        object$divideEdgesBy))
  
}
