#' @title Specify parameters for fitting a POUMM model. 
#' 
#' @name specifyPOUMM_PMM
#' 
#' @description Specification and validation of POUMM/PMM settings.
#' specifyPOUMM sets default POUMM settings. 
#' Read the "Specification"-section in \code{vignette("poumm-first-steps")} for 
#' guidelines and examples on how to use this function. 
#'
#' @param zMin,zMean,zMax,zVar,zSD summary statistics of the observed values 
#'   used for constructing default parameter values and limits;
#' @param tMin,tMean,tMax summary statistics of the root-tip disntances in the
#'   phylogenetic tree.
#' @param parMapping An R-function that can handle, both, a numeric vector 
#'   or a numeric matrix as argument. This function should transform the input 
#'   vector or each row-vector (if the input is matrix) into a (row-)vector of 
#'   the POUMM parameters alpha, theta, sigma, sigmae, g0. For a vector input 
#'   the function should return a vector with named elements alpha, theta, 
#'   sigma, sigmae, g0. For a matrix input the function should return a matrix 
#'   with the same number of rows and columns alpha, theta, sigma, sigmae, g0. 
#'   Only finite non-negative values are allowed for alpha, sigma, and sigmae. 
#'   Returning Inf, -Inf, NA or NaN for any of these parameters will result in 
#'   an error during likelihood calculation. Only finite numerical values are 
#'   allowed for theta. The parameter
#'   g0 is treated in a special way and can assume either a finite numerical 
#'   value or one of NA or NaN. If g0 = finite value, this value is used together
#'   with the corresponding values of alpha, theta, sigma, and sigmae for 
#'   likelihood calcuation. If g0 = NA (meaing value Not Avaiable), the value of
#'   g0 is calculated analytically during likelihood calculation in order to 
#'   maximise one of the following: \cr
#'   \enumerate{
#'    \item if a normal prior for g0 was specified (see g0Prior), 
#'       \eqn{pdf(z | \alpha, \theta, \sigma, \sigma_e, g0, tree) x prior(g0)}. 
#'    \item otherwise, \eqn{pdf(z | \alpha, \theta, \sigma, \sigma_e, g0, tree)}. 
#'   }
#'   If g0 = NaN (meaning Not a Number), then the likelihood is marginalized w.r.t. 
#'    the g0's prior distribution (see g0Prior), i.e. the likelihood returned is:
#'    \eqn{pdf(z | \alpha, \theta, \sigma, \sigma_e, tree) = Integral(pdf(z|\alpha,\theta,\sigma,\sigma_e,g0) x pdf(g0) d g0; g0 from -\infty to +\infty) }
#'    In this case (g0=NaN), if g0Prior is not specified, 
#'    it is assumed that g0Prior is the stationary OU normal distribution with 
#'    mean, theta, and variance, varOU(Inf, alpha, sigma). \cr
#'  Examples: \cr
#'   
#'  \preformatted{ 
#'  # Default for POUMM: identity for alpha, theta, sigma, sigmae, NA for g0.
#'  parMapping = function(par) {
#'    if(is.matrix(par)) {
#'      atsseg0 <- cbind(par[, 1:4, drop = FALSE], NA) 
#'      colnames(atsseg0) <- c("alpha", "theta", "sigma", "sigmae", "g0")
#'    } else {
#'      atsseg0 <- c(par[1:4], NA) 
#'      names(atsseg0) <- c("alpha", "theta", "sigma", "sigmae", "g0")
#'    }
#'    atsseg0
#'  }
#'
#'  # Default for PMM: fixed values for alpha and theta; identity for 
#'  # sigma and sigmae, NA for g0.
#'  parMapping = function(par) {
#'    if(is.matrix(par)) {
#'      atsseg0 <- cbind(0, 0, par[, 1:2, drop = FALSE], NA) 
#'      colnames(atsseg0) <- c("alpha", "theta", "sigma", "sigmae", "g0")
#'    } else {
#'      atsseg0 <- c(0, 0, par[1:2], NA) 
#'      names(atsseg0) <- c("alpha", "theta", "sigma", "sigmae", "g0")
#'    }
#'    atsseg0
#'  }
#' }
#' @param parLower,parUpper two named numeric vectors of the same length
#'   indicating the boundaries of the search region for the ML-fit. Calling
#'   parMapping on parLower and parUpper should result in appropriate values of
#'   the POUMM parameters alpha, theta, sigma sigmae and g0. Examples: \cr
#' \preformatted{
#' # Default for POUMM:
#' parLower = c(alpha = 0, theta = zMin - 2 * (zMax - zMin), sigma = 0, sigmae = 0)
#' parUpper = c(alpha = 50, theta = zMax + 2 * (zMax - zMin), 
#'              sigma = poumm::sigma(H2 = .99, alpha = 50, sigmae = zSD, t = tMean), 
#'              sigmae = zSD)
#'
#' # Default for PMM:
#' parLower = c(sigma = 0, sigmae = 0)
#' parUpper = c(sigma = poumm::sigma(H2 = .99, alpha = 0, sigmae=zSD, t = tMean), 
#'              sigmae = sd(z))
#' }
#' 
#' @param g0Prior 
#' @param control List of parameters passed on to optim in the ML-fit, default 
#'   list(factr=1e9), see ?optim.
#'
#' @param parInitMCMC a function(chainNo, fitML) returning an initial state of
#'   an MCMC as a vector. The argument fitML can be used to specify an initial 
#'   state, close to a previously found likelihood optimum. Example: \cr
#' \preformatted{ 
#'  # Default for POUMM:
#'  parInitMCMC = function(chainNo, fitML) {
#'    if(!is.null(fitML)) {
#'      parML <- c(fitML$par[c('alpha', 'theta', 'sigma', 'sigmae')])
#'    } else {
#'      parML <- NULL
#'    }
#'    
#'    init <- rbind(
#'      c(alpha = 0, theta = 0, sigma = 1, sigmae = 0),
#'      parML,
#'      c(alpha = 0, theta = 0, sigma = 1, sigmae = 1)
#'    )
#'    
#'    init[(chainNo - 1) \%\% nrow(init) + 1, ]
#'  }
#'
#'  # Default for PMM:
#'  parInitMCMC = function(chainNo, fitML = NULL) {
#'    if(!is.null(fitML)) {
#'      parML <- c(fitML$par[c('sigma', 'sigmae')])
#'    } else {
#'      parML <- NULL
#'    }
#'    
#'    init <- rbind(
#'      c(sigma = 1, sigmae = 0),
#'      parML,
#'      c(sigma = 1, sigmae = 1)
#'    )
#'    
#'    init[(chainNo - 1) \%\% nrow(init) + 1, ]
#'  }
#' }
#' 
#' @param parPriorMCMC A function of a numeric parameter-vector returning the 
#' log-prior for this parameter vector. Example: \cr
#' 
#' \preformatted{
#' # Default for POUMM:
#'  parPriorMCMC = function(par) {
#'    dexp(par[1], rate = .01, TRUE) + 
#'      dnorm(par[2], zMean, 10 * zSD, TRUE) +
#'      dexp(par[3],  rate = .0001, TRUE) + 
#'      dexp(par[4], rate = .01, TRUE)
#'  }
#'  
#'  # Default for PMM:
#'  parPriorMCMC = function(par) {
#'    dexp(par[1],  rate = .0001, TRUE) + 
#'    dexp(par[2], rate = .01, TRUE)
#'  }
#' }
#' 
#' @param parScaleMCMC Numeric matrix indicating the initial jump-distribution 
#'   matrix for the MCMC fit. Default for POUMM is diag(4); for PMM - diag(2)
#' @param nSamplesMCMC Integer indicating the length of each MCMC chain. Defaults to 1e5.
#' @param nAdaptMCMC Integer indicating how many initial MCMC iterations should 
#'   be used for adaptation of the jump-distribution matrix (see details in 
#'   ?POUMM). Defaults to nSamplesMCMC meaning continuous adaptation throughout
#'   the MCMC. 
#' @param thinMCMC Integer indicating the thinning interval of the mcmc-chains. 
#'   Defaults to 100.
#' @param accRateMCMC numeric between 0 and 1 indicating the target 
#'   acceptance rate of the Metropolis sampling (see details in ?POUMM). Default 0.01.
#' @param gammaMCMC controls the speed of adaption. Should be between 0.5 and
#'   1. A lower gamma leads to faster adaption. Default value is 0.5.
#' @param nChainsMCMC integer indicating the number of chains to run. 
#'   Defaults to 3 chains, from which the first one is a sample from the prior 
#'   distribution (see samplePriorMCMC).
#' @param samplePriorMCMC Logical indicating if sampling from the prior 
#'   should be done for the first chain (see nChainsMCMC). This is useful to 
#'   compare mcmc's for an overlap between prior and posterior distributions. 
#'   Default is TRUE.
#' @param parallelMCMC Logical indicating whether the MCMC chains should be run 
#'   in parallel. Setting this option to TRUE results in using 
#'   \code{foreach::foreach() \%dopar\% { }} construct for the MCMC fit. In order for
#'   parallel execution to be done, you should create a computing cluster and
#'   register it as parallel back-end (see example in package vignette and the 
#'   web-page https://github.com/tobigithub/R-parallel/wiki/R-parallel-Setups).
#'   
#' @return A named list to be passed as a spec argument to POUMM or PMM. 
#' 
#' @export
specifyPOUMM <- function(zMin, zMean, zMax, zVar, zSD = sqrt(zVar), 
                         tMin, tMean, tMax, 
                         parMapping = NULL, 
                         parLower = NULL, parUpper = NULL, 
                         g0Prior = NULL,
                         parInitML = NULL,
                         control = NULL,
                         
                         parPriorMCMC = NULL, 
                         parInitMCMC = NULL, 
                         parScaleMCMC = NULL,
                         nSamplesMCMC = 1e5, nAdaptMCMC = nSamplesMCMC, 
                         thinMCMC = 100, 
                         accRateMCMC = .01, gammaMCMC = 0.5, nChainsMCMC = 3, 
                         samplePriorMCMC = TRUE,
                         parallelMCMC = FALSE
                         ) {
  spec <- list(parMapping = parMapping, 
               parLower = parLower, parUpper = parUpper, 
               g0Prior = g0Prior,
               parInitML = parInitML,
               control = control,
               
               parPriorMCMC = parPriorMCMC, 
               parInitMCMC = parInitMCMC, 
               parScaleMCMC = parScaleMCMC,
               nSamplesMCMC = nSamplesMCMC, nAdaptMCMC = nAdaptMCMC, 
               thinMCMC = thinMCMC, accRateMCMC = accRateMCMC, 
               gammaMCMC = gammaMCMC, nChainsMCMC = nChainsMCMC, 
               samplePriorMCMC = samplePriorMCMC,
               parallelMCMC = parallelMCMC)
  
  specDefault <- list(
    parMapping = function(par) {
      if(is.matrix(par)) {
        atsseg0 <- cbind(par[, 1:4, drop = FALSE], NA) 
        colnames(atsseg0) <- c("alpha", "theta", "sigma", "sigmae", "g0")
      } else {
        atsseg0 <- c(par[1:4], NA) 
        names(atsseg0) <- c("alpha", "theta", "sigma", "sigmae", "g0")
      }
      atsseg0
    },
    
    parLower = c(alpha = 0, theta = zMin - 2 * (zMax - zMin), 
                 sigma = 0, sigmae = 0), 
    
    parUpper = c(alpha = 50, theta = zMax + 2 * (zMax - zMin), 
                 sigma = poumm::sigma(H2 = .99, alpha = 50, sigmae = zSD, t = tMean), 
                 sigmae = zSD),
    
    g0Prior = NULL, 
    
    parInitML = NULL,
    
    control = list(factr = 1e8),
    
    parPriorMCMC = function(par) {
      dexp(par[1], rate = .1, TRUE) + 
        dnorm(par[2], zMean, 5 * zSD, TRUE) +
        dexp(par[3],  rate = .1, TRUE) + 
        dexp(par[4], rate = .1, TRUE)  
      
    },
    
    parInitMCMC = function(chainNo, fitML = NULL) {
      if(!is.null(fitML)) {
        parML <- c(fitML$par[c('alpha', 'theta', 'sigma', 'sigmae')])
        names(parML) <- c('alpha', 'theta', 'sigma', 'sigmae')
      } else {
        parML <- NULL
      }
      
      init <- rbind(
        c(alpha = 0, theta = 0, sigma = 1, sigmae = 0),
        parML,
        c(alpha = 0, theta = 0, sigma = 1, sigmae = 1)
      )
      
      init[(chainNo - 1) %% nrow(init) + 1, ]
    },
    
    parScaleMCMC = diag(4),
    
    nSamplesMCMC = nSamplesMCMC, nAdaptMCMC = nAdaptMCMC, 
    thinMCMC = thinMCMC, accRateMCMC = accRateMCMC, gammaMCMC = gammaMCMC, 
    
    nChainsMCMC = nChainsMCMC, samplePriorMCMC = samplePriorMCMC
  )
  
  if(is.list(spec)) {
    for(name in names(specDefault)) {
      if(is.null(spec[[name]]) & !is.null(specDefault[[name]])) {
        spec[[name]] <- specDefault[[name]]
      } 
    } 
  } else {
    spec <- specDefault
  }
  
  validateSpecPOUMM(spec)
  spec
}

#' @describeIn specifyPOUMM_PMM Specify parameter for fitting a PMM model
#' 
specifyPMM <- function(zMin, zMean, zMax, zVar, zSD = sqrt(zVar), 
                       tMin, tMean, tMax, 
                       parMapping = NULL, 
                       parLower = NULL, parUpper = NULL, 
                       g0Prior = NULL,
                       parInitML = NULL,
                       control = NULL,
                       parPriorMCMC = NULL, 
                       parInitMCMC = NULL, 
                       parScaleMCMC = NULL,
                       nSamplesMCMC = 1e5, nAdaptMCMC = nSamplesMCMC, 
                       thinMCMC = 100, 
                       accRateMCMC = .01, gammaMCMC = 0.5, nChainsMCMC = 3, 
                       samplePriorMCMC = TRUE,
                       parallelMCMC = FALSE) {

  spec <- list(parMapping = parMapping, 
               parLower = parLower, parUpper = parUpper, 
               g0Prior = g0Prior,
               parInitML = parInitML,
               control = control,
               parPriorMCMC = parPriorMCMC, 
               parInitMCMC = parInitMCMC, 
               parScaleMCMC = parScaleMCMC,
               nSamplesMCMC = nSamplesMCMC, nAdaptMCMC = nAdaptMCMC, 
               thinMCMC = thinMCMC, accRateMCMC = accRateMCMC, 
               gammaMCMC = gammaMCMC, nChainsMCMC = nChainsMCMC, 
               samplePriorMCMC = samplePriorMCMC,
               parallelMCMC = parallelMCMC)
  
  specDefault <- specifyPOUMM(
    zMin, zMean, zMax, zVar, zSD, tMin, tMean, tMax, 
    parMapping = function(par) {
      if(is.matrix(par)) {
        atsseg0 <- cbind(0, 0, par[, 1:2, drop = FALSE], NA) 
        colnames(atsseg0) <- c("alpha", "theta", "sigma", "sigmae", "g0")
      } else {
        atsseg0 <- c(0, 0, par[1:2], NA) 
        names(atsseg0) <- c("alpha", "theta", "sigma", "sigmae", "g0")
      }
      atsseg0
    },
    
    parLower = c(sigma = 0, sigmae = 0),
    parUpper = c(
      sigma = poumm::sigma(H2 = .99, alpha = 0, sigmae = zSD, t = tMean), 
      sigmae = zSD),
    
    parPriorMCMC = function(par) {
      dexp(par[1],  rate = .1, TRUE) +
        dexp(par[2], rate = .1, TRUE)
    },

    parInitMCMC = function(chainNo, fitML = NULL) {
      if(!is.null(fitML)) {
        parML <- c(fitML$par[c('sigma', 'sigmae')])
      } else {
        parML <- NULL
      }
      
      init <- rbind(
        c(sigma = 1, sigmae = 0),
        parML,
        c(sigma = 1, sigmae = 1)
      )
      
      init[(chainNo - 1) %% nrow(init) + 1, ]
    },
    parScaleMCMC = diag(2)
  )
  
  if(is.list(spec)) {
    for(name in names(specDefault)) {
      if(is.null(spec[[name]]) & !is.null(specDefault[[name]])) {
        spec[[name]] <- specDefault[[name]]
      } 
    } 
  } else {
    spec <- specDefault
  }
  
  validateSpecPOUMM(spec)
  spec
}

######### Validate specification ###########
validateSpecPOUMM <- function(spec) {
  with(spec, {
    if(nChainsMCMC <= 0) {
      stop(paste0("nChainsMCMC should be greater or equal to 1 but was ", 
                  nChainsMCMC, 
                  ". To disable MCMC-fit, set doMCMC = FALSE in the POUMM call."))
    }
    if(any(c(nSamplesMCMC, nAdaptMCMC, thinMCMC) <= 0)) {
      stop("nSamplesMCMC, nAdaptMCMC and thinMCMC should be positive integers.")
    }
    if(accRateMCMC <= 0 | accRateMCMC > 1) {
      stop("accRateMCMC should be a positive number between 0 and 1.")
    } 
    if(gammaMCMC < 0.5 | gammaMCMC > 1) {
      stop("gammaMCMC should be between 0.5 and 1.")
    }
    if(nChainsMCMC == 1 & samplePriorMCMC) {
      warning("Only one MCMC to be run which samples from the prior.")
    }
    
    
    listParInitMCMC <- lapply(1:spec$nChainsMCMC, function(i) parInitMCMC(i))
    parNames <- names(listParInitMCMC[[1]])
    parLen <- length(listParInitMCMC[[1]])
    
    if(!is.null(parInitML)) {
      if(!is.vector(parInitML) | length(parInitML) != parLen | 
         !identical(names(parInitML), parNames)) {
        stop("parInitML should be a vector of the same length and with the same names as each vector returned by parInitMCMC.")
      }
    }
    
    if(length(parLower) != parLen | length(parUpper) != parLen | 
       !identical(parNames, names(parLower)) | !identical(parNames, names(parUpper))) {
      warning("parInitMCMC returns vectors of different length or with different names from parLower and parUpper.")
    } 
    
    for(i in 1:spec$nChainsMCMC) {
      parInitMCMC <- listParInitMCMC[[i]]
      if(is.null(parInitMCMC) | !is.vector(parInitMCMC) | length(parInitMCMC) != parLen |
         is.null(names(parInitMCMC)) | !identical(names(parInitMCMC), parNames)) {
        print(parInitMCMC)
        stop(paste0("parInitMCMC should return a named vector with the same length and names for each chain index.", " Check parInitMCMC(",i,")."))
      }
      
      atsseg0 <- parMapping(parInitMCMC)
      if(!is.vector(atsseg0) | is.null(names(atsseg0)) | length(atsseg0) != 5 |
         !identical(names(atsseg0), c("alpha", "theta", "sigma", "sigmae", "g0")) |
         atsseg0[1] < 0 | atsseg0[3] < 0 | atsseg0[4] < 0 | any(is.na(atsseg0[1:4]))) {
        print("parInitMCMC : ")
        print(parInitMCMC)
        print(" mapped to : ")
        print(atsseg0)
        stop(paste0("When given a vector, parMapping should return a vector ",
                    "with named elements alpha>=0, theta, sigma>=0, sigmae>=0 and g0.", 
                    "Only g0 can be NA or NaN which has a specical meaning (described in doc for parMapping). Check parInitMCMC(",i,")."))
      }
      
      prior <- try(parPriorMCMC(parInitMCMC), silent = TRUE)
      if(class(prior) != "numeric" | is.na(prior)) {
        print("parInitMCMC : ")
        print(parInitMCMC)
        stop(paste0("parPriorMCMC should return a finite number but returned:", 
                    toString(prior)))
      }
      
      parInitMat <- matrix(parInitMCMC, nrow = 1)
      colnames(parInitMat) <- names(parInitMCMC)
      atsseg0Mat <- parMapping(parInitMat)
      
      
      if(!is.matrix(atsseg0Mat) | is.null(colnames(atsseg0Mat)) | 
         length(atsseg0Mat) != 5 |
         !identical(colnames(atsseg0Mat), c("alpha", "theta", "sigma", "sigmae", "g0")) |
         any(atsseg0Mat[, 1] < 0) | any(atsseg0Mat[, 3] < 0) | 
         any(atsseg0Mat[, 4] < 0) | any(is.na(atsseg0Mat[, 1:4]))) {
        
        cat("parInitMCMC : \n")
        print(parInitMat)
        cat(" mapped to : \n")
        print(atsseg0Mat)
        stop(paste0("When given a matrix, parMapping should return a matrix ",
                    "with named columns alpha>=0, theta, sigma>=0, sigmae>=0 and g0.", 
                    "Only g0 can be NA or NaN which has a specical meaning ",
                    "(described in doc for parMapping)."))
      }
    }
  })
}
