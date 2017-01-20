
# Implementation of the POUMM likelihood and heritability estimators

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
#' @param ... additional arguments passed to the dVGivenTreeOU function 
#'   (?dVGivenTreeOU for details).
#' @param spec A named list specifying how the ML and MCMC fit should be done.
#'   See ?specifyPOUMM.
#' @param doMCMC logical: should a MCMC fit be performed. This is TRUE by 
#'   default. Note, however, that, to obtain meaningful estimates MCMC may need 
#'   to run for hours, especially on bigger trees. An MCMC fit provides a sample
#'   from the posterior distribution of the parameters given a prior 
#'   distribution and the data. Unlike the ML-fit, it allows to estimate 
#'   confidence intervals for the estimated parameters.
#' @param parallelMCMC Logical: should the MCMC chains be run in parallel. 
#'
#' @param verbose,debug Logical flags indicating whether to print informative 
#'  and/or debug information on the standard output (both are set to to FALSE by
#'  default).
#'   
#'   
#' @details The PMM function fits the PMM model to the tree and data
#'   by fixing the parameter alpha to 0.
#'   
#' @return For POUMM, an object of class 'POUMM'. For PMM, an object of class 'PMM'.
#'   
#' @importFrom stats var sd rnorm dnorm dexp rexp dunif runif 
NULL

#' @describeIn POUMM_PMM Maximum likelihood and Bayesian fit of the POUMM
#' 
#' @export
POUMM <- function(
  z, tree, zName = 'z', treeName = 'tree', 
  parDigits = 6, usempfr = 0, 
  ..., 
  spec = NULL, doMCMC = TRUE,
  verbose = FALSE, debug=FALSE) {
  
  ###### Verify input data ######
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
  
  # ######## tree branch-length transformation ##########
  # if(divideEdgesBy != 1) {
  #   tree$edge.length <- tree$edge.length / divideEdgesBy
  # }
  # 
  # 
  # ######## filtering of tree-tips ########
  # if(!is.null(removeTips) & length(removeTips) > 0) {
  #   if(!is.character(removeTips)) {
  #     stop("Remove tips should be a character vector denoting tips to be ",
  #          "removed from the tree before the model fit.")
  #   } else {
  #     pos <- match(removeTips, tree$tip.label)
  #     pos <- pos[!is.na(pos)]
  #     if(verbose) {
  #       cat('Removing tips', toString(tree$tip.label[pos]), '.')
  #     }
  #     if(length(pos)>0) {
  #       names(z) <- tree$tip.label
  #       tree <- ape::drop.tip(tree, pos)
  #       z <- z[tree$tip.label]
  #       result$N <- length(tree$tip.label)
  #     }
  #   }
  # }
  
  ######## Caching pruneInfo for faster likelihood calculations
  pruneInfo <- pruneTree(tree)
  
  tTips <- nodeTimes(tree, tipsOnly = TRUE)
  
  ######## Default POUMM spec ###########
  spec <- do.call(specifyPOUMM, 
                  c(list(zMin = min(z), zMean = mean(z), zMax = max(z), 
                         zVar = var(z), tMin = min(tTips), tMean = mean(tTips),
                         tMax = max(tTips)), spec))
  
  result <- list(z = z[1:length(tree$tip.label)], tree = tree, 
                 N = length(tree$tip.label), 
                 tMax = max(nodeTimes(tree, tipsOnly = TRUE)),
                 tMean = mean(nodeTimes(tree, tipsOnly = TRUE)), 
                 spec = spec, 
                 ...)
  
  
  
  # define a loglik function to be called during likelihood maximization as well
  # as MCMC-sampling. 
  # The argument par is a named numeric vector. This vector is mapped to the 
  # POUMM parameters alpha, theta, sigma, sigmae and g0 using the function 
  # spec$parMapping. If g0 = NA the loglik
  # funciton finds the value of g0 that maximizes the likelihood, i.e. 
  # p(z | tree, alpha, theta, sigma, sigmae, g0), or if g0Prior specifiies a 
  # normal distribution with mean g0Prior$mean and variance g0Prior$var, 
  # p(z | tree, alpha, theta, sigma, sigmae, g0) x N(g0 | g0Prior$mean,
  # g0Prior$var). 
  loglik <- function(par, memo = NULL) {
    if(parDigits >= 0) {
      par <- as.vector(round(par, parDigits))
    }
    
    atsseg0 <- spec$parMapping(par)
    
    val <- 
      do.call(likPOUMMGivenTreeVTips, c(list(z, tree), as.list(atsseg0), list(
        g0Prior = spec$g0Prior, log = TRUE, pruneInfo = pruneInfo, 
        usempfr = usempfr, ...)))
    
    if(is.na(val) | is.infinite(val)) {
      val <- -1e100
      attr(val, "g0") <- atsseg0[5]
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
          do.call(likPOUMMGivenTreeVTips, c(list(z, tree), as.list(atsseg0), list(
            g0Prior = spec$g0Prior, log = TRUE, pruneInfo = pruneInfo, 
            usempfr = usempfr + 1, ...)))
        
        if(verbose) {
          cat('. New: ', val, '.\n')
        }  
      }
    }
    
    val
  }
  
  result$loglik <- loglik
  
  if(verbose) {
    print('Performing ML-fit...')
  }
  fitML <- do.call(maxLikPOUMMGivenTreeVTips, 
                   c(list(loglik = loglik, verbose = verbose, debug = debug), 
                     spec))
  
  if(verbose) {
    cat("max loglik from ML: \n")
    print(fitML$value)
    cat("parameters at max loglik from ML: \n")
    print(fitML$par)
  }
  
  result[['fitML']] <- fitML
  
  if(doMCMC) {
    if(verbose) {
      print('Performing MCMC-fit...')
    }
    
    fitMCMC <- do.call(
      mcmcPOUMMGivenPriorTreeVTips, 
      c(list(loglik = loglik, fitML = fitML, verbose = verbose, debug = debug),
        spec))
    
    if(verbose) {
      cat("max loglik from MCMCs: \n")
      print(fitMCMC$valueMaxLoglik)
      cat("parameters at max loglik from MCMCs: \n")
      print(fitMCMC$parMaxLoglik)
    }
    
    # if the max likelihood from the MCMC is bigger than the max likelihood 
    # from the ML-fit, it is likely that the ML fit got stuck in a local 
    # optimum, so we correct it.
    if(fitML$value < fitMCMC$valueMaxLoglik) {
      parInitML <- as.vector(fitMCMC$parMaxLoglik)
      names(parInitML) <- names(spec$parLower)
      
      if(all(c(parInitML >= spec$parLower, parInitML <= spec$parUpper))) {
        if(verbose) {
          cat("The MCMC-fit found a better likelihood than the ML-fit. Performing ML-fit starting from the MCMC optimum.")
        }
        
        spec[["parInitML"]] <- parInitML
        fitML2 <- do.call(maxLikPOUMMGivenTreeVTips, 
                         c(list(loglik = loglik, verbose = verbose), spec))
        
        result[['fitML']] <- fitML2
        
      } else {
        message <- "The MCMC-fit found a better likelihood outside of the search-region of the ML-fit."
        cat(message)
        warning(message)
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
  parDigits = 6, usempfr = 0,
  ...,
  spec = NULL, doMCMC = TRUE, 
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
  
  tTips <- nodeTimes(tree, tipsOnly = TRUE)
  
  zMin <- min(z); zMean <- mean(z); zMax <- max(z); zVar <- var(z); zSD <- sd(z);
  tMin <- min(tTips); tMean <- mean(tTips); tMax <- max(tTips);
  
  spec <- do.call(specifyPMM, 
                  c(list(zMin, zMean, zMax, zVar, zSD, tMin, tMean, tMax), spec))
  
  
  res <- POUMM(
    z = z, tree = tree, zName = zName, treeName = treeName, 
    parDigits = parDigits, usempfr = usempfr, ..., 
    spec = spec, doMCMC = doMCMC, 
    verbose = verbose)
  class(res) <- c('PMM', class(res))
  res
}

#' Extract maximum likelihood and degrees of freedom from a fitted POUMM model
#' @param object An object of class POUMM.
#' @export
logLik.POUMM <- function(object, ...) {
  if("POUMM" %in% class(object)) {
    lik <- object$fitML$value
    attr(lik, "df") <- length(object$spec$parLower)
    attr(lik, "nobs") <- object$N
    class(lik) <- "logLik"
    lik
  } else {
    stop("logLik.POUMM called on non POUMM-object.")
  }
}

#' Extract maximum likelihood fitted parameters (coefficients) from a fitted 
#' POUMM  model.
#' 
#' @param object An object of class POUMM.
#' @param mapped Logical indicating whether the standard POUMM parameters should
#' also be extracted.
#' @details The parameters extracted are the ones used as input to the model's
#' parMapping function. 
#' 
#' @return A named vector with the fitted parameters of the model.
#' @export
coef.POUMM <- function(object, mapped = FALSE) {
  if("POUMM" %in% class(object)) {
    pars <- object$fitML$par
    g0 <- attr(object$fitML$value, "g0")
    if(mapped) {
      pars <- object$spec$parMapping(pars)
    }
    if("g0" %in% names(pars)) {
      pars["g0"] <- g0
    }
      
    pars
  } else {
    stop("coef.POUMM called on non POUMM-object.")
  }
}

#' Extract maximum likelihood expected genotypic values at the tips of a tree, 
#' to which a POUMM model has been previously fitted
#' @param object An object of class POUMM.
#' @param g0 A number specifying the root-genotypic value. It defaults to NA, 
#'   which would cause the maximum likelihood value .
#' @param vCov A logical indicating whether a list with the genotypic values and their variance covariance matrix should be returned or only a vector of the genotypic values (default is FALSE).
#' @return If vCov == TRUE, a list with elements g - the genotypic values and 
#' vCov - the variance-covariance matrix of these values for the specific tree, 
#' observed values z and POUMM ML-fit. If vCov == FALSE, only the vector of genotypic values corresponding to the tip-labels in the tree is returned.
#' @export
fitted.POUMM <- function(object, g0 = coef.POUMM(object, mapped=TRUE)['g0'], vCov=FALSE) {
  if("POUMM" %in% class(object)) {
    if(is.nan(g0)) {
      warning("Genotypic values cannot be inferred for g0=NaN; Read documentaton for parMapping in ?specifyPOUMM and use a parMapping function that sets a finite value or NA for g0.")
    }
    p <- coef(object, mapped = TRUE)
    gList <- gPOUMM(object$z, object$tree, g0, 
                    p["alpha"], p["theta"], p["sigma"], p["sigmae"])
    g = as.vector(gList$mu.g.poumm)
    names(g) <- object$tree$tip.label
    if(vCov) {
      list(g = g, vCov = gList$V.g.poumm)
    } else {
      g
    }
  } else {
    stop("fitted.POUMM called on non POUMM-object.")
  }
}


#' Extract maximum likelihood environmental contributions (residuals) at the tips of a tree, to which a POUMM model has been fitted.
#' @param object An object of class POUMM.
#' @param g0 A number specifying the root-genotypic value. It defaults to NA, 
#'   which would cause the maximum likeliheood value .
#' @param vCov A logical indicating whether a list with the genotypic values and their variance covariance matrix should be returned or only a vector of the genotypic values (default is FALSE).
#' @return The vector of e-values (residuals) corresponding to the tip-labels in the tree.
#' @export
residuals.POUMM <- function(object, g0 = coef.POUMM(object, mapped=TRUE)['g0']) {
  if("POUMM" %in% class(object)) {
    if(is.nan(g0)) {
      warning("Environmental contributions cannot be inferred for g0=NaN; Read documentaton for parMapping in ?specifyPOUMM and use a parMapping function that sets a finite value or NA for g0.")
    }
    e <- object$z - fitted(object, g0)
    names(e) <- object$tree$tip.label
    e
  } else {
    stop("residuals.POUMM called on a non POUMM-object.")
  }
}

#' Plots of a POUMM-fit
#' @param object An object of class POUMM.
#' @param type A character indicating the type of plot(s) to be generated.
#'   Defaults to "MCMC", resulting in a trace and density plot for the selected
#'   statistics (see argument stat).
#' @param doPlot Logical indicating whether a plot should be printed on the 
#'   currently active graphics device or whether to return a list of ggplot 
#'   objects for further processing. Defaults to TRUE.
#' @param interactive Logical indicating whether the user should press a key 
#'   before generating a next plot (when needed to display two or more plots).
#'   Defaults to TRUE. Meaningless if 
#'   doPlot = FALSE.
#' @param stat A character vector with the names of statistics to be plotted.
#'   These should be names from the stats-list (see argument statFunctions).
#'   Defaults to c("alpha", "theta", "sigma", "sigmae", "H2tMean", "H2tInf").
#' @param chain A vector of integers indicating the chains to be plotted. 
#' @param startMCMC,endMCMC,thinMCMC Integers used to extract a sample from the 
#'   MCMC-chain; passed to summary().
#' @param statFunctions Named list of statistics functions; passed to summary().
#' 
#' @return If doPlot==FALSE, a named list containing a member called data of
#'   class data.table and several members of class ggplot.
#' @export
plot.POUMM <- 
  function(object, type=c("MCMC"), 
           doPlot = TRUE, interactive = TRUE,
           stat=c("alpha", "theta", "sigma", "sigmae", "g0", "H2tMean"),
           chain=NULL,
           startMCMC = NA, endMCMC = NA, thinMCMC = 1000, 
           statFunctions = statistics(object)) {
  
  if("POUMM" %in% class(object)) {
    summ <- summary(object, mode = "expert", 
                    startMCMC = startMCMC, endMCMC = endMCMC, 
                    thinMCMC = thinMCMC, stats = statFunctions)
    plot(summ)
  } else {
    stop("plot.POUMM called on a non POUMM-object.")
  }
}