#' Summarize the results of a P(OU)MM-fit
#' 
#' @param object a POUMM object returned by POUMM or PMM functions.
#' @param ... Not used, but declared for consistency with the generic method summary.
#' @param startMCMC,endMCMC integers indicating the range of the MCMC chains
#' to be used for the analysis (excluding the initial warm-up phase)
#' @param thinMCMC thinning interval of the MCMC chain to avoid strong 
#' autocorrelation between sampled elements;
#' @param stats a named list of functions of the form function(par) { number },
#' which are called for each sample of each mcmc chain in object. Defaults to 
#' a call of statistics(object) returning a list of statistics functions relevant for 
#' the type of object (either a POUMM or a PMM). See also statistics.
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
                          stats = statistics(object),
                          mode = c('short', 'long', 'expert')) {
  mode <- tolower(mode)
  
  tipTimes <- nodeTimes(object$tree, tipsOnly = TRUE)
  tMax <- max(tipTimes)
  tMean <- mean(tipTimes)
  
  parLower <- matrix(object$spec$parLower, nrow = 1)
  parUpper <- matrix(object$spec$parUpper, nrow = 1)
  parML <- matrix(object$fitML$par, nrow = 1)
  
  anlist <- lapply(1:length(stats), function(i) {
    data.table(stat = names(stats)[i], estML = stats[[i]](parML))
  })
  
  anlist <- c(anlist, list(
    data.table(stat = "logpost", estML = NA),
    data.table(stat = "loglik", estML = object$fitML$value),
    data.table(stat = "g0", estML = attr(object$fitML$value, "g0"))
  ))
  
  an.ML <- rbindlist(anlist)
  an.ML[, N:=object$N]
  setcolorder(an.ML, c('stat', 'N', 'estML'))
  
  if(!is.null(object$fitMCMC)) {
    if(is.na(startMCMC)) {
      startMCMC <- object$spec$nSamplesMCMC / 10
    } 
    if(is.na(endMCMC)) {
      endMCMC <- object$spec$nSamplesMCMC
    } 
    
    anlist <- lapply(1:length(stats), function(i) {
      analyseMCMCs(object$fitMCMC$chains, 
                   stat=stats[[i]], statName=names(stats)[i],
                   start=startMCMC, end=endMCMC, thin=thinMCMC, 
                   as.dt=TRUE)
    })
    
    anlist <- c(anlist, list(
      analyseMCMCs(object$fitMCMC$chains, 
                   stat=NULL, statName='logpost',
                   start=startMCMC, end=endMCMC, thin=thinMCMC, 
                   as.dt=TRUE),
      analyseMCMCs(object$fitMCMC$chains, 
                   stat=NULL, statName='loglik', logprior=object$spec$parPrior,
                   start=startMCMC, end=endMCMC, thin=thinMCMC, 
                   as.dt=TRUE),
      analyseMCMCs(object$fitMCMC$chains, 
                   stat=NULL, statName='g0', logprior=object$spec$parPrior,
                   start=startMCMC, end=endMCMC, thin=thinMCMC, 
                   as.dt=TRUE)
    ))
    
    an.MCMC <- rbindlist(anlist)
    an.MCMC[, samplePriorMCMC:=c(object$spec$samplePriorMCMC, 
                                 rep(FALSE, object$spec$nChainsMCMC - 1))]
    
    if(mode[1]!='expert') {
      an.MCMC <- an.MCMC[, 
                         list(
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
                         by=list(stat, samplePriorMCMC)]
    }
    
    if(mode[1] == 'short') {
      an.MCMC <- an.MCMC[
        samplePriorMCMC == FALSE, 
        list(stat, HPD, ESS, G.R.)]
    } else if(mode[1] == 'long') {
      an.MCMC <- an.MCMC[
        samplePriorMCMC == FALSE, 
        list(stat, Mean, HPD, HPD50, start, end, 
             thin = thinMCMC, ESS, G.R., nChains, mcmc, chain)]
    } else if(mode[1] == 'expert') {
      an.MCMC <- an.MCMC[, list(stat, samplePriorMCMC, 
                                Mean, HPD, HPD50, start, end, thin = thinMCMC,
                                ESS, mcmc = mcs, chain)]
    } else {
      warning(paste('mode should be one of "short", "long" or "expert", but was', mode[1], '.'))
    }
    
  } else {
    an.MCMC <- NULL
  }
  
  if(mode[1]%in%c('short', 'long')) {
    if(!is.null(an.ML) & !is.null(an.MCMC)) {
      res <- merge(an.ML, an.MCMC, by = 'stat', all = TRUE, sort = FALSE)
      res[sapply(HPD, is.null), HPD:=list(list(as.double(c(NA, NA))))]
      if(mode[1] == 'long') 
        res[sapply(HPD50, is.null), HPD50:=list(list(as.double(c(NA, NA))))]
    } else if(!is.null(an.ML)) {
      res <- an.ML
    } else if(!is.null(an.MCMC)) {
      res <- an.MCMC
    }
  } else {
    res <- list(spec = object$spec, startMCMC = startMCMC, endMCMC = endMCMC,
                thinMCMC = thinMCMC,
                ML = an.ML, MCMC = an.MCMC)
  }
  
  class(res) <- c('POUMM.summary', class(res))
  res
}


#' Extract statistics from sampled or inferred parameters of a 
#' POUMM/PMM fit
#' @param object An object of class "POUMM" or "PMM".
#'
#' @details This is a generic method.
#' @export
statistics <- function(object) {
  UseMethod('statistics')
}

#' @describeIn statistics Relevant statistics from the sampled parameters of a
#'   POUMM fit
#' @export
statistics.POUMM <- function(object) {
  listPar <- sapply(1:length(object$spec$parLower), function(i) {
    name <- names(object$spec$parLower)[i]
    stat <- eval(
      parse(text=paste0("list(", name, " = function(par) par[, ", i , "])"))
      )
  })
  listOtherStats <- list(
    H2e = function(par) H2e(z = object$z,
                            sigmae = object$spec$parMapping(par)[, 'sigmae']),
    
    H2tInf = function(par) H2(alpha = object$spec$parMapping(par)[, 'alpha'],
                              sigma = object$spec$parMapping(par)[, 'sigma'],
                              sigmae = object$spec$parMapping(par)[, 'sigmae'],
                              t = Inf),
    
    H2tMax = function(par) H2(alpha = object$spec$parMapping(par)[, 'alpha'],
                              sigma = object$spec$parMapping(par)[, 'sigma'],
                              sigmae = object$spec$parMapping(par)[, 'sigmae'],
                              t = object$tMax),
    
    H2tMean = function(par) H2(alpha = object$spec$parMapping(par)[, 'alpha'],
                               sigma = object$spec$parMapping(par)[, 'sigma'],
                               sigmae = object$spec$parMapping(par)[, 'sigmae'],
                               t = object$tMean),
    
    alpha = function(par) object$spec$parMapping(par)[, 'alpha'],
    theta = function(par) object$spec$parMapping(par)[, 'theta'],
    sigma = function(par) object$spec$parMapping(par)[, 'sigma'],
    sigmae = function(par) object$spec$parMapping(par)[, 'sigmae'],
    sigma2 = function(par) object$spec$parMapping(par)[, 'sigma']^2,
    sigmae2 = function(par) object$spec$parMapping(par)[, 'sigmae']^2,
    
    sigmaG2tMean = function(par) varOU(alpha = object$spec$parMapping(par)[, 'alpha'],
                                       sigma = object$spec$parMapping(par)[, 'sigma'],
                                       t = object$tMean),
    
    sigmaG2tMax = function(par) varOU(alpha = object$spec$parMapping(par)[, 'alpha'],
                                      sigma = object$spec$parMapping(par)[, 'sigma'],
                                      t = object$tMax),
    
    sigmaG2tInf = function(par) varOU(alpha = object$spec$parMapping(par)[, 'alpha'],
                                      sigma = object$spec$parMapping(par)[, 'sigma'],
                                      t = Inf))
  c(listPar, listOtherStats[setdiff(names(listOtherStats), names(listPar))])
}

