#' Summarize the results of a POUMM-fit
#' 
#' @param object a POUMM object returned by POUMM-function (see ?POUMM).
#' @param ... Not used, but declared for consistency with the generic method summary.
#' @param startMCMC,endMCMC integers indicating the range of the MCMC chains
#' to be used for the analysis (excluding the initial warm-up phase)
#' @param thinMCMC thinning interval of the MCMC chain to avoid strong 
#' autocorrelation between sampled elements;
#' @param stats a named list of functions of the form function(par) { number },
#' which are called for each sample of each mcmc chain in object. Defaults to 
#' a call of statistics(object) returning a list of statistics functions relevant for 
#' the object. See also statistics.
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
  
  # declare global variables to avoid CRAN CHECK NOTES "no visible binding":
  N <- estML <- samplePriorMCMC <- HPD <- HPD50 <- ESS <- HPDUpperFiltered <- 
    HPDLowerFiltered <- value <- HPDUpper <- HPDLower <- it <- Mean <- mcs <- 
    ESS <- nChains <- chain <- G.R. <- stat <- NULL
  
  mode <- tolower(mode)
  
  tipTimes <- nodeTimes(object$pruneInfo$tree, tipsOnly = TRUE)
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
                   stat = stats[[i]], statName = names(stats)[i],
                   start = startMCMC, end = endMCMC, thinMCMC = thinMCMC, 
                   as.dt = TRUE)
    })
    
    anlist <- c(anlist, list(
      analyseMCMCs(object$fitMCMC$chains, 
                   stat=NULL, statName='logpost',
                   start = startMCMC, end=endMCMC, thinMCMC = thinMCMC, 
                   as.dt = TRUE),
      analyseMCMCs(object$fitMCMC$chains, 
                   stat = NULL, statName='loglik', logprior=object$spec$parPrior,
                   start = startMCMC, end = endMCMC, thinMCMC = thinMCMC, 
                   as.dt = TRUE)
    ))
    
    if( !("g0" %in% names(stats)) ) {
      anlist <- c(anlist, list(
        analyseMCMCs(object$fitMCMC$chains, 
                     stat=NULL, statName='g0', logprior=object$spec$parPrior,
                     start=startMCMC, end=endMCMC, thinMCMC=thinMCMC, 
                     as.dt=TRUE))
      )
    }
    
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
             thin = thinMCMC, ESS, G.R., nChains, mcmc)]
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
  
  #r.sq <- r.squared(object)
  #attr(res, 'r.squared') <- r.sq
  #attr(res, 'adj.r.squared') <- adj.r.squared(object, r.sq = r.sq)
  class(res) <- c('summary.POUMM', class(res))
  res
}

#' Plot a summary of a POUMM fit
#' @param x An object of class POUMM.
#' @param type A character indicating the type of plot(s) to be generated.
#'   Defaults to "MCMC", resulting in a trace and density plot for the selected
#'   statistics (see argument stat). Currently, only 'MCMC' type is supported.
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
#' @param doZoomIn (type MCMC only) A logical value indicating whether the 
#'   produced plots should have a limitation on the x-axis according to an 
#'   expression set in zoomInFilter (see below). Default value is FALSE.
#' @param zoomInFilter A character string which evaluates as logical value. If 
#'   doZoomIn is set to TRUE, this filter is applied to each point in each MCMC
#'   chain and the data-point is filtered out if it evaluates to FALSE. This 
#'   allows to zoomIn the x-axis of density plots but should be use with caution,
#'   since filtering out points from the MCMC-sample can affect the kernel densities.
#'   Default value is 
#'    paste0("(stat %in% c('H2e','H2tMean','H2tInf','H2tMax') |",
#'           " (value >= HPDLower & value <= HPDUpper))"). 
#' @param ... Not used; included for compatibility with the generic function plot.
#' 
#' @return The function returns nothing if doPlot=TRUE; Otherwise, it returns a
#'  list with a data.table corresponding to the data passed to ggplot and two
#'  ggplot objects corresponding to an MCMC-trace and a posterior density plot.
#'  
#' @import ggplot2
#' @import methods
#'  
#' @export
plot.summary.POUMM <- function(
  x, type=c("MCMC"), 
  doPlot = TRUE, interactive = TRUE,
  stat=c("alpha", "theta", "sigma", "sigmae", "g0", "H2tMean"),
  chain=NULL,
  doZoomIn = FALSE,
  zoomInFilter = paste0("(stat %in% c('H2e','H2tMean','H2tInf','H2tMax') |",
                        " (value >= HPDLower & value <= HPDUpper))"), ...) {
  
  # declare global variables to avoid CRAN CHECK NOTES "no visible binding":
  N <- estML <- samplePriorMCMC <- HPD <- HPD50 <- ESS <- HPDUpperFiltered <- 
    HPDLowerFiltered <- value <- HPDUpper <- HPDLower <- it <- Mean <- mcs <- 
    ESS <- nChains <- chain <- G.R. <- NULL
  
  if(class(x) == "summary.POUMM" & !is.null(x$MCMC)) {
    .stat <- stat
    .chain <- chain
    
    data <- merge(x$ML, x$MCMC, by = "stat")
    
    data <- data[
      { 
        if(!is.null(.stat)) {stat %in% .stat} else TRUE 
      } & { 
        if(!is.null(.chain)) {chain %in% .chain} else TRUE
      }]
    
    setkey(data, stat)
    
    data <- data[list(.stat)]
    
    data <- data[{ 
      if(!is.null(.stat)) {stat %in% .stat} else TRUE 
    } & {
      if(!is.null(.chain)) {chain %in% .chain} else TRUE 
    }, list(
      N, estML, 
      samplePriorMCMC, 
      HPDLower = sapply(HPD, function(.) .[1]),
      HPDUpper = sapply(HPD, function(.) .[2]),
      HPD50Lower = sapply(HPD50, function(.) .[1]),
      HPD50Upper = sapply(HPD50, function(.) .[2]), 
      ESS,
      value = unlist(mcmc), 
      it = seq(x$startMCMC, by = x$thinMCMC, along.with = mcmc[[1]])), 
    by = list(stat = factor(stat), chain = factor(chain))]
    
    data <- data[!doZoomIn | eval(parse(text=zoomInFilter))]
    
    data[, HPDUpperFiltered:=min(max(value), unique(HPDUpper)), 
         list(stat = factor(stat), chain = factor(chain))]
    
    data[, HPDLowerFiltered:=max(min(value), unique(HPDLower)), 
         list(stat = factor(stat), chain = factor(chain))]
    
    if(type == "MCMC") {
      traceplot <- ggplot(data) + 
        geom_line(aes(x=it, y=value, col = chain)) + 
        facet_wrap(~stat, scales = "free")
      
      densplot <- ggplot(data) + 
        geom_density(aes(x=value, col = chain)) + 
        geom_segment(aes(x=HPDLowerFiltered, xend=HPDUpperFiltered, 
                         y=0, yend=0, col = chain)) +
        geom_point(aes(x=estML, y=0)) +
        facet_wrap(~stat, scales = "free")
      
      if(doPlot) {
        print(traceplot) 
        if(interactive) {
          print("Press Enter to see the next plot")
          scan("", what = "character", nlines = 1)
        }
        print(densplot) 
      } else {
        list(data = data, traceplot = traceplot, densplot = densplot)  
      }
    }
  } else {
    stop("plot.summary.POUMM called on a non summary.POUMM-object or a missing MCMC element. Verify that summary.POUMM has been called with mode = 'expert'")
  }
}

#' Extract statistics from sampled or inferred parameters of a 
#' POUMM fit
#' @param object An object of class "POUMM".
#'
#' @details This is a generic method.
#' @export
statistics <- function(object) {
  UseMethod('statistics')
}

#' @describeIn statistics Relevant statistics from the sampled parameters of a
#'   POUMM fit
#'  
#' @export
statistics.POUMM <- function(object) {
  listPar <- sapply(1:length(object$spec$parLower), function(i) {
    name <- names(object$spec$parLower)[i]
    stat <- eval(
      parse(text=paste0("list(", name, " = function(par) par[, ", i , "])"))
      )
  })
  listOtherStats <- list(
    H2e = function(par) H2e(z = object$pruneInfo$z,
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

