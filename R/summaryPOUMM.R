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
#'   currently active graphics device or whether only to return a list of plot- 
#'   objects for further processing. Defaults to TRUE.
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
#' @param palette A vector of colors (can be character strings) corresponding to the 
#' different chains (in their order 1 (prior), 2, 3). Defaults to c("#999999", 
#' "#0072B2", "#CC79A7", "#E69F00", "#D55E00", "#56B4E9", "#009E73", "#F0E442"),
#' which is a color-blind friendly.
#' @param ... Not used; included for compatibility with the generic function plot.
#' 
#' @return If doPlot==TRUE, the function returns nothing and produces output on 
#' the current graphics device as a side-effect. Otherwise, the function returns
#' a list of plot-objects: traceplot and densplot.
#'  
#' @import ggplot2
#' @importFrom GGally ggpairs print_if_interactive wrap ggally_text
#' @importFrom stats cor complete.cases
#' @import methods
#'  
#' @examples 
#' \dontrun{
#' library(POUMM)
#' 
#' set.seed(1)
#' 
#' N <- 1000
#' 
#' # create a random non-ultrametric tree of N tips
#' tree <- ape::rtree(N)  
#' 
#' # Simulate the evolution of a trait along the tree
#' z <- rVNodesGivenTreePOUMM(
#'   tree, g0 = 8, alpha = 1, theta = 4, sigma = 1.2, sigmae = .8)
#' 
#' fit <- POUMM(z[1:N], tree, spec = list(nSamplesMCMC = 4e5))
#' 
#' # Summarize the results from the fit in a table:
#' summary(fit)
#' 
#' # Create plots for some of the inferred parameters/statistics:
#' pl <- plot(fit, stat = c("alpha", "theta", "sigma", "sigmae", "H2tMean"), 
#'            doZoomIn = TRUE, 
#'            zoomInFilter = paste("!(stat %in% c('alpha', 'sigma', 'sigmae')) |",
#'                                 "(value >= 0 & value <= 8)"),
#'            doPlot = FALSE)
#' 
#' pl$traceplot
#' pl$densplot
#' }
#' 
#' @export
plot.summary.POUMM <- function(
  x, type=c("MCMC"), 
  doPlot = TRUE, 
  stat=c("alpha", "theta", "sigma", "sigmae", "g0", "H2tMean"),
  chain=NULL,
  doZoomIn = FALSE,
  zoomInFilter = paste0("(stat %in% c('H2e','H2tMean','H2tInf','H2tMax') |",
                        " (value >= HPDLower & value <= HPDUpper))"), 
  palette = c("#999999", "#0072B2", "#CC79A7", "#E69F00", "#D55E00", "#56B4E9", "#009E73", "#F0E442"), 
  ...) {
  
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
    
    .availStats <- data[, as.character(unique(stat))]
    .stat <- .availStats[1]
    dtm <- data[stat == .stat, list(stat, it, chain, value)]
    dtm[, (.stat) := value]
    dtm[, c("stat", "value") := NULL]
    
    for(.stat in .availStats[-1]) {
      dtm2 <- 
        data[stat == .stat, 
             eval(parse(text=paste0("list(it, chain, ", .stat, ' = value)')))]
      dtm <- merge(dtm, dtm2, by=c("it", "chain"), all=TRUE)
    }
    
    names(palette) <- as.character(1:length(palette))
    
    my_ggplot <- function(...) ggplot(...) + 
      scale_color_manual(values = palette) +
      scale_fill_manual(values = palette)
    
    if(type == "MCMC") {
      traceplot <- my_ggplot(data) + 
        geom_line(aes(x=it, y=value, col = chain)) + 
        facet_wrap(~stat, scales = "free")
      
      my_dens <- function(data, mapping, ...) {
        my_ggplot(data = data, mapping=mapping) +
          geom_density(..., alpha = 0.5) 
      }
      
      my_points <- function(data, mapping, ...) {
        my_ggplot(data = data, mapping = mapping) + 
          geom_point(..., alpha = 0.5)
      }
      
      # copied and modified from GGally
      my_cor <- function(data, mapping, alignPercent = 0, method = "pearson", 
                use = "complete.obs", corAlignPercent = NULL, corMethod = NULL, 
                corUse = NULL, ...) {
        
        # global variable name to avoid check note:
        labelp <- NULL
        
        # copied from GGally
        is_date <- function (x) {
          inherits(x, c("POSIXt", "POSIXct", "POSIXlt", "Date"))
        }
        
        # copied from GGally
        eval_data_col <- function (data, aes_col) {
          eval(aes_col, data)
        }
        
        # copied from GGally
        str_c <- function (..., sep = "", collapse = NULL) {
          paste(..., sep = sep, collapse = collapse)
        }
        if (!is.null(corAlignPercent)) {
          stop("'corAlignPercent' is deprecated.  Please use argument 'alignPercent'")
        }
        if (!is.null(corMethod)) {
          stop("'corMethod' is deprecated.  Please use argument 'method'")
        }
        if (!is.null(corUse)) {
          stop("'corUse' is deprecated.  Please use argument 'use'")
        }
        useOptions <- c("all.obs", "complete.obs", "pairwise.complete.obs", 
                        "everything", "na.or.complete")
        use <- pmatch(use, useOptions)
        if (is.na(use)) {
          warning("correlation 'use' not found.  Using default value of 'all.obs'")
          use <- useOptions[1]
        }
        else {
          use <- useOptions[use]
        }
        cor_fn <- function(x, y) {
          cor(x, y, method = method, use = use)
        }
        xCol <- deparse(mapping$x)
        yCol <- deparse(mapping$y)
        
        if (is.numeric(eval_data_col(data, mapping$colour))) {
          stop("ggally_cor: mapping color column must be categorical, not numeric")
        }
        colorCol <- deparse(mapping$colour)
        singleColorCol <- ifelse(is.null(colorCol), NULL, paste(colorCol, 
                                                                collapse = ""))
        if (use %in% c("complete.obs", "pairwise.complete.obs", "na.or.complete")) {
          
          if (length(colorCol) > 0) {
            if (singleColorCol %in% colnames(data)) {
              rows <- complete.cases(data[c(xCol, yCol, colorCol)])
            }
            else {
              rows <- complete.cases(data[c(xCol, yCol)])
            }
          } else {
            rows <- complete.cases(data[c(xCol, yCol)])
          }
          if (any(!rows)) {
            total <- sum(!rows)
            if (total > 1) {
              warning("Removed ", total, " rows containing missing values")
            }
            else if (total == 1) {
              warning("Removing 1 row that contained a missing value")
            }
          }
          data <- data[rows, ]
        }
        xVal <- data[[xCol]]
        yVal <- data[[yCol]]
        if (length(names(mapping)) > 0) {
          for (i in length(names(mapping)):1) {
            tmp_map_val <- deparse(mapping[names(mapping)[i]][[1]])
            if (tmp_map_val[length(tmp_map_val)] %in% colnames(data)) 
              mapping[[names(mapping)[i]]] <- NULL
            if (length(names(mapping)) < 1) {
              mapping <- NULL
              break
            }
          }
        }
        if (length(colorCol) < 1) {
          colorCol <- "ggally_NO_EXIST"
        }
        if ((singleColorCol != "ggally_NO_EXIST") && (singleColorCol %in% 
                                                      colnames(data))) {
          cord <- as.data.frame(
            as.data.table(data)[, cor_fn(eval(parse(text=(xCol))), eval(parse(text=yCol))), 
                                  by=eval(colorCol)])
          
          colnames(cord)[2] <- "ggally_cor"
          cord$ggally_cor <- signif(as.numeric(cord$ggally_cor), 
                                    3)
          lev <- levels(data[[colorCol]])
          ord <- rep(-1, nrow(cord))
          for (i in 1:nrow(cord)) {
            for (j in seq_along(lev)) {
              if (identical(as.character(cord[i, colorCol]), 
                            as.character(lev[j]))) {
                ord[i] <- j
              }
            }
          }
          cord <- cord[order(ord[ord >= 0]), ]
          cord$label <- str_c(cord[[colorCol]], ": ", cord$ggally_cor)
          xmin <- min(xVal, na.rm = TRUE)
          xmax <- max(xVal, na.rm = TRUE)
          xrange <- c(xmin - 0.01 * (xmax - xmin), xmax + 0.01 * 
                        (xmax - xmin))
          ymin <- min(yVal, na.rm = TRUE)
          ymax <- max(yVal, na.rm = TRUE)
          yrange <- c(ymin - 0.01 * (ymax - ymin), ymax + 0.01 * 
                        (ymax - ymin))
          p <- ggally_text(label = str_c("Correlation: ", ""), 
                           mapping = mapping, xP = 0.5, yP = 0.9, 
                           xrange = xrange, yrange = yrange, color = "black", 
                           ...) + theme(legend.position = "none")
          xPos <- rep(alignPercent, nrow(cord)) * diff(xrange) + 
            min(xrange, na.rm = TRUE)
          yPos <- seq(from = 0.9, to = 0.2, length.out = nrow(cord) + 
                        1)
          yPos <- yPos * diff(yrange) + min(yrange, na.rm = TRUE)
          yPos <- yPos[-1]
          cordf <- data.frame(xPos = xPos, yPos = yPos, labelp = cord$label, 
                              colorCol = cord[[colorCol]])
          cordf$labelp <- factor(cordf$labelp, levels = cordf$labelp)
          p <- p + geom_text(data = cordf, 
                             aes(x = xPos, y = yPos, label = labelp, 
                                 color = colorCol), hjust = "left", ...) + 
            scale_color_manual(values = palette) + 
            theme(panel.grid = element_blank(), 
                  panel.grid.major = element_blank(), 
                  panel.background = element_rect(color="grey"))
            
          p
        } else {
          xmin <- min(xVal, na.rm = TRUE)
          xmax <- max(xVal, na.rm = TRUE)
          xrange <- c(xmin - 0.01 * (xmax - xmin), xmax + 0.01 * 
                        (xmax - xmin))
          ymin <- min(yVal, na.rm = TRUE)
          ymax <- max(yVal, na.rm = TRUE)
          yrange <- c(ymin - 0.01 * (ymax - ymin), ymax + 0.01 * 
                        (ymax - ymin))
          p <- ggally_text(
            label = paste("Correlation:\n", 
                          signif(cor_fn(xVal, yVal), 3), 
                          sep = "", collapse = ""), mapping, xP = 0.5, 
                          yP = 0.5, xrange = xrange, yrange = yrange, ...) + 
            theme(legend.position = "none")
          p
        }
      }
      
      densplot <- ggpairs(
        as.data.frame(dtm), 
        mapping = aes(fill = chain, color = chain), 
        columns = .availStats, 
        diag = list(continuous = my_dens),
        lower = list(continuous = my_points, 
                     combo = wrap("box_no_facet")),
        upper = list(continuous = my_cor, 
                     combo = wrap("box_no_facet")))
        
      if(doPlot) {
        print(traceplot) 
        if(interactive()) {
          print("Press Enter to see a univariate density plot")
          scan("", what = "character", nlines = 1)
        }
        print(densplot) 
      } else {
        list(traceplot = traceplot, densplot = densplot)  
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

