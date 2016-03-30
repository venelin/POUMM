NULL

#' Find a maximum likelihood fit of the POUMM model
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
                     divideEdgesBy=1, control=list(), verbose=FALSE, ...) {
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
  
  pruneInfo=pruneTree(tree)
  
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
        if(verbose&(count%%100==0|log)) {
          cat('Count calls=', count, '; parMemo=', toString(round(parMemo, 6)), 
              '; valMemo=', valMemo, ';   ### par=', toString(round(params, 6)), 
              '; value=', val, '\n', sep='')
        }
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
    value <- -do.call(likVTreeOU, 
                      c(list(v, tree), as.list(parFixedNoAlpha), as.list(par), list(log=T, distgr=distgr), 
                        list(impl='R5', pruneInfo=pruneInfo), list(...)))
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
      #cat('optimFixedAlpha called on alpha=', alpha, ', parFixed=', toString(parFixed), ', parInit=', parInit, ', ...=', toString(list(...)), '\n')
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

