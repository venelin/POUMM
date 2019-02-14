#' Cached objects for the POUMM vignettes and examples
#'
#' A list containing a simulated tree, trait-values and POUMM objects (model 
#' fits). To use these objects in examples you can load them into the global 
#' workspace with the command: 
#' `data(vignetteCachedResults); list2env(vignetteCachedResults, globalenv());`.
#'
#' @format This is a list containing the following named elements:
#' \describe{
#'   \item{g, z, e }{ numeric vectors of simulated genotypic values, phenotypic values and measurement errors. }
#'   \item{tree }{ a simulated phylogenetic tree. }
#'   \item{fitPOUMM, fitPOUMM2, fitH2tMean }{ POUMM fit objects to tree and z. }
#' }
"vignetteCachedResults"
