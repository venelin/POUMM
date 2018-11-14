library(testthat)
library(POUMM)

# test disabled in release-versions (CRAN) (takes too long)
if(POUMMIsADevRelease(numVersionComponents = 3)) {
context("POUMM fit")
  
  set.seed(1)
  
  N <- 100
  tree <- ape::rtree(N)  
  
  g0 <- 16
  alpha <- 2
  theta <- 4
  sigma <- .2
  sigmae <- .7
  se = rep(0, N) #rexp(N, 1/.01)
  z <- rVNodesGivenTreePOUMM(tree, g0, alpha, theta, sigma, sqrt(sigmae^2+se^2))
  
  source("generate-data.R")
  fit1 <- POUMM(z = z[1:N], tree = tree, 
                spec = specifyPOUMM(z[1:N], tree, nSamplesMCMC = 10000), 
                verbose=TRUE)
  
  if(interactive()) {
    summary(fit1)
    summary(fit1, mode="long")
    summary(fit1, mode="expert")
    
    plot(fit1)
  }

}