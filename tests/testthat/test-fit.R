library(testthat)
library(POUMM)

# test disabled in release-versions (CRAN) (takes too long)
if(POUMMIsADevRelease()) {
context("POUMM fit")
  
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  
  N <- 100
  tree <- ape::rtree(N)  
  
  g0 <- 16
  alpha <- 2
  theta <- 4
  sigma <- .2
  sigmae <- .7
  se = rep(0, N) #rexp(N, 1/.01)
  z <- rVNodesGivenTreePOUMM(tree, g0, alpha, theta, sigma, sqrt(sigmae^2+se^2))
  
  test_that("fit passes without errors", {
    expect_is(fit1 <- POUMM(z = z[1:N], tree = tree, 
                            spec = specifyPOUMM(z[1:N], tree, nSamplesMCMC = 10000), 
                            verbose=FALSE), "POUMM")
    
    # end value not changed warnings from coda
    expect_warning(smmShort <- summary(fit1))
    expect_is(smmShort, "summary.POUMM")
    
    expect_warning(smmLong <- summary(fit1, mode="long"))
    expect_is(smmLong, "summary.POUMM")
    
    expect_warning(smmExpert <- summary(fit1, mode="expert"))
    expect_is(smmExpert, "summary.POUMM")
    
    expect_warning(plList <- plot(fit1, interactive = FALSE, doPlot = FALSE))
    expect_is(plList, "list")
  })
  
}