library(testthat)
library(POUMM)
library(mvtnorm)

EPS <- 10^-6

#' An "algebraic" implementation of the POUMM likelihood calculation based on
#' multivariate normal density function
#' This is an internal function used for tests only.
dVTipsGivenTreePOUMMg0Alg <- function(z, tree, alpha, theta, sigma, sigmae, g0) {
  tanc <- ape::vcv(tree)
  tipTimes <- diag(tanc)
  N <- length(tree$tip.label)
  
  D0 <- sapply(tipTimes, function(t) exp(-alpha*t))
  D1 <- 1-D0
  Vt <- covVTipsGivenTreePOUMM(tree, alpha, sigma, 0, tanc)
  Ve <- diag(sigmae^2, N, N)
  
  muz <- D0 * g0  +  D1 * theta
  
  mvtnorm::dmvnorm(z, muz, Vt+Ve, log=TRUE)
}

context("POUMM Likelihood")
set.seed(1)

N <- 100
tree <- ape::rtree(N)  

g0 <- 16
alpha <- 2
theta <- 4
sigma <- .2
sigmae <- .8
z <- rVNodesGivenTreePOUMM(tree, g0, alpha, theta, sigma, sigmae)

pruneInfo <- POUMM:::pruneTree(tree, z)

# pruneInfo2 <- POUMM::pruneTree(tree, z)
# 
# pruneInfo$integrator$reorderEdges()
# pruneInfo$integrator$edge

test_that(
  "fastLik(R) vs algebraic lik", {
    expect_true(EPS > abs(likPOUMMGivenTreeVTips(z[1:N], tree, 0, theta, sigma, sigmae, g0) -
                               dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, 0, theta, sigma, sigmae, g0)))
    expect_true(EPS > abs(likPOUMMGivenTreeVTips(z[1:N], tree, 0, theta, sigma, 0, g0) -
                               dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, 0, theta, sigma, 0, g0)))
    expect_true(EPS > abs(likPOUMMGivenTreeVTips(z[1:N], tree, alpha, theta, sigma, sigmae, g0) -
                               dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, alpha, theta, sigma, sigmae, g0)))
    expect_true(EPS > abs(likPOUMMGivenTreeVTips(z[1:N], tree, alpha, theta, sigma, 0, g0) -
                               dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, alpha, theta, sigma, 0, g0)))
  })

test_that(
  "fastLik(C++) vs algebraic lik", {
    expect_true(EPS > abs(likPOUMMGivenTreeVTipsC(pruneInfo$integrator, 0, theta, sigma, sigmae, g0) -
                               dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, 0, theta, sigma, sigmae, g0)))
    expect_true(EPS > abs(likPOUMMGivenTreeVTipsC(pruneInfo$integrator, 0, theta, sigma, 0, g0)  -
                               dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, 0, theta, sigma, 0, g0)))
    expect_true(EPS > abs(likPOUMMGivenTreeVTipsC(pruneInfo$integrator, alpha, theta, sigma, sigmae, g0) -
                               dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, alpha, theta, sigma, sigmae, g0)))
    expect_true(EPS > abs(likPOUMMGivenTreeVTipsC(pruneInfo$integrator, alpha, theta, sigma, 0, g0) -
                               dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, alpha, theta, sigma, 0, g0)))
  })

context("POUMM parameters") 

# At POUMM stationary state (equilibrium, t=Inf)
H2 <- POUMM::H2(alpha = 0.75, sigma = 1, sigmae = 1, t = Inf)     # 0.4
alpha <- POUMM::alpha(H2 = H2, sigma = 1, sigmae = 1, t = Inf)    # 0.75
sigma <- POUMM::sigmaOU(H2 = H2, alpha = 0.75, sigmae = 1, t = Inf) # 1
sigmae <- POUMM::sigmae(H2 = H2, alpha = 0.75, sigma = 1, t = Inf) # 1

test_that("stationary state (t = Inf)",
          expect_equal(c(H2, alpha, sigma, sigmae), c(0.4, 0.75, 1, 1)))

# At finite time t = 0.2
H2 <- POUMM::H2(alpha = 0.75, sigma = 1, sigmae = 1, t = 0.2)     # 0.1473309
alpha <- POUMM::alpha(H2 = H2, sigma = 1, sigmae = 1, t = 0.2)    # 0.75
sigma <- POUMM::sigmaOU(H2 = H2, alpha = 0.75, sigmae = 1, t = 0.2) # 1
sigmae <- POUMM::sigmae(H2  =  H2, alpha = 0.75, sigma = 1, t = 0.2) # 1

test_that("at finite time (t = 0.2)", 
          expect_equal(c(H2, alpha, sigma, sigmae), c(
            POUMM::H2(alpha, sigma, sigmae, 0.2), 0.75, 1, 1)))


context("specifyPOUMM")
test_that("default calls no error", {
  expect_silent(specifyPOUMM())
  expect_silent(specifyPOUMM_ATH2tMeanSe())
  expect_silent(specifyPOUMM_ATH2tMeanSeG0())
  expect_silent(specifyPOUMM_ATS())
  expect_silent(specifyPOUMM_ATSG0())
  expect_silent(specifyPOUMM_ATSSeG0())
  expect_silent(specifyPMM())
  expect_silent(specifyPMM_H2tMeanSe())
  expect_silent(specifyPMM_H2tMeanSeG0())
  expect_silent(specifyPMM_SSeG0())
})


context("POUMM fit")

# for(i in 1:10) {
#   zOU <- POUMM::rVNodesGivenTreePOUMM(tree, g0, alpha, theta, sigma, sigmae)
#   zBM <- POUMM::rVNodesGivenTreePOUMM(tree, g0, 0, 0, sigma, sigmae)
#   
#   test_that("fitPOUMM_on_zBM calls no error:", {
#     expect_silent(fitPMM_on_zBM <<- POUMM::POUMM(zBM[1:N], tree, spec = POUMM::specifyPMM(zBM[1:N], tree), doMCMC = FALSE))
#     expect_output(fitPOUMM_on_zBM <<- POUMM::POUMM(zBM[1:N], tree, doMCMC = TRUE, spec = specifyPOUMM(zBM[1:N], tree, nSamplesMCMC = 2e4)))
#     expect_silent(fitPMM_on_zOU <<- POUMM::POUMM(zOU[1:N], tree, spec = POUMM::specifyPMM(zOU[1:N], tree), doMCMC = FALSE))
#     expect_output(fitPOUMM_on_zOU <<- POUMM::POUMM(zOU[1:N], tree, doMCMC = TRUE, spec = specifyPOUMM(zBM[1:N], tree, nSamplesMCMC = 2e4)))
#   })
#   print("==============================================")
#   print(paste("Test", i))
#   print(lmtest::lrtest(fitPMM_on_zBM, fitPOUMM_on_zBM))
#   print("----------------------------------------------")
#   print(lmtest::lrtest(fitPMM_on_zOU, fitPOUMM_on_zOU))
#   
# }
# 
