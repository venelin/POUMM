library(testthat)
library(POUMM)

EPS <- 10^-9

context("POUMM Likelihood")
set.seed(1)

N <- 2e2
tree <- ape::rtree(N)  

g0 <- 16
alpha <- 2
theta <- 4
sigma <- .2
sigmae <- .8
z <- rVNodesGivenTreePOUMM(tree, g0, alpha, theta, sigma, sigmae)

pruneInfo <- POUMM:::pruneTree(tree, z)

test_that(
  "fastLik(R) vs algebraic lik", {
    expect_true(EPS > abs(likPOUMMGivenTreeVTips(z[1:N], tree, 0, theta, sigma, sigmae, g0) -
                               POUMM:::dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, 0, theta, sigma, sigmae, g0)))
    expect_true(EPS > abs(likPOUMMGivenTreeVTips(z[1:N], tree, 0, theta, sigma, 0, g0) -
                               POUMM:::dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, 0, theta, sigma, 0, g0)))
    expect_true(EPS > abs(likPOUMMGivenTreeVTips(z[1:N], tree, alpha, theta, sigma, sigmae, g0) -
                               POUMM:::dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, alpha, theta, sigma, sigmae, g0)))
    expect_true(EPS > abs(likPOUMMGivenTreeVTips(z[1:N], tree, alpha, theta, sigma, 0, g0) -
                               POUMM:::dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, alpha, theta, sigma, 0, g0)))
  })

test_that(
  "fastLik(C++) vs algebraic lik", {
    expect_true(EPS > abs(POUMM:::likPOUMMGivenTreeVTipsC(pruneInfo$integrator, 0, theta, sigma, sigmae, g0) -
                               POUMM:::dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, 0, theta, sigma, sigmae, g0)))
    expect_true(EPS > abs(POUMM:::likPOUMMGivenTreeVTipsC(pruneInfo$integrator, 0, theta, sigma, 0, g0)  -
                               POUMM:::dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, 0, theta, sigma, 0, g0)))
    expect_true(EPS > abs(POUMM:::likPOUMMGivenTreeVTipsC(pruneInfo$integrator, alpha, theta, sigma, sigmae, g0) -
                               POUMM:::dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, alpha, theta, sigma, sigmae, g0)))
    expect_true(EPS > abs(POUMM:::likPOUMMGivenTreeVTipsC(pruneInfo$integrator, alpha, theta, sigma, 0, g0) -
                               POUMM:::dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, alpha, theta, sigma, 0, g0)))
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
  expect_silent(specifyPOUMM_ATSSeG0())
  expect_silent(specifyPMM())
  expect_silent(specifyPMM_H2tMeanSe())
  expect_silent(specifyPMM_H2tMeanSeG0())
  expect_silent(specifyPMM_SSeG0())
})
