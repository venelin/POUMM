library(testthat)
library(POUMM)
library(mvtnorm)

EPS <- 10^-4

#' An "algebraic" implementation of the POUMM likelihood calculation based on
#' multivariate normal density function
#' This is an internal function used for tests only.
dVTipsGivenTreePOUMMg0Alg <- function(z, tree, alpha, theta, sigma, sigmae, se, g0) {
  tanc <- ape::vcv(tree)
  tipTimes <- diag(tanc)
  N <- length(tree$tip.label)
  
  D0 <- sapply(tipTimes, function(t) exp(-alpha*t))
  D1 <- 1-D0
  Vt <- covVTipsGivenTreePOUMM(tree, alpha, sigma, 0, tanc)
  Ve <- diag(sigmae^2+se^2, N, N)
  
  muz <- D0 * g0  +  D1 * theta
  
  mvtnorm::dmvnorm(z, muz, Vt+Ve, log=TRUE)
}

absCalls <- 0
my_abs <- function(x) {
  absCalls <<- absCalls+1
  cat("abs call #", absCalls, ": x=", x, "\n", sep="")
  abs(x)
}

context("POUMM Likelihood")
set.seed(1)

N <- 100
tree <- ape::rtree(N)  

g0 <- 16
alpha <- 2
theta <- 4
sigma <- .2
sigmae <- .7
se = rexp(N, 1/.01)
z <- rVNodesGivenTreePOUMM(tree, g0, alpha, theta, sigma, sqrt(sigmae^2+se^2))

pruneInfo <- POUMM:::pruneTree(tree, z, se)



test_that(
  "fastLik(R) vs algebraic lik", {
    expect_true(EPS > my_abs(likPOUMMGivenTreeVTips(z[1:N], tree, 0, theta, sigma, sqrt(sigmae^2+se^2), g0) -
                               dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, 0, theta, sigma, sigmae, se, g0)))
    expect_true(EPS > my_abs(likPOUMMGivenTreeVTips(z[1:N], tree, 0, theta, sigma, sqrt(0+se^2), g0) -
                               dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, 0, theta, sigma, 0, se, g0)))
    expect_true(EPS > my_abs(likPOUMMGivenTreeVTips(z[1:N], tree, alpha, theta, sigma, sqrt(sigmae^2+se^2), g0) -
                               dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, alpha, theta, sigma, sigmae, se, g0)))
    expect_true(EPS > my_abs(likPOUMMGivenTreeVTips(z[1:N], tree, alpha, theta, sigma, sqrt(0+se^2), g0) -
                               dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, alpha, theta, sigma, 0, se, g0)))
    expect_true(EPS > my_abs(likPOUMMGivenTreeVTips(z[1:N], tree, alpha, theta, sigma, sqrt(sigmae^2+0), g0) -
                            dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, alpha, theta, sigma, sigmae, 0, g0)))
  })

test_that(
  "fastLik(C++ Armadillo) vs algebraic lik", {
    expect_true(EPS > my_abs(likPOUMMGivenTreeVTipsC(pruneInfo$integrator, 0, theta, sigma, sigmae, g0) -
                               dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, 0, theta, sigma, sigmae, se, g0)))
    expect_true(EPS > my_abs(likPOUMMGivenTreeVTipsC(pruneInfo$integrator, 0, theta, sigma, 0, g0)  -
                               dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, 0, theta, sigma, 0, se, g0)))
    expect_true(EPS > my_abs(likPOUMMGivenTreeVTipsC(pruneInfo$integrator, alpha, theta, sigma, sigmae, g0) -
                               dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, alpha, theta, sigma, sigmae, se, g0)))
    expect_true(EPS > my_abs(likPOUMMGivenTreeVTipsC(pruneInfo$integrator, alpha, theta, sigma, 0, g0) -
                               dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, alpha, theta, sigma, 0, se, g0)))
    expect_true(EPS > my_abs(likPOUMMGivenTreeVTips(z[1:N], tree, alpha, theta, sigma, sigmae, g0) -
                            dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, alpha, theta, sigma, sigmae, 0, g0)))
  })

test_that(
  "fastLik(C++ vectorized) vs algebraic lik", {
    expect_true(EPS > my_abs(likPOUMMGivenTreeVTipsC2(pruneInfo$integrator, 0, theta, sigma, sigmae, g0) -
                            dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, 0, theta, sigma, sigmae, se, g0)))
    expect_true(EPS > my_abs(likPOUMMGivenTreeVTipsC2(pruneInfo$integrator, 0, theta, sigma, 0, g0)  -
                            dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, 0, theta, sigma, 0, se, g0)))
    expect_true(EPS > my_abs(likPOUMMGivenTreeVTipsC2(pruneInfo$integrator, alpha, theta, sigma, sigmae, g0) -
                            dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, alpha, theta, sigma, sigmae, se, g0)))
    expect_true(EPS > my_abs(likPOUMMGivenTreeVTipsC2(pruneInfo$integrator, alpha, theta, sigma, 0, g0) -
                            dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, alpha, theta, sigma, 0, se, g0)))
  })


test_that(
  "fastLik(C++ vectorized and multicore) vs algebraic lik", {
    expect_true(EPS > my_abs(likPOUMMGivenTreeVTipsC4(pruneInfo$integrator, 0, theta, sigma, sigmae, g0) -
                            dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, 0, theta, sigma, sigmae, se, g0)))
    expect_true(EPS > my_abs(likPOUMMGivenTreeVTipsC4(pruneInfo$integrator, 0, theta, sigma, 0, g0)  -
                            dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, 0, theta, sigma, 0, se, g0)))
    expect_true(EPS > my_abs(likPOUMMGivenTreeVTipsC4(pruneInfo$integrator, alpha, theta, sigma, sigmae, g0) -
                            dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, alpha, theta, sigma, sigmae, se, g0)))
    expect_true(EPS > my_abs(likPOUMMGivenTreeVTipsC4(pruneInfo$integrator, alpha, theta, sigma, 0, g0) -
                            dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, alpha, theta, sigma, 0, se, g0)))
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
fit1 <- POUMM(z = z[1:N], tree = tree, 
              spec = specifyPOUMM(z[1:N], tree, nSamplesMCMC = 10000), 
              verbose=TRUE)

context("summary POUMM fit")
summary(fit1)
summary(fit1, mode="long")
summary(fit1, mode="expert")

context("plot POUMM fit") 
plot(fit1)


context("POUMM Likelihood speed")
set.seed(1)

for(N2 in c(100, 1000, 10000)) {
  tree2 <- ape::rtree(N2)  
  
  se2 = rexp(N2, 1/.01)
  z2 <- rVNodesGivenTreePOUMM(tree2, g0, alpha, theta, sigma, sqrt(sigmae^2+se2^2))
  
  pruneInfo <- POUMM:::pruneTree(tree2, z2, se2)
  
  mb <- microbenchmark::microbenchmark(
    R=likPOUMMGivenTreeVTips(z2[1:N2], tree2, 0, theta, sigma, sqrt(sigmae^2+se2^2), g0),
    Armadillo=likPOUMMGivenTreeVTipsC(pruneInfo$integrator, 0, theta, sigma, sigmae, g0),
    SIMD=likPOUMMGivenTreeVTipsC2(pruneInfo$integrator, 0, theta, sigma, sigmae, g0),
    SIMD_OMP=likPOUMMGivenTreeVTipsC4(pruneInfo$integrator, 0, theta, sigma, sigmae, g0),
    times=10
  )
  
  cat("N2=")
  print(N2)
  print(mb)
}

