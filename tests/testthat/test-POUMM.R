library(testthat)
library(POUMM)
library(microbenchmark)

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

# pruneInfo <- POUMM:::pruneTree(tree)
# tmp <- tempfile()
# Rprof(tmp, interval=0.001, line.profiling=TRUE)
# for(i in 1:1000)
#   likPOUMMGivenTreeVTips(z[1:N], tree, 0, theta, sigma, sigmae, g0, pruneInfo = pruneInfo, backend = "C")
# Rprof(NULL)
# summaryRprof(tmp)

EPS <- 10^-9

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
    expect_true(EPS > abs(likPOUMMGivenTreeVTips(z[1:N], tree, 0, theta, sigma, sigmae, g0, backend = "C") -
                               POUMM:::dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, 0, theta, sigma, sigmae, g0)))
    expect_true(EPS > abs(likPOUMMGivenTreeVTips(z[1:N], tree, 0, theta, sigma, 0, g0, backend = "C") -
                               POUMM:::dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, 0, theta, sigma, 0, g0)))
    expect_true(EPS > abs(likPOUMMGivenTreeVTips(z[1:N], tree, alpha, theta, sigma, sigmae, g0, backend = "C") -
                               POUMM:::dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, alpha, theta, sigma, sigmae, g0)))
    expect_true(EPS > abs(likPOUMMGivenTreeVTips(z[1:N], tree, alpha, theta, sigma, 0, g0, backend = "C") -
                               POUMM:::dVTipsGivenTreePOUMMg0Alg(z[1:N], tree, alpha, theta, sigma, 0, g0)))
  })



# test_that("alpha",
#           expect_true(table[stat=='alpha', t.test(est.ML-alpha)$p.value > .01]))
# test_that("theta",
#           expect_true(table[stat=='theta', t.test(est.ML-theta)$p.value > .01]))
# 
# test_that("sigma",
#           expect_true(table[stat=='sigma2',
#                             t.test(sqrt(est.ML)-sigma)$p.value > .01]))
# test_that("sigmae",
#           expect_true(table[stat=='sigmae2',
#                             t.test(sqrt(est.ML)-sigmae)$p.value > .01]))
# 
# test_that("H2tMean", expect_true(table[stat=='H2tMean', t.test(est.ML-H2tMean)$p.value > .01]))
# test_that("H2tInf", expect_true(table[stat=='H2tInf', t.test(est.ML-H2tInf)$p.value > .01]))
# 
# 
