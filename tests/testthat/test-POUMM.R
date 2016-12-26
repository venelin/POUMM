library(testthat)
library(poumm)
library(TreeSim)
library(data.table)
library(foreach)
library(doParallel)

test_check("poumm")

set.seed(1)


# 1. We generate one non-ultrametric tree using the TreeSim package.
# 2. We fix the "true" inital genotypic value, g0, and "true" POUMM parameters.
# 3. We repeatedly simulate trait evolution on the tree under the specified true 
# POUMM parameters and use the POUMM function with doMCMC = FALSE to perform an 
# ML fit only.
# 4. For each POUMM parameter and for the phylogenetic heritability we use t.test
# to check that the sample mean of the ML-values obtained from the simulation
# is not significantly different from the true parameter value.

# 1.
# Number of tips
N <- 4000
tree <- TreeSim::sim.bdsky.stt(N, lambdasky = 2, deathsky = 1, 
                               timesky=c(0, Inf), sampprobsky = 1)[[1]]

tMean <- mean(nodeTimes(tree, tipsOnly = TRUE))

# 2. 
# true POUMM parameters
g0 <- 16
alpha <- 2
theta <- 8
sigma <- 2
sigmae <- 1

# "true" phylogenetic heritability at tMean
H2tMean <- H2(alpha, sigma, sigmae, tMean)

# "true" equilibrium phylogenetic heritability
H2tInf <- H2(alpha, sigma, sigmae, Inf)

# 3.
# number of cores to run ML fit on:
nCores <- 8

# number of simulations replicates
nSim <- 128

cluster <- makeCluster(8)
registerDoParallel(cluster)

simulPOUMM <- foreach(i=1:nSim, .packages = "poumm") %dopar% {
  set.seed(i)
  simul <- list()
  simul$g <- rVNodesGivenTreePOUMM(tree, z0 = g0, alpha, theta, sigma, sigmae = 0)
  simul$e <- rnorm(length(simul$g), 0, sigmae)
  simul$z <- simul$g + simul$e
  simul$fit <- POUMM(z=simul$z[1:N], tree=tree, distgr=g0, parUpper = c(alpha=8, theta=16, sigma = 8, sigmae = 4), doMCMC = FALSE,
                     verbose = TRUE)
  simul
}

stopCluster(cluster)

# 4.
table <- rbindlist(lapply(1:length(simulPOUMM), function(iSim) {
  tab <- summary(simulPOUMM[[iSim]]$fit, start.MCMC = 0)
  tab[, i:=iSim]
  tab[stat=="H2tMean", value.Sim:=var(simulPOUMM[[iSim]]$g) / var(simulPOUMM[[iSim]]$z)]
  tab
}))



test_that("alpha", 
          expect_true(table[stat=='alpha', t.test(est.ML-alpha)$p.value > .01]))
test_that("theta", 
          expect_true(table[stat=='theta', t.test(est.ML-theta)$p.value > .01]))

test_that("sigma", 
          expect_true(table[stat=='sigma2', 
                            t.test(sqrt(est.ML)-sigma)$p.value > .01]))
test_that("sigmae", 
          expect_true(table[stat=='sigmae2', 
                            t.test(sqrt(est.ML)-sigmae)$p.value > .01]))

test_that("H2tMean", expect_true(table[stat=='H2tMean', t.test(est.ML-H2tMean)$p.value > .05]))
test_that("H2tInf", expect_true(table[stat=='H2tInf', t.test(est.ML-H2tInf)$p.value > .05]))

# 
# phytools::plotBranchbyTrait(tree, z, mode = "nodes", 
#                             show.tip.label = FALSE, 
#                             show.node.label = FALSE)
# 
# H2e <- H2e(z = z[1:N], sigmae = sigmae)
# 
# # H2 at equilibrium
# c(H2, poumm::H2(alpha, sigma, sigmae, Inf), H2e)
# 
# fit <- POUMM(z=z[1:N], tree = tree, distgr = 'maxlik', verbose = TRUE, usempfr = -2, doMCMC=FALSE)
# 
# # genotypic values
# g2 <- rVNodesGivenTreePOUMM(tree2, z0 = 8, alpha, theta, sigma, sigmae = 0)
# 
# # adding random environmental deviation
# e2 <- rnorm(length(g2), 0, sigmae)
# 
# z2 <- g2 + e2
# 
# likPOUMMGivenTreeVTips(z2[1:N], tree, alpha, theta, sigma, sigmae, distgr='maxlik')
# 
# fit2 <- POUMM(z=z2[1:N], tree=tree, distgr='maxlik', verbose = TRUE, usempfr = FALSE, doMCMC = TRUE)
# 
# H2e <- H2e(z = z2[1:N], sigmae = sigmae)
# 
# x0 <- 8
# nSteps <- 10000
# dt <- 0.0001
# set.seed(1)
# traj <- rTrajCondOU(x0, dt, 2, 2, 1, steps=nSteps)
# set.seed(1)
# traj2 <- rTrajCondOUDef(x0, dt, 2, 2, 1, steps=nSteps)
# set.seed(1)
# traj3 <- rep(0, nSteps)
# traj3[1] <- x0
# for(i in 1:nSteps) {
#   traj3[i+1] <- rCondOU(1, traj3[i], dt, 2, 2, 1)
# }
# 
# plot(traj, type='l')
# lines(traj2, lty=2, col='red')
# lines(traj3[-1], lty=3, col='yellow')
# 
# all.equal(traj, traj2)
# all.equal(traj, traj3[-1])
# 
# traj-traj2[-1]
# 
# var(rCondOU(1e6, 8, Inf, 2, 2, 1))
