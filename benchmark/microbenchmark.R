library(microbenchmark)
library(diversitree)
library(POUMM)

set.seed(1)

N <- 1e5
tree <- ape::rtree(N)  

g0 <- 16
alpha <- 2
theta <- 4
sigma <- .2
sigmae <- .8
z <- rVNodesGivenTreePOUMM(tree, g0, alpha, theta, sigma, sigmae)


ou.lik <- make.ou(tree, z[1:N], states.sd =  sigmae, with.optimum = TRUE, control = list(method = 'pruning', backend = 'C'))
pruneInfo <- POUMM:::pruneTree(tree, z[1:N])
pruneInfo$integrator$reorderEdges()

microbenchmark(
  likPOUMMGivenTreeVTipsC(pruneInfo$integrator, alpha, theta, sigma, sigmae, g0),
  likPOUMMGivenTreeVTipsC2(pruneInfo$integrator, alpha, theta, sigma, sigmae, g0),
  likPOUMMGivenTreeVTipsC_old(pruneInfo$integrator, alpha, theta, sigma, sigmae, g0),
  likPOUMMGivenTreeVTips(z[1:N], tree, alpha, theta, sigma, sigmae, g0, pruneInfo = pruneInfo),
  ou.lik(c(s2=sigma^2, alpha=alpha, theta=theta)),
  times = 100)

Rprof()
for(i in 1:1000) POUMM:::likPOUMMGivenTreeVTipsC(pruneInfo$integrator, alpha, theta, sigma, sigmae, g0)
Rprof(NULL)
summaryRprof()


spec <- specifyPOUMM(
  z[1:N], tree, 
  parLower=c(alpha = 0.1, theta = 0.1, sigma = 0.1, sigmae = 0.1),
  parPriorMCMC = function(par) {
    dunif(par[1], 0.1, 10, TRUE) + dnorm(par[2], 4, 5, TRUE) + dunif(par[3], 0.1, 10, TRUE) + dunif(par[4], 0.1, 10, TRUE)
  },
  parInitMCMC = function(chainNo, fitML = NULL) {
    init <- rbind(c(alpha = 0.2, theta = 4, sigma = 1, sigmae = 0.2), 
                  c(alpha = 0.2, theta = 3, sigma = 1, sigmae = 1))
    init[(chainNo - 1)%%nrow(init) + 1, ]
  },
  parallelMCMC = FALSE, nSamplesMCMC = 1e5, samplePriorMCMC = FALSE, nChainsMCMC = 1) 

microbenchmark(
  fit1 <- POUMM(
    z[1:N], tree, 
    spec = spec,
    usempfr = 0, useArma = FALSE), 
  fit2 <- POUMM(
    z[1:N], tree, 
    spec = spec, 
    usempfr = 0),
  times = 1)

a <- rnorm(1e5)
b <- abs(a)
c <- a+b
print(c[1:5])
A <- matrix(rnorm(1e4), 1e2, 1e2)
B <- t(2.0* abs(A))
microbenchmark(a+b, a*b, log(b), a^2, A+B, A*B, log(B), times=100)
