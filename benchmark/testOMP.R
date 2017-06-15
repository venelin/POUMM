library(microbenchmark)
library(diversitree)
library(geiger)
library(OUwie)
library(bayou)
library(POUMM)
library(data.table)

set.seed(1)
nTests <- 100 # number of executions in each microbenchmark
Ns <- 10^(1:3)

g0 <- 16
alpha <- 2
theta <- 4
sigma <- .2
sigmae <- .8

a <- rexp(nTests)
th <- runif(nTests, 0, 4)
s <- rexp(nTests)
se <- rexp(nTests)
se[1:(nTests/2)] <- sigmae
th[1:(nTests/2)] <- mean(z[1:N])

times <- values <- trees <- NULL


for(treeType in c("ultrametric", "non-ultrametric")) {
  for(N in Ns) {
    cat("Tree-type:", treeType, "; tree-size:", N, "\n")
    
    if(treeType == "ultrametric") {
      tree <- ape::rcoal(N) 
    } else {
      tree <- ape::rtree(N) 
    }
    
    z <- rVNodesGivenTreePOUMM(tree, g0, alpha, theta, sigma, sigmae)
    zNamed <- z[1:N]
    names(zNamed) <- tree$tip.label
    
    dataOUwie <- data.table(Genus_species = tree$tip.label, Reg = 1, X = z[1:N])
    
    lik.geiger.C <- geiger:::bm.lik(tree, z[1:N], SE = NA, model="OU")
    
    
    lik.divtree.R <- make.ou(tree, z[1:N], states.sd =  sigmae, with.optimum = TRUE, 
                             control = list(method = 'pruning', backend = 'R'))
    
    lik.divtree.C <- make.ou(tree, z[1:N], states.sd =  sigmae, with.optimum = TRUE, 
                             control = list(method = 'pruning', backend = 'C'))
    
    # I tried to test bayou but it generates an index-out-of-out error. 
    # it is hard to test since there are no code examples in the doc.
    # bayou.cache <- bayou:::.prepare.ou.univariate(tree, zNamed, 0)
    
    # bayou.lik(pars = c(1, 1,1,1,1), tree, zNamed, model = "OU")
    # 
    # bayou.mcmc(tree, zNamed, prior = function(pars) {1})
    # 
    # OU.lik(pars = c(1, 1,1,1), tree, zNamed, model = "OU")
    # 
    
    
    pruneInfo <- pruneTree(tree, z[1:N])
    
    res_OUwie <- res_geiger_C <- res_divtree_R <- res_divtree_C <- 
      res_POUMM_R <- res_POUMM_C_arma <- res_POUMM_C_vec <- res_POUMM_C_par <- 
      res_diversitree_C <- res_diversitree_R <- res_geiger_C <- rep(0.0, nTests)
    
    counter_OUwie <- counter_geiger_C <- counter_divtree_R <- counter_divtree_C <- 
      counter_POUMM_R <- 
      counter_POUMM_C_arma <- counter_POUMM_C_vec <- counter_POUMM_C_par <- 1
    
    
    mb <- microbenchmark(
      "geiger: C++" = {
        i <- counter_geiger_C
        counter_geiger_C <- counter_geiger_C + 1
        res_geiger_C[i] <- lik.geiger.C(c(alpha = a[i], sigsq = s[i]^2, SE = se[i]))
      },
      "diversitree: R" = {
        i <- counter_divtree_R
        counter_divtree_R <- counter_divtree_R + 1
        res_divtree_R[i] <- lik.divtree.R(c(s2=s[i]^2, alpha=a[i], theta=th[i]))
      },
      "diversitree: C++" = {
        i <- counter_divtree_C
        counter_divtree_C <- counter_divtree_C + 1
        res_divtree_C[i] <- lik.divtree.C(c(s2=s[i]^2, alpha=a[i], theta=th[i]))
      },
      # very slow: skipping
      # "OUwie: R" = {
      #   i <- counter_OUwie
      #   counter_OUwie
      #   res <- OUwie.fixed(phy = tree, data = dataOUwie, model=c("OU1"), simmap.tree=FALSE, scaleHeight=FALSE,
      #               clade=NULL, alpha=a[i], sigma.sq=s[i]^2 , theta=th[i], quiet=TRUE)
      #   res_OUwie[i] <- res$loglik
      # },
      "POUMM: R" = {
        i <- counter_POUMM_R
        counter_POUMM_R <- counter_POUMM_R + 1
        res_POUMM_R[i] <- likPOUMMGivenTreeVTips(
          z = z[1:N], tree = tree, 
          alpha = a[i], theta = th[i], sigma = s[i], sigmae = se[i], pruneInfo = pruneInfo)
      },
      "POUMM: C++, Armadillo" = {
        i <- counter_POUMM_C_arma
        counter_POUMM_C_arma <- counter_POUMM_C_arma + 1
        res_POUMM_C_arma[i] <- likPOUMMGivenTreeVTipsC(
          alpha = a[i], theta = th[i], sigma = s[i], sigmae = se[i], 
          integrator = pruneInfo$integrator)
      },
      "POUMM: C++" = {
        i <- counter_POUMM_C_vec
        counter_POUMM_C_vec <- counter_POUMM_C_vec + 1
        res_POUMM_C_vec[i] <- likPOUMMGivenTreeVTipsC2(
          alpha = a[i], theta = th[i], sigma = s[i], sigmae = se[i], 
          integrator = pruneInfo$integrator)
      },
      "POUMM: C++, omp" = {
        i <- counter_POUMM_C_par
        counter_POUMM_C_par <- counter_POUMM_C_par + 1
        res_POUMM_C_par[i] <- likPOUMMGivenTreeVTipsC4(
          alpha = a[i], theta = th[i], sigma = s[i], sigmae = se[i], 
          integrator = pruneInfo$integrator)
      },
      times = nTests
    )
    
    values <- 
      rbind(values, 
            data.table(treeType = treeType, N = N, testId = 1:nTests, 
                       alpha = a, theta = th, sigma = s, sigmae = se, 
                       "geiger: C++" = res_geiger_C, 
                       "diversitree: R" = res_divtree_R, 
                       "diversitree: C++" = res_divtree_C,  
                       "POUMM: R" = res_POUMM_R, 
                       "POUMM: C++, Armadillo" = res_POUMM_C_arma, 
                       "POUMM: C++" = res_POUMM_C_vec, 
                       "POUMM: C++, omp" = res_POUMM_C_par))
    times <- 
      rbind(times, 
            cbind(data.table(treeType = treeType, N = N, as.data.table(print(mb)))))
    trees <- c(trees, list(tree))
  } 
}

