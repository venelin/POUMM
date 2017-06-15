library(microbenchmark)
library(diversitree)
library(geiger)
library(POUMM)
library(data.table)
library(phylolm)


args <- commandArgs(trailingOnly = TRUE)

compilerInfo <- as.character(args[1])
nCores <- as.integer(args[2])
tryNo <- as.integer(args[3])
useCachedTreesAndValues <- as.logical(args[4])

if(R.version[['os']]=='linux-gnu') {
  # this only works on linux
  cpuInfo <- system("cat /proc/cpuinfo | grep 'model name' | uniq", intern = TRUE)
} else {
  # this only works on mac OS x
  cpuInfo <- system("sysctl -a -n machdep.cpu.brand_string", intern = TRUE)
}

print(cpuInfo)

set.seed(1)
nTests <- 100 # number of executions in each microbenchmark
Ns <- 10^(1:5)

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

times <- values <- trees <- NULL

treeNo <- 0

if(!useCachedTreesAndValues) {
  trees <- list()
  zs <- list()
} else {
  # load the trees and zs lists lists
  load('TreesAndValues.RData')
}

# perform the benchmark on two ultrametric and two non-ultrametric trees.
for(treeType in c("ultrametric", "ultrametric", "non-ultrametric", "non-ultrametric")) {
  for(N in Ns) {
    treeNo <- treeNo + 1
    
    cat("Tree-type:", treeType, "; tree-size:", N, "\n")
    
    if(!useCachedTreesAndValues) {
      cat("Generating tree...\n")
      if(treeType == "ultrametric") {
        tree <- ape::rcoal(N) 
      } else {
        tree <- ape::rtree(N) 
      }
      cat("Generating trait values...\n")
      z <- rVNodesGivenTreePOUMM(tree, g0, alpha, theta, sigma, sigmae)
      
      trees[[treeNo]] <- tree
      zs[[treeNo]] <- z
      save(trees, zs, file='TreesAndValues.RData')
    } else {
      cat("Using cached tree and trait-values...\n")
      tree <- trees[[treeNo]]
      z <- zs[[treeNo]]
    }
    
    th[1:(nTests/2)] <- mean(z[1:N])
    
    
    # perform other package execution only for single CPU and only 
    # in compiler == "icpc" case. All other cases are unnecessary duplicates 
    # since they don't influence the geiger and the diversitree implementations.
    if(nCores == 1 & compilerInfo == "icpc-omp-for-simd") {
      cat("Creating cache structures for diversitree and geiger...\n")
      lik.geiger.C <- geiger:::bm.lik(tree, z[1:N], SE = NA, model="OU")
      
      
      lik.divtree.R <- 
        make.ou(tree, z[1:N], states.sd =  sigmae, with.optimum = TRUE, 
                control = list(method = 'pruning', backend = 'R'))
      
      lik.divtree.C <- 
        make.ou(tree, z[1:N], states.sd =  sigmae, with.optimum = TRUE, 
                control = list(method = 'pruning', backend = 'C'))
    }
    cat("Creating cache structures for POUMM...\n")
    pruneInfo <- pruneTree(tree, z[1:N])
    
    res_phylolm <- res_geiger_C <- res_divtree_R <- res_divtree_C <- 
      res_POUMM_R <- res_POUMM_C_arma <- res_POUMM_C_par <- 
      res_diversitree_C <- res_diversitree_R <- res_geiger_C <- rep(0.0, nTests)
    
    counter_phylolm <- counter_geiger_C <- counter_divtree_R <- counter_divtree_C <- 
      counter_POUMM_R <- counter_POUMM_C_arma <- counter_POUMM_C_par <- 1
    
    counter_gc <- 1;
    
    if(nCores == 1 & compilerInfo == "icpc-omp-for-simd") {
      mbSerial <- microbenchmark(
        "gc" = {
          if(counter_gc %% 10 == 0) {
            gc()
          }
          counter_gc <-  counter_gc + 1
        },
        "phylolm: C++" = {
          i <- counter_phylolm
          counter_phylolm <- counter_phylolm + 1
          res_phylolm[i] <- OU1d.loglik(
            trait = z[1:N], phy = tree, model = "OUfixedRoot", 
            parameters = list(ancestral.state = 0, 
                              alpha = a[i], sigma2 = s[i]^2, 
                              sigma2_error = se[i]^2, 
                              optimal.value = th[i]))
        },
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
        "POUMM: R" = {
          i <- counter_POUMM_R
          counter_POUMM_R <- counter_POUMM_R + 1
          res_POUMM_R[i] <- likPOUMMGivenTreeVTips(
            z = z[1:N], tree = tree, 
            alpha = a[i], theta = th[i], sigma = s[i], sigmae = se[i], 
            pruneInfo = pruneInfo)
        },
        "POUMM: C++, Armadillo" = {
          i <- counter_POUMM_C_arma
          counter_POUMM_C_arma <- counter_POUMM_C_arma + 1
          res_POUMM_C_arma[i] <- likPOUMMGivenTreeVTipsC(
            alpha = a[i], theta = th[i], sigma = s[i], sigmae = se[i], 
            integrator = pruneInfo$integrator)
        },
        times = nTests)
    } else {
      mbSerial <- NULL
    }
    
    counter_gc <- 1;
    
    mbParallel <- microbenchmark(
      "gc" = {
        if(counter_gc %% 10 == 0) {
          gc()
        }
        counter_gc <-  counter_gc + 1
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
    
    if(nCores == 1 & compilerInfo == "icpc-omp-for-simd") {
      values <- 
        rbind(values, 
              data.table(compilerInfo = compilerInfo, 
                         cpuInfo=cpuInfo, nCores = nCores, tryNo = tryNo,
                         treeType = treeType, treeNo = treeNo, N = N, 
                         testId = 1:nTests, 
                         alpha = a, theta = th, sigma = s, sigmae = se, 
                         "geiger: C++" = res_geiger_C, 
                         "diversitree: R" = res_divtree_R, 
                         "diversitree: C++" = res_divtree_C,  
                         "POUMM: R" = res_POUMM_R, 
                         "POUMM: C++, Armadillo" = res_POUMM_C_arma, 
                         "POUMM: C++, omp" = res_POUMM_C_par))
    } else {
      values <- 
        rbind(values, 
              data.table(compilerInfo = compilerInfo, 
                         cpuInfo=cpuInfo, nCores = nCores, tryNo = tryNo,
                         treeType = treeType, treeNo = treeNo, N = N, 
                         testId = 1:nTests, 
                         alpha = a, theta = th, sigma = s, sigmae = se, 
                         "geiger: C++" = as.double(NA), 
                         "diversitree: R" = as.double(NA), 
                         "diversitree: C++" = as.double(NA),  
                         "POUMM: R" = as.double(NA), 
                         "POUMM: C++, Armadillo" = as.double(NA), 
                         "POUMM: C++, omp" = res_POUMM_C_par))
    }
    
    times <- 
      rbind(times, 
            cbind(data.table(compilerInfo = compilerInfo, 
                             cpuInfo=cpuInfo, nCores = nCores, tryNo = tryNo,
                             treeType = treeType, N = N, 
                             as.data.table(print(mbParallel, unit="ms")))))
    if(nCores == 1 & compilerInfo == "icpc-omp-for-simd") {
      times <- 
        rbind(times, 
              cbind(data.table(compilerInfo = compilerInfo, 
                               cpuInfo=cpuInfo, nCores = nCores, tryNo = tryNo,
                               treeType = treeType, N = N, 
                               as.data.table(print(mbSerial, unit="ms")))))
    }
    trees <- c(trees, list(tree))
  } 
}

save(values, times, 
     file = paste0("Results_", compilerInfo, "_", nCores, "_cores_", tryNo, ".RData"))
