
<!-- README.md is generated from README.Rmd. Please edit that file -->
Why should I use it?
====================

The Phylogenetic Ornstein-Uhlenbeck Mixed Model (POUMM) allows to estimate the phylogenetic heritability of continuous traits, to test hypotheses of neutral evolution versus stabilizing selection, to quantify the strength of stabilizing selection, to estimate measurement error and to make predictions about the evolution of a phenotype and phenotypic variation in a population. The POUMM package provides an easy and efficient way to perform this variety of analyses on large macro-evolutionary or epidemic trees. It implements a fast-likelihood calculation algorithm that makes it possible to perform MCMC-sampling with millions of iterations within an hour. Then, it provides a number of standard generic functions such as logLik, plot, summary allowing to quickly assess the quality of the model fit and to answer relevant biological questions about the data.

How do I use it?
================

Here is a quick example on how to use the package on a simulated tree and data:

``` r
# number of tips in the tree
N <- 500
# a random tree (a phylo object)
tr <- ape::rtree(N)
# simulate trait values along the tree
z <- POUMM::rVNodesGivenTreePOUMM(tr, 0, 2, 3, 1, 1)[1:N]
# fitting the POUMM to the tree and data
fit <- POUMM::POUMM(z, tr, spec = POUMM::specifyPOUMM(nSamplesMCMC = 5e4))
# plotting the inferred parameters (MCMC trace, and posterior densities):
plot(fit)
# summary as a table
summary(fit)
# information criterions based on maximum likelihood
AIC(fit)
BIC(fit)
# parameters maximizing the likelihood:
coef(fit)
# log-likelihood 
logLik(fit)
# most-likely genotypic values at the tips
fitted(fit)
# residuals represent the difference between the genotypic values and the 
# observed values. 
plot(resid(fit))
abline(h=0)

# fit PMM to the same data and do a likelihood ratio test
fitPMM <- POUMM::POUMM(z, tr, spec = POUMM::specifyPMM(nSamplesMCMC = 5e4))
lmtest::lrtest(fitPMM, fit)
```

You can find more examples in the package vignette entitled "A User Guide to The POUMM R-package" and in the help-pages, i.e. ?POUMM, ?specifyPOUMM, ?summary.POUMM, ?plot.POUMM.

How do I get it?
================

Install the newest package version from CRAN.

``` r
install.packages("POUMM")
```

This will install all 3rd party dependencies with one exception: the package for high precision floating point arithmetics (Rmpfr). While the POUMM can be run without Rmpfr installed, installing this package is recommended to prevent numerical instability in some extreme cases, i.e. very small positive values for the POUMM parameters alpha, sigma, and sigmae. Currently, Rmpfr is listed as "Suggests" in POUMM's DESCRIPTION, since it's installation has been problematic on some systems (i.e. Linux).
