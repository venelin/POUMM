---
title: "NEWS about the POUMM R-package"
author: "Venelin Mitov"
date: "04 January, 2017"
output: html_document
---

# POUMM 2.0
* Using SPLiTTree for the pruning likelihood calculation.
* Automatic web-site generation using the R-package pkgdown.
* Project source-code published on github: https://github.com/venelin/POUMM.git.
* Written vignette "Interpreting the POUMM".

# POUMM 1.3.2
* Fixing test-errors on r-patched-solaris-x86 and notes on r-devel-linux-x86_64-fedora-clang.
* Fixed memory errors reported by ASAN and valgrind.
* Improved default values for the parameter sigmae.

# POUMM 1.3
* Parallel likelihood calculation: (tested on Linux using Intel complier v16.0.0).
* Added possibility to specify known measurement standard error for each tip of the tree. 
(see `?POUMM`).
* Added a bivariate densplot of the posterior samples displaying a pairs-plot, i.e. a 
matrix of plots with marginal densities at the diagonal, scatter plots on the
lower triangle and correlation on the upper (see `?plot.summary.POUMM`);
possibility to change the color-palette used in the plot with color-blindness
friendly default colors.
* New default parameter settings: now the ML-limits (parUpper and parLower) as well as the
priors for the parameters alpha and sigma are scaled according to the length of the tree (tMean)
and the empirical variance of the trait (zVar) (see ?specifyPOUMM). 
Please note that the default parameter settings may change in future releases. 
* Improved ML-fit: now the optimization is run 50 times starting from different points in the 
search region. The ML-fit is still kept deterministic, i.e. consecutive runs should find the same
optimum.
* Fixed a dependency on the adaptMCMC package (now added as an import).
* Update references.

# POUMM 1.2.1
* First version of the package released on CRAN.

