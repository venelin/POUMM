## Test environments
* local OS X install, R 3.3.3 and R 3.3.2
* x86_64-slackware-linux-gnu, R 3.2.2
* x86_64-w64-mingw32 (64-bit) using devtools::build_win() (see note below)

For devtools::build_win() the packaging and installation succeeds, and the example-code
works but the check failed because the vignette-builder has been killed after 20 
minutes:
* checking re-building of vignette outputs ... ERROR
Check process probably crashed or hung up for 20 minutes ... killed

On my laptop (Mac Book pro Retina late 2013) the code in the vignettes takes around
30 minutes to executed. It is possible to speed-up the vignette generation by including 
a cache-file"UserGuideCache.RData" in the package bundle. However, this will increase 
the size of the source-bundle to nearly 5MB (current size is 500k). I'm not sure 
what CRAN policy treats this issue with long-running code in R-vignettes.

## R CMD check on Mac OS X:
There were no ERRORs or WARNINGs. 

There were NOTEs:
* checking DESCRIPTION meta-information ... NOTE
Package listed in more than one of Depends, Imports, Suggests, Enhances:
  ‘Rcpp’
Not sure here, so I'm following examples found on the web. 
  
* checking dependencies in R code ... NOTE
package 'methods' is used but not declared

I couldn't understand why this is happening but since it is only a note, 
I'm ignoring it.

* checking R code for possible problems ... NOTE
Undefined global functions or variables:
  ESS G.R. HPD HPD50 HPDLower HPDLowerFiltered HPDUpper
  HPDUpperFiltered J Mean N chain estML it mcs nChains samplePriorMCMC stat value
  
This comes from data.table[, value:=expression] operators in my code. It's not
really an issue.
  

## Downstream dependencies
Not applicable since this is the first release of the package.