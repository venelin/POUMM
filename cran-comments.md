## Test environments
* local OS X install, R 3.3.2
* x86_64-slackware-linux-gnu, R 3.2.2

## R CMD check results
There were no ERRORs or WARNINGs. 

There were 1 NOTE:
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
  
* checking Rd line widths ... NOTE
Rd file 'plot.POUMM.Rd':
  \usage lines wider than 90 characters:
       zoomInFilter = "(stat %in% c('H2e','H2tMean','H2tInf','H2tMax') | (value >= HPDLower & value <= HPDUpper))",

Rd file 'plot.summary.POUMM.Rd':
  \usage lines wider than 90 characters:
       zoomInFilter = "(stat %in% c('H2e','H2tMean','H2tInf','H2tMax') | (value >= HPDLower & value <= HPDUpper))",


## Downstream dependencies
Not applicable since this is the first release of the package.