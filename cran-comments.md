## Test environments
* local OS X install, R 3.3.3 and R 3.3.2
* x86_64-slackware-linux-gnu, R 3.2.2
* x86_64-w64-mingw32 (64-bit) using devtools::build_win()

## R CMD check on Mac OS X:
There were no ERRORs, WARNINGs or NOTEs. 

## Downstream dependencies
Not applicable.

## ASAN and valgrind errors now should be fixed
The errors were due to read-only access to an address right after the last element of an Armadillo vector. Previously, there were no runtime crashes since the access was only read-only. Now this should be fixed.