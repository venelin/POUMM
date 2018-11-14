## Previous errors  

All previous errors that caused the Archiving of the package in February 2018 have been fixed. These were:

```
CRAN repository db overrides:
    X-CRAN-Comment: Archived on 2018-02-01 as check errors were not
  
    Errors included undeclared use of package(s) in vignettes and failing
      corrected despite reminders.
      tests, and using 'microbenchmark' unconditionally.
  
  Size of tarball: 5326045 bytes
```

## Test environments
### local OS X install, R 3.5.1: 
1 NOTE :

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Venelin Mitov <vmitov@gmail.com>’

* checking installed package size ... NOTE
  installed size is  6.6Mb
  sub-directories of 1Mb or more:
    libs   5.5Mb
```

### devtools::check_win_release()

1 Note

### devtools::check_win_devel() 

1 Note

### travis ci:

1 Note

Results available here: https://travis-ci.org/venelin/POUMM

### devtools::check_rhub()
- 1 NOTE on Windows server (same as above)
- 1 Note on Ubuntu (same as above)
- Error on Fedora Linux:
    
```
Build ID:	POUMM_2.1.1.03.tar.gz-33b7de80573a4ebbbce1780b21d85cd5
Platform:	Fedora Linux, R-devel, clang, gfortran
Submitted:	21 minutes 5.9 seconds ago
Build time:	20 minutes 44.1 seconds

* installing *source* package "POUMM" ...
** libs
/usr/bin/clang++ -std=gnu++11 -I"/opt/R-devel/lib64/R/include" -DNDEBUG  -I"/home/docker/R/Rcpp/include" -I/usr/local/include    -fpic  -g -O2  -c Rcpp.cpp -o Rcpp.o
/usr/bin/clang++ -std=gnu++11 -shared -L/opt/R-devel/lib64/R/lib -L/usr/local/lib64 -o POUMM.so Rcpp.o -L/opt/R-devel/lib64/R/lib -lR
installing to /home/docker/POUMM.Rcheck/POUMM/libs
** R
** inst
** preparing package for lazy loading
Error in loadNamespace(j <- i[[1L]], c(lib.loc, .libPaths()), versionCheck = vI[[j]]) : 
  there is no package called "lattice"
ERROR: lazy loading failed for package "POUMM"
* removing "/home/docker/POUMM.Rcheck/POUMM"
```

- Error on Debian Linux:
    
```
Build ID:	POUMM_2.1.1.03.tar.gz-ff431dbd2ae84f18bb9ccc13fc26bb71
Platform:	Debian Linux, R-devel, GCC ASAN/UBSAN
Submitted:	24 minutes 13.5 seconds ago
Build time:	24 minutes 2 seconds

* 
About to run xvfb-run san.sh POUMM_2.1.1.03.tar.gz
* installing to library ‘/home/docker/R’
* installing *source* package ‘POUMM’ ...
** libs
g++ -fsanitize=undefined,bounds-strict -fno-omit-frame-pointer -std=gnu++11 -I"/usr/local/lib/R/include" -DNDEBUG  -I"/home/docker/R/Rcpp/include" -I/usr/local/include    -fpic  -g -O2 -Wall -pedantic -mtune=native -c Rcpp.cpp -o Rcpp.o
In file included from AbcPOUMM.h:28:0,
                 from Rcpp.cpp:3:
./SPLITT.h:2008:0: warning: ignoring #pragma omp simd [-Wunknown-pragmas]
     _PRAGMA_OMP_SIMD
 
./SPLITT.h:2020:0: warning: ignoring #pragma omp simd [-Wunknown-pragmas]
     _PRAGMA_OMP_SIMD
 
./SPLITT.h:2030:0: warning: ignoring #pragma omp simd [-Wunknown-pragmas]
     _PRAGMA_OMP_SIMD
 
./SPLITT.h:2039:0: warning: ignoring #pragma omp simd [-Wunknown-pragmas]
     _PRAGMA_OMP_SIMD
 
./SPLITT.h:2046:0: warning: ignoring #pragma omp simd [-Wunknown-pragmas]
       _PRAGMA_OMP_SIMD
 
./SPLITT.h:2069:0: warning: ignoring #pragma omp parallel [-Wunknown-pragmas]
 #pragma omp parallel
 
./SPLITT.h:2071:0: warning: ignoring #pragma omp for [-Wunknown-pragmas]
   _PRAGMA_OMP_FOR_SIMD
 
./SPLITT.h:2081:0: warning: ignoring #pragma omp barrier [-Wunknown-pragmas]
 #pragma omp barrier
 
./SPLITT.h:2084:0: warning: ignoring #pragma omp for [-Wunknown-pragmas]
     _PRAGMA_OMP_FOR_SIMD
 
./SPLITT.h:2094:0: warning: ignoring #pragma omp barrier [-Wunknown-pragmas]
 #pragma omp barrier
 
./SPLITT.h:2097:0: warning: ignoring #pragma omp for [-Wunknown-pragmas]
       _PRAGMA_OMP_FOR_SIMD
 
./SPLITT.h:2112:0: warning: ignoring #pragma omp parallel [-Wunknown-pragmas]
 #pragma omp parallel
 
./SPLITT.h:2121:0: warning: ignoring #pragma omp for [-Wunknown-pragmas]
   _PRAGMA_OMP_FOR_SIMD
 
./SPLITT.h:2130:0: warning: ignoring #pragma omp for [-Wunknown-pragmas]
     _PRAGMA_OMP_FOR_SIMD
 
./SPLITT.h:2155:0: warning: ignoring #pragma omp parallel [-Wunknown-pragmas]
 #pragma omp parallel
 
./SPLITT.h:2194:0: warning: ignoring #pragma omp parallel [-Wunknown-pragmas]
 #pragma omp parallel
 
./SPLITT.h:2196:0: warning: ignoring #pragma omp for [-Wunknown-pragmas]
   _PRAGMA_OMP_FOR_SIMD
 
./SPLITT.h:2206:0: warning: ignoring #pragma omp for [-Wunknown-pragmas]
     _PRAGMA_OMP_FOR_SIMD
 
./SPLITT.h:2219:0: warning: ignoring #pragma omp parallel [-Wunknown-pragmas]
 #pragma omp parallel
 
./SPLITT.h:2221:0: warning: ignoring #pragma omp for [-Wunknown-pragmas]
   _PRAGMA_OMP_FOR_SIMD
 
./SPLITT.h:2229:0: warning: ignoring #pragma omp for [-Wunknown-pragmas]
     _PRAGMA_OMP_FOR_SIMD
 
./SPLITT.h:2240:0: warning: ignoring #pragma omp parallel [-Wunknown-pragmas]
 #pragma omp parallel
 
./SPLITT.h:2249:0: warning: ignoring #pragma omp for [-Wunknown-pragmas]
   _PRAGMA_OMP_FOR_SIMD
 
./SPLITT.h:2259:0: warning: ignoring #pragma omp barrier [-Wunknown-pragmas]
 #pragma omp barrier
 
./SPLITT.h:2262:0: warning: ignoring #pragma omp for [-Wunknown-pragmas]
       _PRAGMA_OMP_FOR_SIMD
 
./SPLITT.h:2270:0: warning: ignoring #pragma omp simd [-Wunknown-pragmas]
       _PRAGMA_OMP_SIMD
 
./SPLITT.h:2283:0: warning: ignoring #pragma omp simd [-Wunknown-pragmas]
         _PRAGMA_OMP_SIMD
 
./SPLITT.h:2300:0: warning: ignoring #pragma omp parallel [-Wunknown-pragmas]
 #pragma omp parallel
 
./SPLITT.h:2309:0: warning: ignoring #pragma omp for [-Wunknown-pragmas]
   _PRAGMA_OMP_FOR_SIMD
 
./SPLITT.h:2319:0: warning: ignoring #pragma omp barrier [-Wunknown-pragmas]
 #pragma omp barrier
 
./SPLITT.h:2322:0: warning: ignoring #pragma omp for [-Wunknown-pragmas]
         _PRAGMA_OMP_FOR_SIMD
 
./SPLITT.h:2331:0: warning: ignoring #pragma omp simd [-Wunknown-pragmas]
         _PRAGMA_OMP_SIMD
 
./SPLITT.h:2345:0: warning: ignoring #pragma omp parallel [-Wunknown-pragmas]
 #pragma omp parallel
 
./SPLITT.h:2354:0: warning: ignoring #pragma omp for [-Wunknown-pragmas]
   _PRAGMA_OMP_FOR_SIMD
 
./SPLITT.h:2363:0: warning: ignoring #pragma omp barrier [-Wunknown-pragmas]
 #pragma omp barrier
 
./SPLITT.h:2366:0: warning: ignoring #pragma omp for [-Wunknown-pragmas]
       _PRAGMA_OMP_FOR_SIMD
 
./SPLITT.h:2383:0: warning: ignoring #pragma omp simd [-Wunknown-pragmas]
       _PRAGMA_OMP_SIMD
 
./SPLITT.h:2525:0: warning: ignoring #pragma omp simd [-Wunknown-pragmas]
     _PRAGMA_OMP_SIMD
 
./SPLITT.h:2539:0: warning: ignoring #pragma omp simd [-Wunknown-pragmas]
     _PRAGMA_OMP_SIMD
 
./SPLITT.h:2548:0: warning: ignoring #pragma omp simd [-Wunknown-pragmas]
       _PRAGMA_OMP_SIMD
 
./SPLITT.h:2556:0: warning: ignoring #pragma omp parallel [-Wunknown-pragmas]
 #pragma omp parallel
 
./SPLITT.h:2565:0: warning: ignoring #pragma omp for [-Wunknown-pragmas]
   _PRAGMA_OMP_FOR_SIMD
 
./SPLITT.h:2574:0: warning: ignoring #pragma omp for [-Wunknown-pragmas]
       _PRAGMA_OMP_FOR_SIMD
 
./SPLITT.h: In instantiation of ‘void SPLITT::PostOrderTraversal<TraversalSpecification>::TraverseTreeMultiThreadLoopVisits() [with TraversalSpecification = SPLITT::AbcPOUMM<SPLITT::OrderedTree<unsigned int, double> >]’:
./SPLITT.h:1885:79:   required from ‘void SPLITT::PostOrderTraversal<TraversalSpecification>::TraverseTree(SPLITT::PostOrderTraversal<TraversalSpecification>::ModeType) [with TraversalSpecification = SPLITT::AbcPOUMM<SPLITT::OrderedTree<unsigned int, double> >; SPLITT::PostOrderTraversal<TraversalSpecification>::ModeType = SPLITT::PostOrderMode]’
./SPLITT.h:583:5:   required from ‘SPLITT::TraversalTask<TraversalSpecification>::StateType SPLITT::TraversalTask<TraversalSpecification>::TraverseTree(const ParameterType&, SPLITT::uint) [with TraversalSpecification = SPLITT::AbcPOUMM<SPLITT::OrderedTree<unsigned int, double> >; SPLITT::TraversalTask<TraversalSpecification>::StateType = std::vector<double>; SPLITT::TraversalTask<TraversalSpecification>::ParameterType = std::vector<double>; SPLITT::uint = unsigned int]’
Rcpp.cpp:86:51:   required from here
./SPLITT.h:2114:8: warning: variable ‘tid’ set but not used [-Wunused-but-set-variable]
   uint tid;
        ^~~
g++ -fsanitize=undefined,bounds-strict -fno-omit-frame-pointer -std=gnu++11 -shared -L/usr/local/lib/R/lib -L/usr/local/lib -o POUMM.so Rcpp.o -L/usr/local/lib/R/lib -lR
installing to /home/docker/R/POUMM/libs
** R
** inst
** preparing package for lazy loading
Warning: S3 methods ‘network::as.network.phylo’, ‘igraph::as.igraph.phylo’, ‘network::as.network.evonet’, ‘igraph::as.igraph.evonet’ were declared in NAMESPACE but not found
Warning: S3 method ‘xts::as.xts.data.table’ was declared in NAMESPACE but not found
../../src/tbbmalloc/backref.cpp:154:21: runtime error: index 1 out of bounds for type 'BackRefBlock *[1]'
../../src/tbb/scheduler.cpp:326:23: runtime error: member call on address 0x7f0addb8ff40 which does not point to an object of type 'task'
0x7f0addb8ff40: note: object has invalid vptr
 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              invalid vptr
../../src/tbb/scheduler.cpp:354:32: runtime error: member call on address 0x7f0addb8ff40 which does not point to an object of type 'task'
0x7f0addb8ff40: note: object has invalid vptr
 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              invalid vptr
../../src/tbb/scheduler.cpp:1166:14: runtime error: member call on address 0x7f0addb8ff40 which does not point to an object of type 'task'
0x7f0addb8ff40: note: object has invalid vptr
 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              invalid vptr
../../src/tbb/scheduler.cpp:1168:14: runtime error: member call on address 0x7f0addb8ff40 which does not point to an object of type 'task'
0x7f0addb8ff40: note: object has invalid vptr
 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              invalid vptr
../../src/tbb/scheduler.h:560:34: runtime error: member call on address 0x7f0addb8ff40 which does not point to an object of type 'task'
0x7f0addb8ff40: note: object has invalid vptr
 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              invalid vptr
../../src/tbb/scheduler.h:560:34: runtime error: member call on address 0x7f0addb8ff40 which does not point to an object of type 'task'
0x7f0addb8ff40: note: object has invalid vptr
 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              invalid vptr
../../src/tbb/scheduler.h:560:34: runtime error: member call on address 0x7f0addb8ff40 which does not point to an object of type 'task'
0x7f0addb8ff40: note: object has invalid vptr
 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              invalid vptr
** help
*** installing help indices
** building package indices
** installing vignettes
** testing if installed package can be loaded
Warning: S3 methods ‘network::as.network.phylo’, ‘igraph::as.igraph.phylo’, ‘network::as.network.evonet’, ‘igraph::as.igraph.evonet’ were declared in NAMESPACE but not found
Warning: S3 method ‘xts::as.xts.data.table’ was declared in NAMESPACE but not found
../../src/tbbmalloc/backref.cpp:154:21: runtime error: index 1 out of bounds for type 'BackRefBlock *[1]'
../../src/tbb/scheduler.cpp:326:23: runtime error: member call on address 0x7f6767537f40 which does not point to an object of type 'task'
0x7f6767537f40: note: object has invalid vptr
 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              invalid vptr
../../src/tbb/scheduler.cpp:354:32: runtime error: member call on address 0x7f6767537f40 which does not point to an object of type 'task'
0x7f6767537f40: note: object has invalid vptr
 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              invalid vptr
../../src/tbb/scheduler.cpp:1166:14: runtime error: member call on address 0x7f6767537f40 which does not point to an object of type 'task'
0x7f6767537f40: note: object has invalid vptr
 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              invalid vptr
../../src/tbb/scheduler.cpp:1168:14: runtime error: member call on address 0x7f6767537f40 which does not point to an object of type 'task'
0x7f6767537f40: note: object has invalid vptr
 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              invalid vptr
../../src/tbb/scheduler.h:560:34: runtime error: member call on address 0x7f6767537f40 which does not point to an object of type 'task'
0x7f6767537f40: note: object has invalid vptr
 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              invalid vptr
../../src/tbb/scheduler.h:560:34: runtime error: member call on address 0x7f6767537f40 which does not point to an object of type 'task'
0x7f6767537f40: note: object has invalid vptr
 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              invalid vptr
../../src/tbb/scheduler.h:560:34: runtime error: member call on address 0x7f6767537f40 which does not point to an object of type 'task'
0x7f6767537f40: note: object has invalid vptr
 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              invalid vptr
* DONE (POUMM)
Running tests

R Under development (unstable) (2018-06-20 r74924) -- "Unsuffered Consequences"
Copyright (C) 2018 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> tools:::.runPackageTestsR()

  Running ‘test-all.R’

R Under development (unstable) (2018-06-20 r74924) -- "Unsuffered Consequences"
Copyright (C) 2018 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> library(testthat)
> 
> test_check("POUMM")
Loading required package: POUMM
Loading required package: Rcpp
Registered S3 methods overwritten by 'ggplot2':
  method         from 
  [.quosures     rlang
  c.quosures     rlang
  print.quosures rlang
../../src/tbbmalloc/backref.cpp:154:21: runtime error: index 1 out of bounds for type 'BackRefBlock *[1]'
../../src/tbb/scheduler.cpp:326:23: runtime error: member call on address 0x7f18aae5ff40 which does not point to an object of type 'task'
0x7f18aae5ff40: note: object has invalid vptr
 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              invalid vptr
../../src/tbb/scheduler.cpp:354:32: runtime error: member call on address 0x7f18aae5ff40 which does not point to an object of type 'task'
0x7f18aae5ff40: note: object has invalid vptr
 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              invalid vptr
../../src/tbb/scheduler.cpp:1166:14: runtime error: member call on address 0x7f18aae5ff40 which does not point to an object of type 'task'
0x7f18aae5ff40: note: object has invalid vptr
 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              invalid vptr
../../src/tbb/scheduler.cpp:1168:14: runtime error: member call on address 0x7f18aae5ff40 which does not point to an object of type 'task'
0x7f18aae5ff40: note: object has invalid vptr
 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              invalid vptr
../../src/tbb/scheduler.h:560:34: runtime error: member call on address 0x7f18aae5ff40 which does not point to an object of type 'task'
0x7f18aae5ff40: note: object has invalid vptr
 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              invalid vptr
../../src/tbb/scheduler.h:560:34: runtime error: member call on address 0x7f18aae5ff40 which does not point to an object of type 'task'
0x7f18aae5ff40: note: object has invalid vptr
 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              invalid vptr
../../src/tbb/scheduler.h:560:34: runtime error: member call on address 0x7f18aae5ff40 which does not point to an object of type 'task'
0x7f18aae5ff40: note: object has invalid vptr
 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              invalid vptr

Attaching package: 'POUMM'

The following object is masked from 'package:stats':

    simulate

  generate 10000 samples 
  generate 10000 samples 
  generate 10000 samples 
../../src/tbb/scheduler.h:560:34: runtime error: member call on address 0x7f18aae5ff40 which does not point to an object of type 'task'
0x7f18aae5ff40: note: object has invalid vptr
 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              invalid vptr
/home/docker/R/RcppParallel/include/tbb/task.h:749:49: runtime error: member call on address 0x7f18aae6a600 which does not point to an object of type 'scheduler'
0x7f18aae6a600: note: object is of type 'tbb::internal::custom_scheduler<tbb::internal::IntelSchedulerTraits>'
 00 00 00 00  b0 09 6a ab 18 7f 00 00  00 00 00 00 00 00 00 00  60 f6 e6 aa 18 7f 00 00  60 f6 e6 aa
              ^~~~~~~~~~~~~~~~~~~~~~~
              vptr for 'tbb::internal::custom_scheduler<tbb::internal::IntelSchedulerTraits>'
../../src/tbb/custom_scheduler.h:405:64: runtime error: member call on address 0x7f18aae5ff40 which does not point to an object of type 'task'
0x7f18aae5ff40: note: object has invalid vptr
 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              invalid vptr
══ testthat results  ═══════════════════════════════════════════════════════════
OK: 31 SKIPPED: 0 FAILED: 0
Warning messages:
1: S3 methods 'network::as.network.phylo', 'igraph::as.igraph.phylo', 'network::as.network.evonet', 'igraph::as.igraph.evonet' were declared in NAMESPACE but not found 
2: S3 method 'xts::as.xts.data.table' was declared in NAMESPACE but not found 
> 
> proc.time()
   user  system elapsed 
 68.076   0.180  68.345 
Running examples
> tools:::.createExdotR('POUMM', system.file(package = 'POUMM'))
  Extracting from parsed Rd's ...
> 
> 
Loading required package: Rcpp
Warning: S3 methods ‘network::as.network.phylo’, ‘igraph::as.igraph.phylo’, ‘network::as.network.evonet’, ‘igraph::as.igraph.evonet’ were declared in NAMESPACE but not found
Warning: S3 method ‘xts::as.xts.data.table’ was declared in NAMESPACE but not found
Registered S3 methods overwritten by 'ggplot2':
  method         from 
  [.quosures     rlang
  c.quosures     rlang
  print.quosures rlang
../../src/tbbmalloc/backref.cpp:154:21: runtime error: index 1 out of bounds for type 'BackRefBlock *[1]'
../../src/tbb/scheduler.cpp:326:23: runtime error: member call on address 0x7fc5d41cff40 which does not point to an object of type 'task'
0x7fc5d41cff40: note: object has invalid vptr
 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              invalid vptr
../../src/tbb/scheduler.cpp:354:32: runtime error: member call on address 0x7fc5d41cff40 which does not point to an object of type 'task'
0x7fc5d41cff40: note: object has invalid vptr
 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              invalid vptr
../../src/tbb/scheduler.cpp:1166:14: runtime error: member call on address 0x7fc5d41cff40 which does not point to an object of type 'task'
0x7fc5d41cff40: note: object has invalid vptr
 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              invalid vptr
../../src/tbb/scheduler.cpp:1168:14: runtime error: member call on address 0x7fc5d41cff40 which does not point to an object of type 'task'
0x7fc5d41cff40: note: object has invalid vptr
 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              invalid vptr
../../src/tbb/scheduler.h:560:34: runtime error: member call on address 0x7fc5d41cff40 which does not point to an object of type 'task'
0x7fc5d41cff40: note: object has invalid vptr
 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              invalid vptr
../../src/tbb/scheduler.h:560:34: runtime error: member call on address 0x7fc5d41cff40 which does not point to an object of type 'task'
0x7fc5d41cff40: note: object has invalid vptr
 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              invalid vptr
../../src/tbb/scheduler.h:560:34: runtime error: member call on address 0x7fc5d41cff40 which does not point to an object of type 'task'
0x7fc5d41cff40: note: object has invalid vptr
 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              invalid vptr

Attaching package: ‘POUMM’

The following object is masked from ‘package:stats’:

    simulate

[1] 0.2517647
[1] 0.25
../../src/tbb/scheduler.h:560:34: runtime error: member call on address 0x7fc5d41cff40 which does not point to an object of type 'task'
0x7fc5d41cff40: note: object has invalid vptr
 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              invalid vptr
/home/docker/R/RcppParallel/include/tbb/task.h:749:49: runtime error: member call on address 0x7fc5d41da600 which does not point to an object of type 'scheduler'
0x7fc5d41da600: note: object is of type 'tbb::internal::custom_scheduler<tbb::internal::IntelSchedulerTraits>'
 00 00 00 00  b0 f9 a0 d4 c5 7f 00 00  00 00 00 00 00 00 00 00  60 f6 1d d4 c5 7f 00 00  60 f6 1d d4
              ^~~~~~~~~~~~~~~~~~~~~~~
              vptr for 'tbb::internal::custom_scheduler<tbb::internal::IntelSchedulerTraits>'
../../src/tbb/custom_scheduler.h:405:64: runtime error: member call on address 0x7fc5d41cff40 which does not point to an object of type 'task'
0x7fc5d41cff40: note: object has invalid vptr
 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  00 00 00 00
              ^~~~~~~~~~~~~~~~~~~~~~~
              invalid vptr
Time elapsed:  3.04 0.068 3.112 0 0 
null device 
          1 
Running vignette code
> tools::buildVignettes(dir = '.', tangle = TRUE)
> 
> 
```

## Downstream dependencies
Not applicable.

