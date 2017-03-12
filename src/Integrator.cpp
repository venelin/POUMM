#include <RcppArmadillo.h>
#include <algorithm>
#include <R_ext/Rdynload.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Needed to compile on WINDOWS:
typedef unsigned int uint;

// BEGIN: Needed for r-devel (R 3.4)
void R_init_POUMM(DllInfo *info) {
  /* Register routines,
  allocate resources. */
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}

void R_unload_POUMM(DllInfo *info) {
    /* Release resources. */
}
// END Needed for r-devel (R 3.4)

class Integrator {
  uvec ZERO;
  uvec ONE;
  uvec TWO;
  
  uint nLevels;
  uint M;
  uint N;
  
  arma::vec zReord;
  arma::vec tReord;
  arma::uvec eReord;
  
  arma::uvec nodesIndex;
  arma::uvec unIndex;
  
  arma::mat abcMat;
  arma::vec a, b, c;
  
  uint count_abc_calls;
  
public:
  // Default constructor; 
  Integrator(): ZERO(1), ONE(1), TWO(1) {
    ZERO.fill(0); ONE.fill(1); TWO.fill(2);
  };
  
  // should be called after constructing an Integrator object
  void setPruningInfo(
      arma::vec z_, arma::umat edge_, arma::vec t_,
      uint M_, uint N_,
      arma::uvec endingAt_,
      arma::uvec nodesVector_, 
      arma::uvec nodesIndex_,
      arma::uvec unVector_, 
      arma::uvec unIndex_) {
    
    count_abc_calls = 0;
    
    M = M_; N = N_; 
    
    // this is the actual value of nLevels (not corrected for 0 indexing)
    nLevels = nodesIndex_.n_elem - 1;
    
    a = b = c = vec(M);
    
    // values corrected for 0-based indexing
    nodesIndex = nodesIndex_ - 1; 
    unIndex = unIndex_ - 1; 
    
    reorderEdges(edge_ - 1, t_, z_, 
                 endingAt_ - 1, nodesVector_ - 1, nodesIndex, 
                 unVector_ - 1, unIndex, M, nLevels);
  };
  
  void multiReplace(arma::uvec &x, arma::uvec const& a, arma::uvec const& b) {
    arma::uvec ind = sort_index(x);
    arma::uvec xInd = x(ind);
    std::pair<arma::uvec::iterator, arma::uvec::iterator> bounds;
    for(int i = 0; i < a.n_elem; ++i) {
      bounds = std::equal_range(xInd.begin(), xInd.end(), a.at(i));
      if(bounds.first != bounds.second) {
        uint first = bounds.first - xInd.begin();
        uint last = bounds.second - xInd.begin() - 1;
        x.elem(ind.subvec(first, last)).fill(b.at(i));
      }
    }
  }
  
  void reorderEdges(arma::umat const& edge, arma::vec const& t, 
                    arma::vec const& z,
                    arma::uvec const& endingAt,
                    arma::uvec const& nodesVector,
                    arma::uvec const& nodesIndex, 
                    arma::uvec const& unVector,
                    arma::uvec const& unIndex,
                    uint M, uint nLevels) {
    
    arma::uvec parents = edge.col(0);
    parents.replace(N, 2*M-1);
    eReord = parents;
    tReord = t;
    zReord = z;
    
    // edges pointing to tips
    uvec es = endingAt.elem(nodesVector(span(nodesIndex(0) + 1, nodesIndex(1))));
    
    uint unJ = 0;
    
    uint lenUnAll = 0;
    while(lenUnAll != es.n_elem) {
      uvec un = unVector(span(unIndex(unJ) + 1, unIndex(unJ + 1)));

      uvec edgeEnds = edge(es(un), ONE);
      uvec edgeEndsNew = arma::regspace<uvec>(unIndex(unJ) + 1, unIndex(unJ + 1));
      multiReplace(parents, edgeEnds, M+edgeEndsNew);
      
      eReord.subvec(unIndex(unJ) + 1, unIndex(unJ + 1)) = parents(es(un));
      multiReplace(eReord, edgeEnds, M+edgeEndsNew);
      
      tReord.subvec(unIndex(unJ) + 1, unIndex(unJ + 1)) = t.elem(es(un));
      zReord.subvec(unIndex(unJ) + 1, unIndex(unJ + 1)) = z.elem(edgeEnds);
      
      ++unJ;
      lenUnAll += un.n_elem;
    }

    // edges pointing to internal nodes
    for(int i = 1; i < nLevels; ++i) {
      es = endingAt.elem(nodesVector(span(nodesIndex(i) + 1, nodesIndex(i + 1))));
      
      uint lenUnAll = 0;
      while(lenUnAll != es.n_elem) {
        uvec un = unVector(span(unIndex(unJ) + 1, unIndex(unJ + 1)));

        uvec edgeEnds = edge(es(un), ONE);
        uvec edgeEndsNew = arma::regspace<uvec>(unIndex(unJ) + 1, unIndex(unJ + 1));
        multiReplace(parents, edgeEnds, M+edgeEndsNew);
        
        eReord.subvec(unIndex(unJ) + 1, unIndex(unJ + 1)) = parents(es(un));
        multiReplace(eReord, edgeEnds, M+edgeEndsNew);
        
        tReord.subvec(unIndex(unJ) + 1, unIndex(unJ + 1)) = t.elem(es(un));
        
        ++unJ;
        lenUnAll += un.n_elem;
      }
    }
    
    eReord = eReord - M;
  };
  
  arma::vec abc(double alpha, double theta, double sigma, double sigmae) {
    
    ++count_abc_calls;
    
    // time for main-loop
    //auto start = std::chrono::high_resolution_clock::now();
    
    double sigma2 = sigma*sigma, logsigma = log(sigma), 
      logsigmae = log(sigmae), sigmae2 = sigmae*sigmae;
    
    arma::vec talpha = tReord * alpha;
    arma::vec etalpha = exp(talpha);
    arma::vec e2talpha = etalpha % etalpha;
    arma::vec fe2talpha(M);
    if(alpha != 0) {
      fe2talpha = alpha / (1 - e2talpha);
    } else {
      fe2talpha = -0.5 / tReord;
    }
    
    // needed to re-initialize abcMat;
    a.fill(0); b.fill(0); c.fill(0);
    
    // edges pointing to tips
    uint eFirst = nodesIndex(0) + 1, eLast = nodesIndex(1);
    
    vec gutalphasigma2(M);
    
    gutalphasigma2(span(eFirst, eLast)) = e2talpha(span(eFirst, eLast)) + 
      ((-0.5 / sigmae2) * sigma2) / fe2talpha(span(eFirst, eLast));
    
    // (possibly reordered) shifted tip values
    vec z1 = zReord(span(eFirst, eLast)) - theta ;
    
    vec z1z1 = z1 % z1;
    
    if(sigmae != 0) {
      // integration over g1 including e1 = z1 - g1
      c(span(eFirst, eLast)) = -0.5 * log(gutalphasigma2(span(eFirst, eLast))) -
        //0.25 * sigma2 * pow(z1 / sigmae2, 2) /
        0.25 * sigma2 * z1z1 / (sigmae2*sigmae2) /
          (fe2talpha(span(eFirst, eLast)) - alpha + (-0.5 / sigmae2) * sigma2) +
            //talpha(span(eFirst, eLast)) + (-0.5 * (M_LN_2PI  + z1 % (z1 / sigmae2)) - logsigmae) ;
      talpha(span(eFirst, eLast)) + (-0.5 * (M_LN_2PI  + z1z1 / sigmae2) - logsigmae) ;
      
      b(span(eFirst, eLast)) = (etalpha(span(eFirst, eLast)) % (z1 / sigmae2)) / 
        gutalphasigma2(span(eFirst, eLast)) ;
      
      a(span(eFirst, eLast)) = (-0.5 / sigmae2) / gutalphasigma2(span(eFirst, eLast));
    } else {
      // integration over g1 with e1 = 0
      a(span(eFirst, eLast)) = fe2talpha(span(eFirst, eLast)) / sigma2;
      
      b(span(eFirst, eLast)) = -2 * etalpha(span(eFirst, eLast)) % z1 % a(span(eFirst, eLast));
      
      c(span(eFirst, eLast)) = talpha(span(eFirst, eLast)) + 0.5 * log(-fe2talpha(span(eFirst, eLast))) -
        M_LN_SQRT_PI - logsigma + pow(etalpha(span(eFirst, eLast)) % z1, 2) % a(span(eFirst, eLast));
    }
    
    uint unJ = 0;
    //update parent abcs
    uint lenUnAll = 0;
    while(lenUnAll != eLast - eFirst + 1) {
      a(eReord(span(unIndex(unJ) + 1, unIndex(unJ + 1)))) += 
        a(span(unIndex(unJ) + 1, unIndex(unJ + 1)));
      b(eReord(span(unIndex(unJ) + 1, unIndex(unJ + 1)))) += 
        b(span(unIndex(unJ) + 1, unIndex(unJ + 1)));
      c(eReord(span(unIndex(unJ) + 1, unIndex(unJ + 1)))) += 
        c(span(unIndex(unJ) + 1, unIndex(unJ + 1)));
      
      lenUnAll +=  unIndex(unJ + 1) - unIndex(unJ);
      ++unJ;
    }
    
    // edges pointing to internal nodes
    for(int i = 1; i < nLevels; ++i) {
      uint eFirst = nodesIndex(i) + 1, eLast = nodesIndex(i+1);
      
      // edges pointing to internal nodes, for which all children nodes have been visited
      gutalphasigma2(span(eFirst, eLast)) = e2talpha(span(eFirst, eLast)) + (a(span(eFirst, eLast)) * sigma2) / fe2talpha(span(eFirst, eLast));
      
      c(span(eFirst, eLast)) = -0.5 * log(gutalphasigma2(span(eFirst, eLast))) -
        0.25 * sigma2 * pow(b(span(eFirst, eLast)), 2) /
          (fe2talpha(span(eFirst, eLast)) - alpha + a(span(eFirst, eLast)) * sigma2) +
            talpha(span(eFirst, eLast)) + c(span(eFirst, eLast));
      
      b(span(eFirst, eLast)) = (etalpha(span(eFirst, eLast)) % b(span(eFirst, eLast))) / 
        gutalphasigma2(span(eFirst, eLast));
      
      a(span(eFirst, eLast)) = a(span(eFirst, eLast)) / gutalphasigma2(span(eFirst, eLast));
      
      lenUnAll = 0;
      while(lenUnAll != eLast - eFirst + 1) {
        a(eReord(span(unIndex(unJ) + 1, unIndex(unJ + 1)))) += 
          a(span(unIndex(unJ) + 1, unIndex(unJ + 1)));
        b(eReord(span(unIndex(unJ) + 1, unIndex(unJ + 1)))) += 
          b(span(unIndex(unJ) + 1, unIndex(unJ + 1)));
        c(eReord(span(unIndex(unJ) + 1, unIndex(unJ + 1)))) += 
          c(span(unIndex(unJ) + 1, unIndex(unJ + 1)));
        
        lenUnAll +=  unIndex(unJ + 1) - unIndex(unJ);
        ++unJ;
      }
      
    }
    
    arma::vec res(3);
    res.at(0) = a.at(M-1);
    res.at(1) = b.at(M-1);
    res.at(2) = c.at(M-1);
    return res;
  };
  
  uint get_count_abc_calls() const {
    return count_abc_calls;
  }
};

RCPP_MODULE(IntegratorPOUMM) {
  class_<Integrator>( "Integrator" )
  .constructor()
  .method( "setPruningInfo", &Integrator::setPruningInfo )
  .method( "abc", &Integrator::abc )
  .property("count_abc_calls", &Integrator::get_count_abc_calls )
  ;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
library(POUMM)

tree <- ape::rtree(10)
z <- rnorm(10)
moduleIntegratorPOUMM <- Rcpp::Module( "IntegratorPOUMM", "POUMM" )
Integrator <- moduleIntegratorPOUMM$Integrator

integrator <- Integrator$new(z, tree$edge, tree$edge.length)

prI <- POUMM:::pruneTree(tree, z)
*/
