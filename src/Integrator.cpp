#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

class Integrator {
  uvec ZERO;
  uvec ONE;
  uvec TWO;
  
  uint nLevels;
  uint M;
  uint N;
  arma::vec z;
  arma::vec t;
  arma::umat edge;
  arma::uvec endingAt;
  arma::uvec nodesVector;
  arma::uvec nodesIndex;
  arma::uvec unVector;
  arma::uvec unIndex;
  
  arma::mat pif;
  
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
    
    z = vec(z_); t = vec(t_); 
    
    edge = umat(edge_ - 1);
      
    pif = mat(M, 3);
    
    
    // values corrected for 0-based indexing
    endingAt = endingAt_ - 1; 
    nodesVector = nodesVector_ - 1;
    nodesIndex = nodesIndex_ - 1; 
    unIndex = unIndex_ - 1; 
    unVector = unVector_ - 1;
  };
  
  arma::rowvec abc(double alpha, double theta, double sigma, double sigmae) {
    
    ++count_abc_calls;
    
    // time for main-loop
    //auto start = std::chrono::high_resolution_clock::now();
    
    double sigma2 = sigma*sigma, logsigma = log(sigma), 
      logsigmae = log(sigmae), sigmae2 = sigmae*sigmae;
    
    arma::vec talpha = t * alpha;
    arma::vec etalpha = exp(talpha);
    arma::vec e2talpha = etalpha % etalpha;
    arma::vec fe2talpha(M);
    if(alpha != 0) {
      fe2talpha = alpha / (1 - e2talpha);
    } else {
      fe2talpha = -0.5 / t;
    }
    
    
    // needed to re-initialize pif;
    pif.fill(0);
    
    // edges pointing to tips
    uvec es = endingAt.elem(nodesVector(span(nodesIndex(0) + 1, nodesIndex(1))));
    uvec edgeEnds = edge(es, ONE);
    
    vec gutalphasigma2(M);
    
    gutalphasigma2.elem(es) = e2talpha.elem(es) + ((-0.5 / sigmae2) * sigma2) / fe2talpha.elem(es);
    
    // (possibly reordered) shifted tip values
    vec z1 = z.elem(edgeEnds) - theta ;
    
    if(sigmae != 0) {
      // integration over g1 including e1 = z1 - g1
      pif(edgeEnds, TWO) = -0.5 * log(gutalphasigma2.elem(es)) -
        0.25 * sigma2 * pow(z1 / sigmae2, 2) /
          (fe2talpha.elem(es) - alpha + (-0.5 / sigmae2) * sigma2) +
            talpha.elem(es) + (-0.5 * (M_LN_2PI  + z1 % (z1 / sigmae2)) - logsigmae) ;
      
      pif(edgeEnds, ONE) = (etalpha.elem(es) % (z1 / sigmae2)) / gutalphasigma2.elem(es) ;
      
      pif(edgeEnds, ZERO) = (-0.5 / sigmae2) / gutalphasigma2.elem(es);
    } else {
      // integration over g1 with e1 = 0
      pif(edgeEnds, ZERO) = fe2talpha.elem(es) / sigma2;
      
      pif(edgeEnds, ONE) = -2 * etalpha.elem(es) % z1 % pif(edgeEnds, ZERO);
      
      pif(edgeEnds, TWO) = talpha.elem(es) + 0.5 * log(-fe2talpha.elem(es)) -
        M_LN_SQRT_PI - logsigma + pow(etalpha.elem(es) % z1, 2) % pif(edgeEnds, ZERO);
    }
    
    uint unJ = 0;
    
    //update parent pifs
    uint lenUnAll = 0;
    while(lenUnAll != es.n_elem) {
      uvec un = unVector(span(unIndex(unJ) + 1, unIndex(unJ + 1)));
      ++unJ;
      lenUnAll += un.n_elem;
      
      pif.rows(edge(es(un), ZERO)) += pif.rows(edge(es(un), ONE)) ;
    }
    
    // edges pointing to internal nodes
    for(int i = 1; i < nLevels; ++i) {
      es = endingAt.elem(nodesVector(span(nodesIndex(i) + 1, nodesIndex(i + 1))));
      edgeEnds = edge(es, ONE);
      
      // edges pointing to internal nodes, for which all children nodes have been visited
      gutalphasigma2.elem(es) = e2talpha.elem(es) + (pif(edgeEnds, ZERO) * sigma2) / fe2talpha.elem(es);
      
      pif(edgeEnds, TWO) = -0.5 * log(gutalphasigma2.elem(es)) -
        // use the formula log(a+c) = log(a) + log(1+c/a), where
        // a=e2talpha, c = (pif[edgeEnds, 0] * sigma2) / fe2talpha.elem(es), b = e.
        0.25 * sigma2 * pow(pif(edgeEnds, ONE), 2) /
          (fe2talpha.elem(es) - alpha + pif(edgeEnds, ZERO) * sigma2) +
            talpha.elem(es) + pif(edgeEnds, TWO);
      
      pif(edgeEnds, ONE) = (etalpha.elem(es) % pif(edgeEnds, ONE)) / gutalphasigma2.elem(es);
      
      pif(edgeEnds, ZERO) = pif(edgeEnds, ZERO) / gutalphasigma2.elem(es);
      
      //update parent pifs
      lenUnAll = 0;
      while(lenUnAll != es.n_elem) {
        uvec un = unVector(span(unIndex(unJ) + 1, unIndex(unJ + 1)));
        ++unJ;
        lenUnAll += un.n_elem;
        
        pif.rows(edge(es(un), ZERO)) += pif.rows(edge(es(un), ONE)) ;
      }
    }
    
    
    //auto elapsed = std::chrono::high_resolution_clock::now() - start;

    //long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();

    //cout<<microseconds<<endl;


    
    return pif.row(N);
  };
  
  arma::mat get_pif() const {
    return pif;
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
  .property( "pif", &Integrator::get_pif )
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
