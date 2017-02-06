#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat mainLoopLikPOUMM_C(
    uint nLevels, uint M, uint N,
    arma::vec z,
    double alpha, double theta, double sigma, double sigmae, double logsigmae,
    double sigma2, double logsigma, double sigmae2, double loge2, double logpi,
    arma::vec talpha, arma::vec etalpha, arma::vec e2talpha, arma::vec fe2talpha,
    arma::umat edge, arma::uvec endingAt,
    arma::uvec nodesVector, arma::uvec nodesIndex,
    arma::uvec unVector, arma::uvec unIndex) {
  
  // correct for 0-based indexing:
  edge -= 1;
  endingAt -= 1;
  nodesVector -= 1;
  nodesIndex -= 1;
  unIndex -= 1;
  unVector -= 1;

  // the following is causing a compilation error "expected expression" in RStudio, but works well in xCode!.
  //uvec ZERO({1});
  //uvec ONE({1});
  //uvec TWO({2});
  
  uvec ZERO(1);
  uvec ONE(1);
  uvec TWO(1);
  
  ZERO.fill(0);
  ONE.fill(1);
  TWO.fill(2);
  
  mat pif(M, 3);
  pif.fill(0);
  
  vec gutalphasigma2 = e2talpha + ((-0.5 / sigmae2) * sigma2) / fe2talpha;
  
  // edges pointing to tips
  uvec es = endingAt.elem(nodesVector(span(nodesIndex(0) + 1, nodesIndex(1))));
  uvec edgeEnds = edge(es, ONE);
  
  if(sigmae != 0) {
    // (possibly reordered) shifted tip values
    vec z1 = z.elem(edgeEnds) - theta ;
    
    pif(edgeEnds, TWO) = -0.5 * log(gutalphasigma2.elem(es)) -
      0.25 * sigma2 * pow(z1 / sigmae2, 2) /
        (fe2talpha.elem(es) - alpha + (-0.5 / sigmae2) * sigma2) +
          talpha.elem(es) + (-0.5 * (loge2 + logpi  + z1 % (z1 / sigmae2)) - logsigmae) ;
    
    pif(edgeEnds, ONE) = (etalpha.elem(es) % (z1 / sigmae2)) / gutalphasigma2.elem(es) ;
    
    pif(edgeEnds, ZERO) = (-0.5 / sigmae2) / gutalphasigma2.elem(es);
  } else {
    // no environmental deviation
    // (possibly reordered) shifted tip values
    vec g1 = z.elem(edgeEnds) - theta;
    
    pif(edgeEnds, ZERO) = fe2talpha.elem(es) / sigma2;    
    
    pif(edgeEnds, ONE) = -2 * etalpha.elem(es) % g1 % pif(edgeEnds, ZERO);
    
    pif(edgeEnds, TWO) = talpha.elem(es) + 0.5 * log(-fe2talpha.elem(es)) -
      0.5*logpi - logsigma + pow(etalpha.elem(es) % g1, 2) % pif(edgeEnds, ZERO);

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
  
  return pif;
}

