#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
#include <math.h>
#include <R_ext/Rdynload.h>

#ifdef _OPENMP

// Need to decide wheter to use '#pragma omp for' or '#pragma omp for simd'
#if _OPENMP >= 201307  // OMP 4.0 or higher

#define _PRAGMA_OMP_FOR_SIMD _Pragma("omp for simd")
#define _PRAGMA_OMP_FOR _Pragma("omp for")
#define _PRAGMA_OMP_SIMD _Pragma("omp simd")

#else // #if _OPENMP >= 201307

#define _PRAGMA_OMP_FOR_SIMD _Pragma("omp for")
#define _PRAGMA_OMP_FOR _Pragma("omp for")
#define _PRAGMA_OMP_SIMD _Pragma("omp simd")

#endif // _OPENMP >= 201307

#include <omp.h>

#else // #ifdef _OPENMP

// the preprocessor directives should simply be ignored at compile-time
#define _PRAGMA_OMP_FOR_SIMD _Pragma("omp for simd")
#define _PRAGMA_OMP_FOR _Pragma("omp for")
#define _PRAGMA_OMP_SIMD _Pragma("omp simd")

#endif // #ifdef _OPENMP




// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;

// Needed to compile on WINDOWS:
typedef unsigned int uint;

// BEGIN: Needed for r-devel (R 3.4)
void R_init_POUMM(DllInfo *info) {
  /* Register routines, allocate resources. */
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}

void R_unload_POUMM(DllInfo *info) {
    /* Release resources. */
}
// END Needed for r-devel (R 3.4)


/* Tuning the schedule at compile time, e.g. : 
 * schedule(static, chunk_size > 1000 ? chunk_size:1000)
 */
//#ifndef _OMP_SCHEDULE
//#define _OMP_SCHEDULE schedule(static)
//#endif // _OMP_SCHEDULE


class Integrator {
  uvec ZERO;
  uvec ONE;
  uvec TWO;
  
  uint nLevels;
  uint M;
  uint N;
  
  arma::vec zReord;
  arma::vec seReord;
  
  bool has_se;
  
  arma::vec sum_se2_sigmae2;
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
      arma::vec z_, arma::vec se_, arma::umat edge_, arma::vec t_,
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
    
    has_se = false;
    
    for(int i = 0; i < se_.n_elem; i++) {
      if(se_[i] > 0) {
        has_se = true;
      }
    }
    
    reorderEdges(edge_ - 1, t_, z_, se_,
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
  
  int check_OPENMP() {
    // thread id and number of threads
    int tid = 0, nthreads = 1;
#ifdef _OPENMP
    Rcpp::Rcout << "Compiled with omp v" << _OPENMP << std::endl; 
    #pragma omp parallel private(tid)
    {
      nthreads = omp_get_num_threads();
      tid = omp_get_thread_num();
      Rcpp::Rcout << "Hello from thread " << tid << " of " << nthreads << std::endl;
    }
#else
    Rcpp::Rcout << "not available" <<endl;
#endif
    return nthreads;
  }
  
  void reorderEdges(arma::umat const& edge, arma::vec const& t, 
                    arma::vec const& z, arma::vec& se,
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
    seReord = se;
    
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
      seReord.subvec(unIndex(unJ) + 1, unIndex(unJ + 1)) = se.elem(edgeEnds);
      
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
  
  arma::vec abc_arma(double alpha, double theta, double sigma, double sigmae) {
    
    ++count_abc_calls;
    
    double sigmae2 = sigmae*sigmae;
    sum_se2_sigmae2 = sigmae2 + seReord % seReord;
    
    double sigma2 = sigma*sigma, logsigma = log(sigma);
    
    arma::vec log_se_total = log(sqrt(sum_se2_sigmae2));
    
    
    
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
      ((-0.5 / sum_se2_sigmae2) * sigma2) / fe2talpha(span(eFirst, eLast));
    
    // (possibly reordered) shifted tip values
    vec z1 = zReord(span(eFirst, eLast)) - theta ;
    
    vec z1z1 = z1 % z1;
    
    if(sigmae != 0 || has_se) {
      // integration over g1 including e1 = z1 - g1
      c(span(eFirst, eLast)) = -0.5 * log(gutalphasigma2(span(eFirst, eLast))) -
        0.25 * sigma2 * z1z1 / (sum_se2_sigmae2 % sum_se2_sigmae2) /
          (fe2talpha(span(eFirst, eLast)) - alpha + (-0.5 / sum_se2_sigmae2) * sigma2) +
      talpha(span(eFirst, eLast)) + (-0.5 * (M_LN_2PI  + z1z1 / sum_se2_sigmae2) - log_se_total) ;
      
      b(span(eFirst, eLast)) = (etalpha(span(eFirst, eLast)) % (z1 / sum_se2_sigmae2)) / 
        gutalphasigma2(span(eFirst, eLast)) ;
      
      a(span(eFirst, eLast)) = (-0.5 / sum_se2_sigmae2) / gutalphasigma2(span(eFirst, eLast));
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
  
  arma::vec abc_omp_simd(double alpha, double theta, double sigma, double sigmae) {
    ++count_abc_calls;
    
    double sigmae2 = sigmae*sigmae;
    sum_se2_sigmae2 = sigmae2 + seReord % seReord;
    
    double sigma2 = sigma*sigma, logsigma = log(sigma);
    
    arma::vec log_se_total = log(sqrt(sum_se2_sigmae2));
    
    a.fill(0); b.fill(0); c.fill(0);
    
    vec& a = this->a;
    vec& b = this->b;
    vec& c = this->c;
    vec& zReord = this->zReord;
    vec& tReord = this->tReord;
    
    vec talpha(M - 1);
    vec etalpha(M - 1);
    vec e2talpha(M - 1);
    vec fe2talpha(M - 1);
    vec gutalphasigma2(M - 1);
    vec z1(N);
    vec z1z1(N);
    
    if(alpha != 0) {
      _PRAGMA_OMP_SIMD
      for(int i = 0; i < M - 1; i++) {
        talpha[i] = tReord[i] * alpha;
        etalpha[i] = exp(talpha[i]);
        e2talpha[i] = etalpha[i] * etalpha[i];
        fe2talpha[i] = alpha / (1 - e2talpha[i]);
      }
    } else {
      _PRAGMA_OMP_SIMD
      for(int i = 0; i < M - 1; i++) {
        talpha[i] = tReord[i] * alpha;
        etalpha[i] = exp(talpha[i]);
        e2talpha[i] = etalpha[i] * etalpha[i];
        fe2talpha[i] = -0.5 / tReord[i];
      }
    }
    
    // edges pointing to tips
    const uint eFirst = nodesIndex[0] + 1;
    const uint eLast = nodesIndex[1];
    
    if(sigmae != 0 || has_se) {
      _PRAGMA_OMP_SIMD
      for(uint i = eFirst; i < eLast + 1; i++) {
        gutalphasigma2[i] = e2talpha[i] + ((-0.5 / sum_se2_sigmae2[i]) * sigma2) / fe2talpha[i];
        z1[i - eFirst] = zReord[i] - theta;
        z1z1[i - eFirst] = z1[i - eFirst] * z1[i - eFirst];
        // integration over g1 including e1 = z1 - g1
        c[i] = -0.5 * log(gutalphasigma2[i]) -
          0.25 * sigma2 * z1z1[i - eFirst] / (sum_se2_sigmae2[i]*sum_se2_sigmae2[i]) /
            (fe2talpha[i] - alpha + (-0.5 / sum_se2_sigmae2[i]) * sigma2) +
              talpha[i] + (-0.5 * (M_LN_2PI  + z1z1[i-eFirst] / sum_se2_sigmae2[i]) - log_se_total[i]);
        b[i] = (etalpha[i] * (z1[i] / sum_se2_sigmae2[i])) / gutalphasigma2[i];
        a[i] = (-0.5 / sum_se2_sigmae2[i]) / gutalphasigma2[i];  
      }
    } else {
      _PRAGMA_OMP_SIMD
      for(uint i = eFirst; i < eLast+1; i++) {
        z1[i - eFirst] = zReord[i] - theta;
        z1z1[i - eFirst] = z1[i - eFirst] * z1[i - eFirst];
        // integration over g1 including e1 = 0
        a[i] = fe2talpha[i] / sigma2;  
        b[i] = -2 * etalpha[i] * z1[i] * a[i];
        c[i] = talpha[i] + 0.5 * log(-fe2talpha[i]) -
          M_LN_SQRT_PI - logsigma + e2talpha[i] * z1z1[i-eFirst] * a[i];
      }
    }
    
    uint unJ = 0;
    //update parent abcs
    uint lenUnAll = 0;
    while(lenUnAll != eLast - eFirst + 1) {
      
      const uint unFirst = unIndex(unJ) + 1;
      const uint unLast = unIndex(unJ + 1);
      
      _PRAGMA_OMP_SIMD
        for(uint i = unFirst; i < unLast + 1; i++) {
          a[eReord[i]] += a[i];
          b[eReord[i]] += b[i];
          c[eReord[i]] += c[i];
        }
        lenUnAll +=  unIndex(unJ + 1) - unIndex(unJ);
      ++unJ;
    }
    
    // edges pointing to internal nodes
    for(int j = 1; j < nLevels; j++) {
      const uint eFirst = nodesIndex(j) + 1;
      const uint eLast = nodesIndex(j+1);
      
      // edges pointing to internal nodes, for which all children nodes have been visited      
      _PRAGMA_OMP_SIMD
        for(uint i = eFirst; i < eLast+1; i++) {
          gutalphasigma2[i] = e2talpha[i] + (a[i] * sigma2) / fe2talpha[i];
          c[i] = -0.5 * log(gutalphasigma2[i]) - 0.25 * sigma2 * b[i] * b[i] /
            (fe2talpha[i] - alpha + a[i] * sigma2) + talpha[i] + c[i];
          b[i] = (etalpha[i] * b[i]) / gutalphasigma2[i];
          a[i] /= gutalphasigma2[i];
        }
        
      lenUnAll = 0;
      while(lenUnAll != eLast - eFirst + 1) {
        const uint unFirst = unIndex(unJ) + 1;
        const uint unLast = unIndex(unJ + 1);
        
        _PRAGMA_OMP_SIMD  
          for(uint i = unFirst; i < unLast + 1; i++) {
            a[eReord[i]] += a[i];
            b[eReord[i]] += b[i];
            c[eReord[i]] += c[i];
          }
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

  arma::vec abc_omp_for_simd(double alpha, double theta, double sigma, double sigmae) {
    ++count_abc_calls;
    
    double sigmae2 = sigmae*sigmae;
    sum_se2_sigmae2 = sigmae2 + seReord % seReord;
    
    double sigma2 = sigma*sigma, logsigma = log(sigma);
    
    arma::vec log_se_total = log(sqrt(sum_se2_sigmae2));
    
    a.fill(0); b.fill(0); c.fill(0);
    
    vec& a = this->a;
    vec& b = this->b;
    vec& c = this->c;
    vec& zReord = this->zReord;
    vec& tReord = this->tReord;
    
    vec talpha(M - 1);
    vec etalpha(M - 1);
    vec e2talpha(M - 1);
    vec fe2talpha(M - 1);
    vec gutalphasigma2(M - 1);
    vec z1(N);
    vec z1z1(N);
    
    uint nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_num_threads();
#endif
    
#pragma omp parallel 
{
  int chunk_size = M/nthreads + 1;
  
  if(alpha != 0) {
    _PRAGMA_OMP_FOR_SIMD
    for(int i = 0; i < M - 1; i++) {
      talpha[i] = tReord[i] * alpha;
      etalpha[i] = exp(talpha[i]);
      e2talpha[i] = etalpha[i] * etalpha[i];
      fe2talpha[i] = alpha / (1 - e2talpha[i]);
    }
  } else {
    _PRAGMA_OMP_FOR_SIMD
    for(int i = 0; i < M - 1; i++) {
      talpha[i] = tReord[i] * alpha;
      etalpha[i] = exp(talpha[i]);
      e2talpha[i] = etalpha[i] * etalpha[i];
      fe2talpha[i] = -0.5 / tReord[i];
    }
  }
  
  // edges pointing to tips
  const uint eFirst = nodesIndex[0] + 1;
  const uint eLast = nodesIndex[1];
  
  chunk_size = (eLast - eFirst + 1)/nthreads + 1;
  
  if(sigmae != 0 || has_se) {
    _PRAGMA_OMP_FOR_SIMD
    for(uint i = eFirst; i < eLast + 1; i++) {
      gutalphasigma2[i] = e2talpha[i] + ((-0.5 / sum_se2_sigmae2[i]) * sigma2) / fe2talpha[i];
      z1[i - eFirst] = zReord[i] - theta;
      z1z1[i - eFirst] = z1[i - eFirst] * z1[i - eFirst];
      // integration over g1 including e1 = z1 - g1
      c[i] = -0.5 * log(gutalphasigma2[i]) -
        0.25 * sigma2 * z1z1[i - eFirst] / (sum_se2_sigmae2[i]*sum_se2_sigmae2[i]) /
          (fe2talpha[i] - alpha + (-0.5 / sum_se2_sigmae2[i]) * sigma2) +
            talpha[i] + (-0.5 * (M_LN_2PI  + z1z1[i-eFirst] / sum_se2_sigmae2[i]) - log_se_total[i]);
      b[i] = (etalpha[i] * (z1[i] / sum_se2_sigmae2[i])) / gutalphasigma2[i];
      a[i] = (-0.5 / sum_se2_sigmae2[i]) / gutalphasigma2[i];  
    }
  } else {
    _PRAGMA_OMP_FOR_SIMD
    for(uint i = eFirst; i < eLast+1; i++) {
      // gutalphasigma2[i] = e2talpha[i] + ((-0.5 / sigmae2) * sigma2) / fe2talpha[i];
      z1[i - eFirst] = zReord[i] - theta;
      z1z1[i - eFirst] = z1[i - eFirst] * z1[i - eFirst];
      // integration over g1 including e1 = 0
      a[i] = fe2talpha[i] / sigma2;  
      b[i] = -2 * etalpha[i] * z1[i] * a[i];
      c[i] = talpha[i] + 0.5 * log(-fe2talpha[i]) -
        M_LN_SQRT_PI - logsigma + e2talpha[i] * z1z1[i-eFirst] * a[i];
    }
  }
  
  uint unJ = 0;
  //update parent abcs
  uint lenUnAll = 0;
  while(lenUnAll != eLast - eFirst + 1) {
    
    const uint unFirst = unIndex(unJ) + 1;
    const uint unLast = unIndex(unJ + 1);
    
    chunk_size = (unLast - unFirst + 1)/nthreads + 1;
    
    _PRAGMA_OMP_FOR_SIMD
      for(uint i = unFirst; i < unLast + 1; i++) {
        a[eReord[i]] += a[i];
        b[eReord[i]] += b[i];
        c[eReord[i]] += c[i];
      }
      lenUnAll +=  unIndex(unJ + 1) - unIndex(unJ);
    ++unJ;
  }
  
  // edges pointing to internal nodes
  for(int j = 1; j < nLevels; j++) {
    const uint eFirst = nodesIndex(j) + 1;
    const uint eLast = nodesIndex(j+1);
    
    chunk_size = (eLast - eFirst + 1)/nthreads + 1;
    
    // edges pointing to internal nodes, for which all children nodes have been visited      
    _PRAGMA_OMP_FOR_SIMD
      for(uint i = eFirst; i < eLast+1; i++) {
        gutalphasigma2[i] = e2talpha[i] + (a[i] * sigma2) / fe2talpha[i];
        c[i] = -0.5 * log(gutalphasigma2[i]) - 0.25 * sigma2 * b[i] * b[i] /
          (fe2talpha[i] - alpha + a[i] * sigma2) + talpha[i] + c[i];
        b[i] = (etalpha[i] * b[i]) / gutalphasigma2[i];
        a[i] /= gutalphasigma2[i];
      }
      
      lenUnAll = 0;
    while(lenUnAll != eLast - eFirst + 1) {
      const uint unFirst = unIndex(unJ) + 1;
      const uint unLast = unIndex(unJ + 1);
      
      chunk_size = (unLast - unFirst + 1)/nthreads + 1;
      
      _PRAGMA_OMP_FOR_SIMD  
        for(uint i = unFirst; i < unLast + 1; i++) {
          a[eReord[i]] += a[i];
          b[eReord[i]] += b[i];
          c[eReord[i]] += c[i];
        }
        lenUnAll +=  unIndex(unJ + 1) - unIndex(unJ);
      ++unJ;
    }
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
  .method( "abc_arma", &Integrator::abc_arma )
  .method( "abc_omp_simd", &Integrator::abc_omp_simd )
  .method( "abc_omp_for_simd", &Integrator::abc_omp_for_simd )
  .method( "check_OPENMP", &Integrator::check_OPENMP )
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
