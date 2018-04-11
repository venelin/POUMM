//   Copyright 2017 Venelin Mitov
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//   http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
//   limitations under the License.
#ifndef ABC_POUMM_H_
#define ABC_POUMM_H_

#include "./SPLITT.h"
#include <iostream>

namespace SPLITT {
template<class NameType>
struct NumericTraitData {
  // use const references to avoid copying of long vectors
  std::vector<NameType> const& names_;
  vec const& z_;
  vec const& se_;
  NumericTraitData(
    std::vector<NameType> const& names,
    vec const& z, vec const& se): names_(names), z_(z), se_(se) {}
};

template<class Tree>
class AbcPOUMM: public TraversalSpecification<Tree> {

public:
  typedef AbcPOUMM<Tree> MyType;
  typedef TraversalSpecification<Tree> BaseType;
  typedef Tree TreeType;
  typedef PostOrderTraversal<MyType> AlgorithmType;
  typedef vec ParameterType;
  typedef NumericTraitData<typename TreeType::NodeType> DataType;
  typedef vec StateType;

  double alpha, theta, sigmae2, sigma2;
  vec z, se;
  vec a, b, c;

  AbcPOUMM(TreeType const& tree, DataType const& input_data):
    BaseType(tree) {

    if(input_data.z_.size() != this->ref_tree_.num_tips() ||
       input_data.se_.size() != this->ref_tree_.num_tips()) {
      std::ostringstream oss;
      oss<<"The vectors z and se must be the same length as the number of tips ("<<
        this->ref_tree_.num_tips()<<"), but were"<<input_data.z_.size()<<
          " and "<<input_data.se_.size()<<" respectively.";
      throw std::invalid_argument(oss.str());
    } else {

      uvec ordNodes = this->ref_tree_.OrderNodes(input_data.names_);
      this->z = At(input_data.z_, ordNodes);
      this->se = At(input_data.se_, ordNodes);
      this->a = vec(this->ref_tree_.num_nodes());
      this->b = vec(this->ref_tree_.num_nodes());
      this->c = vec(this->ref_tree_.num_nodes());
    }
  };

  StateType StateAtRoot() const {
    vec res(3);
    res[0] = a[this->ref_tree_.num_nodes() - 1];
    res[1] = b[this->ref_tree_.num_nodes() - 1];
    res[2] = c[this->ref_tree_.num_nodes() - 1];
    return res;
  };

  void SetParameter(ParameterType const& par) {
    if(par.size() != 4) {
      throw std::invalid_argument(
      "The par vector should be of length 4 with \
      elements corresponding to alpha, theta, sigma and sigmae.");
    }
    if(par[0] < 0 || par[2] < 0 || par[3] < 0) {
      throw std::logic_error("The parameters alpha, sigma and sigmae should be non-negative.");
    }
    this->alpha = par[0];
    this->theta = par[1];
    this->sigma2 = par[2]*par[2];
    this->sigmae2 = par[3]*par[3];
  }

  inline void InitNode(uint i) {
    if(i < this->ref_tree_.num_tips()) {
      double sum_se2_sigmae2 = sigmae2 + se[i]*se[i];
      double z1 = z[i] - theta;
      a[i] = -0.5 / sum_se2_sigmae2;
      b[i] = z1 / sum_se2_sigmae2;
      c[i] = -0.5 * (M_LN_2PI  + z1 * b[i] + log(sum_se2_sigmae2));
    } else {
      a[i] = b[i] = c[i] = 0;
    }
  }

  inline void VisitNode(uint i) {

    double t = this->ref_tree_.LengthOfBranch(i);
    double talpha = t * alpha;
    double etalpha = exp(talpha);
    double e2talpha = etalpha * etalpha;
    double fe2talpha;
    if(alpha != 0) {
      fe2talpha = alpha / (1 - e2talpha);
    } else {
      fe2talpha = -0.5 / t;
    }
    double gutalphasigma2 = e2talpha + (a[i] * sigma2) / fe2talpha;

    c[i] = -0.5 * log(gutalphasigma2) - 0.25 * sigma2 * b[i] * b[i] /
      (fe2talpha - alpha + a[i] * sigma2) + talpha + c[i];
    b[i] = (etalpha * b[i]) / gutalphasigma2;
    a[i] /= gutalphasigma2;
  }

  inline void PruneNode(uint i, uint i_parent) {
    a[i_parent] += a[i];
    b[i_parent] += b[i];
    c[i_parent] += c[i];
  }

};

typedef TraversalTask<
  AbcPOUMM<OrderedTree<uint, double>> > ParallelPruningAbcPOUMM;
}
#endif //ABC_POUMM_H_
