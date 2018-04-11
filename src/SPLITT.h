/*
 *  SPLITT.h
 *  SPLITT
 *
 * Copyright 2017 Venelin Mitov
 *
 * This file is part of SPLITT: a generic C++ library for Serial and Parallel
 * Lineage Traversal of Trees.
 *
 * SPLITT is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * SPLITT is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SPLITT.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * @author Venelin Mitov
 */

#ifndef SPLITT_SPLITT_H_
#define SPLITT_SPLITT_H_

#include <algorithm>
#include <vector>
#include <array>
#include <math.h>
#include <sstream>
#include <limits>
#include <numeric>
#include <chrono>
#include <unordered_map>
#include <mutex>
#include <condition_variable>

#ifdef _OPENMP

// Need to decide wheter to use '#pragma omp for' or '#pragma omp for simd'
#if _OPENMP >= 201307  // OMP 4.0 or higher

#define _PRAGMA_OMP_FOR_SIMD _Pragma("omp for simd")
#define _PRAGMA_OMP_FOR _Pragma("omp for")
#define _PRAGMA_OMP_SIMD _Pragma("omp simd")

#else // #if _OPENMP >= 201307

#define _PRAGMA_OMP_FOR_SIMD _Pragma("omp for")
#define _PRAGMA_OMP_FOR _Pragma("omp for")
#define _PRAGMA_OMP_SIMD  /*_Pragma("omp simd")*/

#endif // _OPENMP >= 201307

#include <omp.h>

#else // #ifdef _OPENMP

// the preprocessor directives should simply be ignored at compile-time
#define _PRAGMA_OMP_FOR_SIMD _Pragma("omp for simd")
#define _PRAGMA_OMP_FOR _Pragma("omp for")
#define _PRAGMA_OMP_SIMD _Pragma("omp simd")

#endif // #ifdef _OPENMP


//' @name SPLITT
//' @title namespace SPLITT;
//' 
//' @description All functions and classes defined in the namespace SPLITT.
//' 
//' [[Rcpp::export]]
namespace SPLITT{

//' @name uint
//' @title typedef unsigned int uint;
//' @description a synonyme for the native type.
//' @family basic-types
typedef unsigned int uint;

//' @name uvec
//' @title typedef \href{http://en.cppreference.com/w/cpp/container/vector}{std::vector}<\link{uint}> uvec;
//' @description a vectof of \link{uint}'s.
//' @family basic-types
typedef std::vector<uint> uvec;

//' @name vec
//' @title typedef \href{http://en.cppreference.com/w/cpp/container/vector}{std::vector}<double> vec;
//' @description a vectof of double's.
//' @family basic-types
typedef std::vector<double> vec;

//' @name bvec
//' @title typedef \href{http://en.cppreference.com/w/cpp/container/vector}{std::vector}<bool> bvec;
//' @description a vectof of bool's.
//' @family basic-types
typedef std::vector<bool> bvec;

const uvec EMPTY_UVEC;

// define an NA constant;
const uint NA_UINT = std::numeric_limits<uint>::max();

template <class VectorClass>
inline std::vector<uint> SortIndices(VectorClass const& v) {

  // initialize original index locations
  std::vector<uint> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),
            [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

template<class VectorValues, class VectorPositions>
inline VectorValues At(VectorValues const& v, VectorPositions const& positions) {
  VectorValues sub;
  sub.resize(positions.size());

  size_t sub_i = 0;
  for(auto pit = positions.begin(); pit != positions.end(); pit++,sub_i++){
    sub[sub_i] = v[*pit];
  }
  return sub;
}

template<class VectorValues>
inline VectorValues At(VectorValues const& v, bvec const& mask) {
  if(mask.size() != v.size()) {
    throw std::length_error("ERR:01001:SPLITT:SPLITT.h:At:: bool vector mask should have the same length as v.");
  }

  size_t res_size = 0;
  for(auto b : mask) if(b) ++res_size;

  VectorValues sub(res_size);

  size_t sub_i = 0;
  for(uint i = 0; i < v.size(); i++){
    if(mask[i]) sub[sub_i++] = v[i];
  }
  return sub;
}

// for each element of x return the index of its first occurence in table or
// NA if the element is not found in table or is equal to NA.
// It is assumed that x does not have duplicated elements or NA elements.
template<class VectorValues, class PosType>
inline std::vector<PosType> Match(
    VectorValues const& x, VectorValues const& table, PosType const& NA) {
  auto minmax_x = std::minmax_element(x.begin(), x.end());
  std::vector<PosType> index(*minmax_x.second - *minmax_x.first + 1, NA);
  for(PosType i = 0; i < table.size(); ++i) {
    if(table[i] >= *minmax_x.first && table[i] <= *minmax_x.second &&
       index[table[i] - *minmax_x.first] == NA) {
      index[table[i] - *minmax_x.first] = i;
    }
  }

  std::vector<PosType> positions(x.size());
  for(size_t i = 0; i < x.size(); ++i) {
    positions[i] = index[x[i] - *minmax_x.first];
  }
  return positions;
}

template<class T>
inline std::vector<T> Seq(T const& first, T const& last) {
  std::vector<T> res(last-first+1);
  std::iota(res.begin(), res.end(), first);
  return res;
}

template<class T>
inline bvec IsNA(std::vector<T> const& x, T const& NA) {
  bvec res(x.size(), true);
  for(uint i = 0; i < x.size(); ++i) {
    if(x[i]==NA) res[i] = TRUE;
  }
  return res;
}

inline bvec IsNA(uvec const& x) {
  return IsNA(x, NA_UINT);
}

template<class T>
inline bvec NotIsNA(std::vector<T> const& x, T const& NA) {
  bvec res(x.size(), true);
  for(uint i = 0; i < x.size(); ++i) {
    if(x[i]==NA) res[i] = false;
  }
  return res;
}

inline bvec NotIsNA(uvec const& x) {
  return NotIsNA(x, NA_UINT);
}

//' @name Tree
//' 
//' @title template<class Node, class Length>class Tree
//' 
//' @description A generic C++ class defining the data structure and basic 
//'   operations with a tree. 
//' 
//' @section Template Arguments:
//' \describe{
//' \item{class Node}{see \code{\link{NodeType}}.}
//' \item{class Length}{see \code{\link{LengthType}}.}
//' }
//' @section Public Methods:
//' \describe{
//' }
//' 
//' @seealso \code{\link{OrderedTree}}
template<class Node, class Length>
class Tree {
public:
//' @name Tree::NodeType 
//' @title typedef Node NodeType;
//' @description A public typedef in class \code{\link{Tree}}. A synonym for the 
//'   template argument \code{Node}. Defines a 
//'   hash-able type such as \code{int} or 
//'   \code{\href{http://en.cppreference.com/w/cpp/string/basic_string}{std::string}}. 
//'   Specifically, this should be a type for which 
//'   \code{\href{http://en.cppreference.com/w/cpp/utility/hash}{std::hash}}
//'   specialization does exist. This is the application-specific node-type. The branches in the tree are defined as couples of nodes 
//'   (branch_start_node, branch_end_node).
//' @details
//'   During the construction of \code{Tree} object, the nodes are assigned 
//'   \code{unsigned int} ids from 0 to M-1 (M being the number of nodes in the
//'   tree). 
//'   
//' @aliases NodeType
//' @seealso \code{\link{Tree}} \code{\link{FindIdOfNode}}
  typedef Node NodeType;
  
//' @name Tree::LengthType
//'
//' @title typedef Length LengthType;
//'
//' @description A public typedef in class \code{\link{Tree}}. A synonym for the template
//'   argument Length. Defines a type that can be associated with a branch. Can
//'   be \code{double}, but also a composite of several attributes on a
//'   branch, such as a \code{double} length and an \code{int} color.
//' @aliases LengthType Tree.LengthType
//' @seealso \code{\link{Tree}}
  typedef Length LengthType;

private:
  // private default constructor
  Tree() {}

protected:
  uint num_tips_;
  uint num_nodes_;
  uvec id_parent_;

  typedef std::unordered_map<NodeType, uint> MapType;
  MapType map_node_to_id_;
  std::vector<NodeType> map_id_to_node_;
  std::vector<LengthType> lengths_;
  std::vector<uvec> id_child_nodes_;

  void init_id_child_nodes() {
    id_child_nodes_ = std::vector<uvec>(this->num_nodes() - this->num_tips());

    // fill child vectors
    for(uint i = 0; i < this->num_nodes() - 1; i++) {
      id_child_nodes_[this->FindIdOfParent(i) - this->num_tips()].push_back(i);
    }
  }

public:
//' @name Tree::Tree
//' 
//' @title Constructor for class Tree
//' 
//' @param branch_start_nodes \code{\href{http://en.cppreference.com/w/cpp/container/vector}{std::vector}<\link{NodeType}> const&}; 
//'   starting node for every branch in the tree.
//' @param branch_end_nodes \code{\href{http://en.cppreference.com/w/cpp/container/vector}{std::vector}<\link{NodeType}> const&}; ending 
//'   node for every branch in the tree; must be the same length as 
//'   branch_start_nodes.
//' @param branch_lengths \code{\href{http://en.cppreference.com/w/cpp/container/vector}{std::vector}<\link{LengthType}> const&}; lengths 
//'   associated with the branches. Pass an empty vector for branch_lengths for 
//'   a tree without branch lengths (i.e. only a topology).
//' 
//' @seealso \code{\link{Tree}} \code{\link{OrderedTree}}
  Tree(std::vector<NodeType> const& branch_start_nodes,
       std::vector<NodeType> const& branch_end_nodes,
       std::vector<LengthType> const& branch_lengths) {

    if(branch_start_nodes.size() != branch_end_nodes.size()) {
      std::ostringstream oss;
      oss<<"ERR:01011:SPLITT:SPLITT.h:Tree::"<<
        " branch_start_nodes and branch_end_nodes should be the same size, but were "
         <<branch_start_nodes.size()<<" and "<<
      branch_end_nodes.size()<<" respectively.";
      throw std::length_error(oss.str());
    }

    // There should be exactly num_nodes_ = number-of-branches + 1
    // distinct nodes
    // This is because each branch can be mapped to its ending node. The +1
    // corresponds to the root node, to which no branch points.
    this->num_nodes_ = branch_start_nodes.size() + 1;


    // we distinguish three types of nodes:
    enum NodeRole { ROOT, INTERNAL, TIP };

    // initially, we traverse the list of branches and order the nodes as they
    // appear during the traversal from 0 to num_nodes_-1. This order is done by
    // incrementing node_id_temp.
    uint node_id_temp = 0;

    std::vector<NodeRole> node_types(num_nodes_, ROOT);

    this->map_id_to_node_.resize(num_nodes_);

    this->map_node_to_id_.reserve(num_nodes_);

    uvec branch_starts_temp(branch_start_nodes.size(), NA_UINT);
    uvec branch_ends_temp(branch_start_nodes.size(), NA_UINT);
    uvec ending_at(num_nodes_, NA_UINT);

    std::vector<typename MapType::iterator> it_map_node_to_id_;
    it_map_node_to_id_.reserve(num_nodes_);


    for(uint i = 0; i < branch_start_nodes.size(); ++i) {
      if(branch_start_nodes[i] == branch_end_nodes[i]) {
        std::ostringstream oss;
        oss<<"ERR:01012:SPLITT:SPLITT.h:Tree:: Found a branch with the same start and end node ("<<
          branch_start_nodes[i]<<"). Not allowed. ";
        throw std::logic_error(oss.str());
      }


      auto it1 = map_node_to_id_.insert(
        std::pair<NodeType, uint>(branch_start_nodes[i], node_id_temp));

      if(it1.second) {
        // node encountered for the first time and inserted in the map_node_to_id_
        map_id_to_node_[node_id_temp] = branch_start_nodes[i];
        if(node_types[node_id_temp] == TIP) {
          node_types[node_id_temp] = INTERNAL;
        }
        branch_starts_temp[i] = node_id_temp;
        it_map_node_to_id_.push_back(it1.first);
        node_id_temp++;
      } else {
        // node encountered in a previous branch
        if(node_types[it1.first->second] == TIP) {
          // the previous encounter of the node was as a branch-end
          node_types[it1.first->second] = INTERNAL;
        } else {
          // do nothing
        }
        branch_starts_temp[i] = it1.first->second;
      }


      auto it2 = map_node_to_id_.insert(std::pair<NodeType, uint>(branch_end_nodes[i], node_id_temp));

      if(it2.second) {
        // node encountered for the first time and inserted in the map_node_to_id
        map_id_to_node_[node_id_temp] = branch_end_nodes[i];

        if(node_types[node_id_temp] == ROOT) {
          // not known if the node has descendants, so we set its type to TIP.
          node_types[node_id_temp] = TIP;
        }
        branch_ends_temp[i] = node_id_temp;
        ending_at[node_id_temp] = i;
        it_map_node_to_id_.push_back(it2.first);
        node_id_temp++;

      } else {
        // node has been previously encountered
        if(ending_at[it2.first->second] != NA_UINT) {
          std::ostringstream oss;
          oss<<"ERR:01013:SPLITT:SPLITT.h:Tree:: Found at least two branches ending at the same node ("<<
            it2.first->first<<"). Check for cycles or repeated branches. ";
          throw std::logic_error(oss.str());
        } else {
          if(node_types[it2.first->second] == ROOT) {
            // the previous enounters of the node were as branch-start -> set
            // the node's type to INTERNAL, because we know for sure that it
            // has descendants.
            node_types[it2.first->second] = INTERNAL;
          }
          branch_ends_temp[i] = it2.first->second;
          ending_at[it2.first->second] = i;
        }
      }
    }

    if(map_node_to_id_.size() != num_nodes_) {
      std::ostringstream oss;
      oss<<"ERR:01014:SPLITT:SPLITT.h:Tree:: The number of distinct nodes ("<<map_node_to_id_.size()<<
        ") should equal the number-of-branches+1 ("<<num_nodes_<<").";
      throw std::logic_error(oss.str());
    }

    auto num_roots = count(node_types.begin(), node_types.end(), ROOT);
    if(num_roots != 1) {
      std::ostringstream oss;
      oss<<"ERR:01015:SPLITT:SPLITT.h:Tree:: There should be exactly one ROOT node, but "<<num_roots<<
        " were found. Check for cycles or for multiple trees.";
      throw std::logic_error(oss.str());
    }

    this->num_tips_ = count(node_types.begin(), node_types.end(), TIP);
    if(num_tips_ == 0) {
      std::ostringstream oss;
      oss<<"ERR:01016:SPLITT:SPLITT.h:Tree:: There should be at least one TIP node, but none"<<
        " was found. Check for cycles.";
      throw std::logic_error(oss.str());
    }

    // assign new ids according to the following convention:
    // tips are numbered from 0 to num_tips_ - 1;
    // internal nodes are numbered from num_tips_ to num_nodes_ - 2;
    // root is numbered num_nodes_ - 1;
    std::vector<uint> node_ids(num_nodes_, NA_UINT);
    uint tip_no = 0, internal_no = num_tips_;
    for(uint i = 0; i < num_nodes_; i++) {
      if(node_types[i] == TIP) {
        node_ids[i] = tip_no;
        tip_no ++;
      } else if(node_types[i] == INTERNAL) {
        node_ids[i] = internal_no;
        internal_no++;
      } else {
        // Here node_types[i] == ROOT should be true
        node_ids[i] = num_nodes_ - 1;
      }
      //map_node_to_id_[map_id_to_node_[i]] = node_ids[i];
      it_map_node_to_id_[i]->second = node_ids[i];
    }


    this->map_id_to_node_ = At(map_id_to_node_, SortIndices(node_ids));

    this->id_parent_ = uvec(num_nodes_ - 1);

    if(branch_lengths.size() == num_nodes_ - 1) {
      this->lengths_ = std::vector<LengthType>(num_nodes_ - 1);
    } else if(branch_lengths.size() != 0) {
      std::ostringstream oss;
      oss<<"ERR:01017:SPLITT:SPLITT.h:Tree:: branch_lengths should be either empty or of size num_nodes_-1 ("<<
        num_nodes_-1<<") but is "<<branch_lengths.size()<<"."<<std::endl;
      throw std::invalid_argument(oss.str());
    }


    if(HasBranchLengths()) {
      for(uint i = 0; i < num_nodes_ - 1; i++) {
        uint branch_start_i = node_ids[branch_starts_temp[i]];
        uint branch_end_i = node_ids[branch_ends_temp[i]];
        id_parent_[branch_end_i] = branch_start_i;
        lengths_[branch_end_i] = branch_lengths[i];
      }
    } else {
      for(uint i = 0; i < num_nodes_ - 1; i++) {
        uint branch_start_i = node_ids[branch_starts_temp[i]];
        uint branch_end_i = node_ids[branch_ends_temp[i]];
        id_parent_[branch_end_i] = branch_start_i;
      }
    }

    init_id_child_nodes();
  }

//' @name Tree::num_nodes
//' @title Number of nodes in a tree
//'  
//' @section Signature:
//' \code{\link{uint} Tree::num_nodes() const;}
//'  
//' @return the numbef of nodes in the tree, including tips, internal nodes and the root.
//' @seealso \code{\link{Tree}} \code{\link{OrderedTree}}
  uint num_nodes() const {
    return num_nodes_;
  }

//' @name Tree::num_tips
//' 
//' @title Number of tips (a.k.a. leaves) in the tree
//' @section Signature:
//' \code{\link{uint} \link{Tree}::num_nodes() const;}
//' @return \code{\link{uint}}, the number of tips in the tree.
//' @seealso \code{\link{Tree}}
  uint num_tips() const {
    return num_tips_;
  }
  
//' @name Tree::HasBranchLengths
//' 
//' @title Does a tree has lengths associated with its branches?
//' @section Signature:
//' \code{bool \link{Tree}::HasBranchLengths() const;}
//' @return \code{bool}, \code{true} if the tree has branch lengths.
//' @aliases HasBranchLengths
//' @family branch-length-methods
//' @seealso \code{\link{Tree}} \code{\link{OrderedTree}}
  bool HasBranchLengths() const {
    return lengths_.size() == id_parent_.size();
  }

//' @name Tree::LengthOfBranch
//' 
//' @title Get the length of a branch ending at node with id \code{i}
//' @section Signature:
//' \code{\link{LengthType} const& \link{Tree}::LengthOfBranch(uint i) const;}
//' @param i \code{\link{uint}}; the id of the end-node for the branch
//' @return \code{\link{LengthType}}; the length associated with the branch ending at node \code{i}.
//' @aliases LengthOfBranch
//' @family branch-length-methods
//' @seealso \code{\link{Tree}} \code{\link{OrderedTree}}
  LengthType const& LengthOfBranch(uint i) const {
    if(i >= lengths_.size()) {
      std::ostringstream oss;
      oss<<"ERR:01021:SPLITT:SPLITT.h:LengthOfBranch:: i is beyond the size of the lengths_ vector."<<
        "Check i and that the tree has branches."<<std::endl;
    }
    return lengths_[i];
  }
  
//' @name Tree::BranchLengths
//' 
//' @title Get a reference to the internal vector of branch lengths.
//' @section Signature:
//' \code{\href{http://en.cppreference.com/w/cpp/container/vector}{std::vector}<\link{LengthType}> const& \link{Tree}::BranchLengths() const;}
//' @return \code{\href{http://en.cppreference.com/w/cpp/container/vector}{std::vector}<\link{LengthType}> const&}; 
//'   a const reference to the internally stored vector of branch lengths, in the order of the end-node ids.
//' @aliases BranchLengths
//' @family branch-length-methods
//' @seealso \code{\link{Tree}} \code{\link{OrderedTree}}
  std::vector<LengthType> const& BranchLengths() const {
    return lengths_;
  }

//' @name Tree::SetLengthOfBranch
//' 
//' @title Set the length of a branch ending at node with id \code{i} to a given \code{value}
//' @section Signature:
//' \code{void \link{Tree}::SetLengthOfBranch(uint i, \link{LengthType} const& value);}
//' @param i \code{\link{uint}}; the id of the end-node of the branch;
//' @param value \code{\link{LengthType} const&}; the new value of to set.
//' @return \code{void}
//' @family branch-length-methods
//' @aliases SetLengthOfBranch
//' @seealso \code{\link{Tree}} \code{\link{OrderedTree}}
  void SetLengthOfBranch(uint i, LengthType const& value) {
    if(!HasBranchLengths()) {
      std::ostringstream oss;
      oss<<"ERR:01031:SPLITT:SPLITT.h:SetLengthOfBranch:: Trying to set a branch length on a tree without branch lengths. "<<
        "Use a SetBranchLengths method to add branch lengths first."<<std::endl;
      throw std::logic_error(oss.str());
    } else if(i >= lengths_.size()) {
      std::ostringstream oss;
      oss<<"i should be smaller than "<<lengths_.size()<<" but was "<<i<<std::endl;
      throw std::out_of_range(oss.str());
    } else {
      lengths_[i] = value;
    }
  }

//' @name Tree::SetBranchLengths
//' 
//' @title Set a new internally stored vector of branch lenthts
//' 
//' @description 
//' @section Signature:
//' \code{void \link{Tree}::SetBranchLengths(\href{http://en.cppreference.com/w/cpp/container/vector}{std::vector}<\link{LengthType}> const& lengths);}
//' @param lengths \code{\href{http://en.cppreference.com/w/cpp/container/vector}{std::vector}<\link{LengthType}> const&}; 
//'   a const reference to a new vector of branch lengths, in the order of the end-node ids. 
//'   The vector should be of length M-1, where M is the 
//'   number of nodes in the tree.
//' @return \code{void}
//' @aliases SetBranchLengths
//' @family branch-length-methods
//' @seealso \code{\link{Tree}} \code{\link{OrderedTree}}
  void SetBranchLengths(std::vector<LengthType> const& lengths) {
    if(lengths.size() != 0 && lengths.size() != num_nodes_ - 1) {
      std::ostringstream oss;
      oss<<"ERR:01041:SPLITT:SPLITT.h:SetBranchLengths:: lengths should be either empty or of size num_nodes_-1 ("<<
        num_nodes_-1<<") but is "<<lengths.size()<<"."<<std::endl;
    } else {
      lengths_ = lengths;
    }
  }

//' @name Tree::SetBranchLengths2
//' 
//' @title Set new branch lenthts to the branches in the order given by their 
//'   application-specific end-nodes.
//' @description If the tree has no branch lengths, the supplied arguments should
//'   be of length M-1, where M is the total number of nodes in the tree (-1, 
//'   because there is no branch leading to the root).
//' @section Signature:
//' \code{void \link{Tree}::SetBranchLengths(\href{http://en.cppreference.com/w/cpp/container/vector}{std::vector}<\link{NodeType}> const& nodes_branch_ends, \href{http://en.cppreference.com/w/cpp/container/vector}{std::vector}<\link{LengthType}> const& lengths);}
//' 
//' @param nodes_branch_ends \code{\href{http://en.cppreference.com/w/cpp/container/vector}{std::vector}<\link{NodeType}> const&}; 
//'   a const reference to a new vector of branch lengths, in the order of the nodes in \code{nodes_branch_ends}.
//' @param lengths \code{\href{http://en.cppreference.com/w/cpp/container/vector}{std::vector}<\link{LengthType}> const&}; 
//'   a const reference to a new vector of branch lengths, in the order of the nodes in \code{nodes_branch_ends}.
//'   
//' @return \code{void}
//' @aliases SetBranchLengths
//' @family branch-length-methods
//' @seealso \code{\link{Tree}} \code{\link{OrderedTree}}
  void SetBranchLengths(std::vector<NodeType> const& nodes_branch_ends,
                        std::vector<LengthType> const& lengths) {
    if(nodes_branch_ends.size() != lengths.size()) {
      throw std::invalid_argument("ERR:01051:SPLITT:SPLITT.h:SetBranchLengths:: The vectors nodes_branch_ends and lengths should be the same size.");
    }
    if( !HasBranchLengths() ) {
      if(nodes_branch_ends.size() != num_nodes_ - 1) {
        std::ostringstream oss;
        oss<<"ERR:01052:SPLITT:SPLITT.h:SetBranchLengths:: Trying to set branch lengths on a tree without such."<<
          " In this case, the vectors nodes_branch_ends and lengths should have an"<<
            "element for each branch but their size is "<<
              nodes_branch_ends.size()<<" (should be "<<num_nodes_ - 1<<")."<<std::endl;
        throw std::invalid_argument(oss.str());
      }
      lengths_ = vec(num_nodes_ - 1);
    }
    std::vector<bool> visited(num_nodes_ - 1, false);
    for(uint i = 0; i < nodes_branch_ends.size(); ++i) {
      uint id = FindIdOfNode(nodes_branch_ends[i]);
      if(i == NA_UINT || i == num_nodes_ - 1) {
        std::ostringstream oss;
        oss<<"ERR:01053:SPLITT:SPLITT.h:SetBranchLengths:: No branch ends at node identified as "<<id<<
          ". Check that nodes_branch_ends correspond to tips or internal nodes (excluding the root)"<<std::endl;
        throw std::logic_error(oss.str());
      } else if(visited[id]) {
        std::ostringstream oss;
        oss<<"ERR:01054:SPLITT:SPLITT.h:SetBranchLengths:: Trying to set the length of the same branch twice. Check nodes_branch_ends for duplicates."<<std::endl;
        throw std::logic_error(oss.str());
      }
      visited[id] = true;
      lengths_[id] = lengths[i];
    }

  }

//' @name Tree::FindNodeWithId
//' @title Get the node with the specified id.
//' @section Signature:
//' \code{\link{NodeType} const& \link{Tree}::FindNodeWithId(uint id) const;}
//' 
//' @description A public method of class \code{\link{Tree}}.
//' 
//' @param id \code{\link{uint}} the id of the node (should be between 0 and M-1, where
//'   M is the number of nodes in hte tree).
//' @family node-methods  
//' 
//' @aliases FindNodeWithId
//' 
//' @seealso \code{\link{Tree}} \code{\link{OrderedTree}}
  NodeType const& FindNodeWithId(uint id) const {
    return map_id_to_node_[id];
  }

//' @name Tree::FindIdOfNode
//' @title Get the internally stored id of a node
//' 
//' @section Signature:
//' \code{uint \link{Tree}::FindIdOfNode(\link{NodeType} const& node) const;}
//' 
//' @description A public method of class \code{\link{Tree}}.
//' 
//' @param node \code{\link{NodeType} const&}; the node;
//' 
//' @return an \code{\link{uint}} from 0 to M-1, where M is the number of nodes 
//'   in the tree.
//' @family node-methods 
//' 
//' @aliases FindIdOfNode
//' 
//' @seealso \code{\link{Tree}} \code{\link{OrderedTree}}
  uint FindIdOfNode(NodeType const& node) const {
    auto it = map_node_to_id_.find(node);
    if(it == map_node_to_id_.end()) {
      return NA_UINT;
    } else {
      return it->second;
    }
  }

//' @name Tree::FindIdOfParent
//' @title Get the parent id of a node with id id_child.
//' @section Signature:
//' \code{uint \link{Tree}::FindIdOfParent(uint id_child) const;}
//' 
//' @param id_child \code{\link{uint}}, the id of the child node. 
//' 
//' @return an \code{\link{uint}} from 0 to M-1, where M is the number of nodes 
//'   in the tree.
//'   
//' @family node-methods
//' 
//' @aliases FindIdOfParent
//'  
//' @seealso \code{\link{Tree}} \code{\link{OrderedTree}}
  uint FindIdOfParent(uint id_child) const {
    return this->id_parent_[id_child];
  }

//' @name Tree::FindChildren
//' @title Get a vector with the ids of the children of a node.
//' @section Signature:
//' \code{uvec const& \link{Tree}::FindChildren(uint i) const;}
//' 
//' @param i \code{\link{uint}}; the id of a node. 
//' 
//' @return a \code{uvec const&}; a const reference to an internally stored vector of the
//'   ids of the children of \code{i}. If \code{i} is a tip, an empty uvec is returned. The
//'   returned vector reference is valid as long as the tree object is not destroyed.
//'   
//' @family node-methods 
//' 
//' @aliases FindChildren
//' 
//' @seealso \code{\link{Tree}} \code{\link{OrderedTree}}
  uvec const& FindChildren(uint i) const {
    if(i < this->num_tips()) {
      return EMPTY_UVEC;
    } else if(i - this->num_tips() < id_child_nodes_.size()) {
      return id_child_nodes_[i - this->num_tips()];
    } else {
      throw std::invalid_argument("ERR:01061:SPLITT:SPLITT.h:FindChildren:: i must be smaller than the number of nodes.");
    }
  }

//' @name Tree::OrderNodes
//' @title Reorder a vector of nodes
//' @section Signature:
//' \code{\link{uvec} \link{Tree}::OrderNodes(\href{http://en.cppreference.com/w/cpp/container/vector}{std::vector}<\link{NodeType}> const& nodes) const;}
//' @param nodes \code{\href{http://en.cppreference.com/w/cpp/container/vector}{std::vector}<\link{NodeType}> const&}; 
//'   a vector of application-specific nodes.
//' @return \code{\link{uvec}}; a vector of positions in \code{nodes} in the order of 
//'   their internally stored ids.
//' @family node-methods
//' @aliases OrderNodes
//' @seealso \code{\link{Tree}} \code{\link{OrderedTree}}
  uvec OrderNodes(std::vector<NodeType> const& nodes) const {
    return OrderNodesPosType(nodes, NA_UINT);
  }
  

//' @name Tree::OrderNodesPosType
//' @title Reorder a vector of nodes (generic w.r.t. the position type).
//' @section Template Arguments:
//' \describe{
//' \item{PosType}{an integer type for the positions in the returned vector.}}
//' @section Signature:
//' \code{template<class PosType> \href{http://en.cppreference.com/w/cpp/container/vector}{std::vector}<PosType> \link{Tree}::OrderNodesPosType(\href{http://en.cppreference.com/w/cpp/container/vector}{std::vector}<\link{NodeType}> const& nodes, PosType const& NA) const;
//' @param nodes \code{\href{http://en.cppreference.com/w/cpp/container/vector}{std::vector}<PosType> \link{Tree}::OrderNodesPosType(\href{http://en.cppreference.com/w/cpp/container/vector}{std::vector}<\link{NodeType}> const&}; 
//'   a vector of application-specific nodes.
//' @param NA \code{PosType const&}; NA value used mainly for the purpose of template-specification.
//' @return \code{\href{http://en.cppreference.com/w/cpp/container/vector}{std::vector}<PosType>}; a vector of positions in nodes in the order of their internally stored ids.
//' the element-type of the returned vector can be specified as a template argument or 
//' dedcued from the type of NA.
//' @aliases OrderNodesPosType
//' @family node-methods
//' @seealso \code{\link{Tree}} \code{\link{OrderedTree}}
  template<class PosType>
  std::vector<PosType> OrderNodesPosType(std::vector<NodeType> const& nodes, PosType const& NA) const {
    uvec ids(nodes.size());
    for(uint i = 0; i < nodes.size(); ++i) {
      auto it = this->map_node_to_id_.find(nodes[i]);
      if(it == this->map_node_to_id_.end()) {
        std::ostringstream oss;
        oss<<"ERR:01071:SPLITT:SPLITT.h:OrderNodesPosType:: At least one of the nodes is not present in the tree ("<<
          nodes[i]<<").";
        throw std::invalid_argument(oss.str());
      } else {
        ids[i] = it->second;
      }
    }
    std::vector<PosType> m = Match(Seq(uint(0), this->num_nodes_ - 1), ids, NA);
    return At(m, NotIsNA(m, NA));
  }
};

//' @name OrderedTree
//' 
//' @title template<class Node, class Length>class OrderedTree: public Tree<Node, Length>
//' 
//' @description A generic C++ class defining the data structure and basic 
//'   operations with a tree. 
//' 
//' @section Template Arguments:
//' \describe{
//' \item{class Node}{see \code{\link{NodeType}}.}
//' \item{class Length}{see \code{\link{LengthType}}.}
//' }
//' @section Public Methods:
//' \describe{
//' }
//' 
//' @seealso \code{\link{OrderedTree}}
template<class Node, class Length>
class OrderedTree: public Tree<Node, Length> {
public:
  typedef Node NodeType;
  typedef Length LengthType;

private:
  // default constructor;
  OrderedTree() {}

protected:
  uvec ranges_id_visit_;
  uvec ranges_id_prune_;

public:

  OrderedTree(
    std::vector<NodeType> const& branch_start_nodes,
    std::vector<NodeType> const& branch_end_nodes,
    std::vector<LengthType> const& branch_lengths):
  Tree<NodeType, LengthType>(branch_start_nodes, branch_end_nodes, branch_lengths),
  ranges_id_visit_(1, 0),
  ranges_id_prune_(1, 0) {

    // insert a fictive branch leading to the root of the tree.
    uvec branch_ends = Seq(uint(0), this->num_nodes_ - 1);

    uvec num_children_remaining(this->num_nodes_, 0);
    for(uint i : this->id_parent_) num_children_remaining[i]++;

    // start by pruning the tips of the tree
    uvec tips_this_level = Seq(uint(0), this->num_tips_ - 1);

    uvec default_pos_vector;
    default_pos_vector.reserve(2);
    std::vector<uvec> pos_of_parent(this->num_nodes_ - this->num_tips_,
                                    default_pos_vector);

    uvec order_branches;
    order_branches.reserve(this->num_nodes_);

    while(tips_this_level[0] != this->num_nodes_ - 1) {
      // while the root has not become a tip itself
      ranges_id_visit_.push_back(
        ranges_id_visit_[ranges_id_visit_.size() - 1] + tips_this_level.size());

      // unique parents at this level
      uvec parents_this_level;
      parents_this_level.reserve(tips_this_level.size());

      uvec tips_next_level;
      tips_next_level.reserve(tips_this_level.size() / 2);

      for(uint i = 0; i < tips_this_level.size(); i++) {
        uint i_parent = this->id_parent_[tips_this_level[i]];
        if(pos_of_parent[i_parent - this->num_tips_].empty()) {
          parents_this_level.push_back(i_parent);
        }
        pos_of_parent[i_parent - this->num_tips_].push_back(i);
      }

      uint num_parents_remaining = parents_this_level.size();
      while( num_parents_remaining ) {

        uint num_parent_updates = 0;
        for(auto i_parent: parents_this_level) {
          if(!pos_of_parent[i_parent - this->num_tips_].empty()) {
            uint i = pos_of_parent[i_parent - this->num_tips_].back();

            num_parent_updates ++;
            order_branches.push_back(tips_this_level[i]);

            num_children_remaining[i_parent]--;
            if(num_children_remaining[i_parent] == 0) {
              tips_next_level.push_back(i_parent);
            }
            pos_of_parent[i_parent - this->num_tips_].pop_back();
            if(pos_of_parent[i_parent - this->num_tips_].empty()) {
              num_parents_remaining--;
            }
          }
        }

        ranges_id_prune_.push_back(
          ranges_id_prune_[ranges_id_prune_.size() - 1] + num_parent_updates);
      }

      tips_this_level = tips_next_level;
    }

    if(this->HasBranchLengths()) {
      this->lengths_ = At(this->lengths_, order_branches);
    }

    uvec id_old = order_branches;
    id_old.push_back(this->num_nodes_ - 1);

    this->id_parent_ = Match(At(this->id_parent_, order_branches),
                             id_old, NA_UINT);

    // update maps
    std::vector<NodeType> map_id_to_node(this->num_nodes_);
    for (uint i = 0; i < this->num_nodes_; i++) {
      map_id_to_node[i] = this->map_id_to_node_[id_old[i]];
      this->map_node_to_id_[map_id_to_node[i]] = i;
    }

    std::swap(this->map_id_to_node_, map_id_to_node);

    this->init_id_child_nodes();
  }

  uint num_levels() const {
    return ranges_id_visit_.size() - 1;
  }

  uint num_parallel_ranges_prune() const {
    return ranges_id_prune_.size() - 1;
  }

  uvec const& ranges_id_visit() const {
    return ranges_id_visit_;
  }

  
  std::array<uint, 2> RangeIdVisitNode(uint i_level) const {
    // double braces required by C++11 standard,
    // http://en.cppreference.com/w/cpp/container/array
    return std::array<uint, 2> {{ranges_id_visit_[i_level],
                                 ranges_id_visit_[i_level+1] - 1}};
  }

  uvec const& ranges_id_prune() const {
    return ranges_id_prune_;
  }

  std::array<uint, 2> RangeIdPruneNode(uint i_step) const {
    return std::array<uint, 2> {{ranges_id_prune_[i_step],
                                 ranges_id_prune_[i_step+1] - 1}};
  }
};

enum PostOrderMode {
  AUTO = 0,
  SINGLE_THREAD_LOOP_POSTORDER = 10,
  SINGLE_THREAD_LOOP_PRUNES = 11,
  SINGLE_THREAD_LOOP_VISITS = 12,
  MULTI_THREAD_LOOP_PRUNES = 21,
  MULTI_THREAD_LOOP_VISITS = 22,
  MULTI_THREAD_LOOP_VISITS_THEN_LOOP_PRUNES = 23,
  MULTI_THREAD_VISIT_QUEUE = 24,
  HYBRID_LOOP_PRUNES = 31,
  HYBRID_LOOP_VISITS = 32,
  HYBRID_LOOP_VISITS_THEN_LOOP_PRUNES = 33
};

inline std::ostream& operator<< (std::ostream& os, PostOrderMode mode) {
  switch(mode) {
  case PostOrderMode::AUTO: os<<"AUTO"; break;
  case PostOrderMode::SINGLE_THREAD_LOOP_POSTORDER: os<<"SINGLE_THREAD_LOOP_POSTORDER"; break;
  case PostOrderMode::SINGLE_THREAD_LOOP_PRUNES: os<<"SINGLE_THREAD_LOOP_PRUNES"; break;
  case PostOrderMode::SINGLE_THREAD_LOOP_VISITS: os<<"SINGLE_THREAD_LOOP_VISITS"; break;
  case PostOrderMode::MULTI_THREAD_LOOP_PRUNES: os<<"MULTI_THREAD_LOOP_PRUNES"; break;
  case PostOrderMode::MULTI_THREAD_LOOP_VISITS: os<<"MULTI_THREAD_LOOP_VISITS"; break;
  case PostOrderMode::MULTI_THREAD_LOOP_VISITS_THEN_LOOP_PRUNES: os<<"MULTI_THREAD_LOOP_VISITS_THEN_LOOP_PRUNES"; break;
  case PostOrderMode::MULTI_THREAD_VISIT_QUEUE: os<<"MULTI_THREAD_VISIT_QUEUE"; break;
  case PostOrderMode::HYBRID_LOOP_PRUNES: os<<"HYBRID_LOOP_PRUNES"; break;
  case PostOrderMode::HYBRID_LOOP_VISITS: os<<"HYBRID_LOOP_VISITS"; break;
  case PostOrderMode::HYBRID_LOOP_VISITS_THEN_LOOP_PRUNES: os<<"HYBRID_LOOP_VISITS_THEN_LOOP_PRUNES"; break;
  };
  return os<< static_cast<int>(mode);
}

template<class TreeType> class VisitQueue {
  std::mutex mutex_;
  std::condition_variable has_a_new_node_;
  
  TreeType const& ref_tree_;
  uvec queue_;
  uvec::iterator it_queue_begin;
  uvec::iterator it_queue_end;
  uvec num_non_visited_children_;
public:
  
  // non-thread safe (called in single-thread mode)
  void Init(uvec const& num_children) {
    std::copy(num_children.begin(), num_children.end(),
              num_non_visited_children_.begin());
    it_queue_begin = queue_.begin();
    it_queue_end = queue_.begin() + ref_tree_.num_tips();
    std::iota(it_queue_begin, it_queue_end, 0);
  }
  
  bool IsTemporarilyEmpty() const {
    return it_queue_begin == it_queue_end & it_queue_end < queue_.end();
  }
  
  // thread-safe
  uint NextInQueue() {
    std::unique_lock<std::mutex> lock(mutex_);
    
    while( IsTemporarilyEmpty() ) {
      has_a_new_node_.wait(lock);
    }
    
    if(it_queue_begin < it_queue_end) {
      uint res = *it_queue_begin;
      ++it_queue_begin;
      return res;
    } else if(it_queue_begin == queue_.end()) {
      // algorithm thread should stop here. all waiting threads should be notified,
      // since no other elements will be inserted in the queue.
      has_a_new_node_.notify_all();
      return ref_tree_.num_nodes();
    } else {
      // algorithm thread continues to check for new node to visit
      // should never execute this
      //std::cout<<"Error returning NA_UINT from VisitQueue."<<std::endl;
      return NA_UINT;
    }
  }
  
  // thread-safe
  // if the parent of i becomes visit-able, it gets inserted in the
  // queue.
  void RemoveVisitedNode(uint i) {
    std::unique_lock<std::mutex> lock(mutex_);
    
    uint i_parent = ref_tree_.FindIdOfParent(i);
    num_non_visited_children_[i_parent - ref_tree_.num_tips()]--;
    if(num_non_visited_children_[i_parent - ref_tree_.num_tips()] == 0) {
      *it_queue_end = i_parent;
      *it_queue_end++;
      has_a_new_node_.notify_one();
    }
  }
  
  // non-thread-safe. should call Init() before using.
  VisitQueue(TreeType const& tree):
  ref_tree_(tree),
  queue_(tree.num_nodes()),
  it_queue_begin(queue_.begin()),
  it_queue_end(queue_.begin()),
  num_non_visited_children_(tree.num_nodes() - tree.num_tips()) {}
  
  // Copy initialization (non-thread-safe)
  VisitQueue(const VisitQueue& other): ref_tree_(other.ref_tree_) {
    auto other_begin = other.queue_.begin();
    queue_ = other.queue_;
    it_queue_begin = queue_.begin() + (other.it_queue_begin  - other_begin);
    it_queue_end = queue_.begin() + (other.it_queue_end  - other_begin);
    num_non_visited_children_ = other.num_non_visited_children_;
  }
};

template<class TraversalSpecification>
class TraversalAlgorithm {
protected:
  typedef typename TraversalSpecification::TreeType TreeType;

  TreeType const& ref_tree_;
  TraversalSpecification& ref_spec_;

  uint num_threads_;

  uvec num_children_;
  VisitQueue<TreeType> visit_queue_;

public:
  TraversalAlgorithm(TreeType const& tree, TraversalSpecification& spec):
  ref_tree_(tree),
  ref_spec_(spec),
  num_children_(tree.num_nodes() - tree.num_tips()),
  visit_queue_(tree) {

#ifdef _OPENMP
#pragma omp parallel
{
  uint tid = omp_get_thread_num();
  // only master thread does this
  if(tid == 0) {
    this->num_threads_ = omp_get_num_threads();
  }
}
#else
this->num_threads_ = 1;
#endif // #ifdef _OPENMP

    for(uint i = tree.num_tips(); i < tree.num_nodes(); i++) {
      num_children_[i - tree.num_tips()] = tree.FindChildren(i).size();
    }
  }

  uint num_threads() const {
    return num_threads_;
  }

  uint VersionOPENMP() const {
#ifdef _OPENMP
    return _OPENMP;
#else
    return 0;
#endif
  }
};

template<class TraversalSpecification>
class PostOrderTraversal: public TraversalAlgorithm<TraversalSpecification> {
  
  // un internal class for rethrowing exception caught within parallel sections.
  // This is inspired from the following stackoverflow discussion:
  // https://stackoverflow.com/questions/11828539/elegant-exceptionhandling-in-openmp
  class ThreadExceptionHandler {
    std::exception_ptr ptr_;
    std::mutex         lock_;
  public:
    ThreadExceptionHandler(): ptr_(nullptr) {}
    // Copy initialization (non-thread-safe)
    ThreadExceptionHandler(const ThreadExceptionHandler& other): ptr_(other.ptr_) {}
      
    // If there is an exception it gets rethrowed and the ptr_ is set to nullptr.
    void Rethrow() {
      if(this->ptr_) {
        std::exception_ptr ptr = this->ptr_;
        // reset ptr_ to nullptr so the handler object is cleaned and reusable
        // after the call to Rethrow().
        this->ptr_ = nullptr;
        std::rethrow_exception(ptr);
      }
    }
    void CaptureException() { 
      std::unique_lock<std::mutex> guard(this->lock_);
      this->ptr_ = std::current_exception(); 
    }
    
    template <typename Function, typename... Parameters>
    void Run(Function f, Parameters... params) {
      try {
        f(params...);
      }
      catch (...) {
        CaptureException();
      }
    }
  };
  ThreadExceptionHandler exception_handler_;  
  
public:
  typedef TraversalAlgorithm<TraversalSpecification> ParentType;

  typedef PostOrderMode ModeType;

  PostOrderTraversal(typename TraversalSpecification::TreeType const& tree,
                     TraversalSpecification& spec): ParentType(tree, spec) { }

  void TraverseTree(ModeType mode) {
    switch(mode) {
    case ModeType::SINGLE_THREAD_LOOP_POSTORDER: TraverseTreeSingleThreadLoopPostorder(); break;
    case ModeType::SINGLE_THREAD_LOOP_PRUNES: TraverseTreeSingleThreadLoopPrunes(); break;
    case ModeType::SINGLE_THREAD_LOOP_VISITS: TraverseTreeSingleThreadLoopVisits(); break;
    case ModeType::MULTI_THREAD_LOOP_PRUNES: TraverseTreeMultiThreadLoopPrunes(); break;
    case ModeType::MULTI_THREAD_LOOP_VISITS_THEN_LOOP_PRUNES: TraverseTreeMultiThreadLoopVisitsThenLoopPrunes(); break;
    case ModeType::MULTI_THREAD_LOOP_VISITS: TraverseTreeMultiThreadLoopVisits(); break;
    case ModeType::MULTI_THREAD_VISIT_QUEUE: TraverseTreeMultiThreadVisitQueue(); break;
    case ModeType::HYBRID_LOOP_PRUNES: TraverseTreeHybridLoopPrunes(); break;
    case ModeType::HYBRID_LOOP_VISITS_THEN_LOOP_PRUNES: TraverseTreeHybridLoopVisitsThenLoopPrunes(); break;
    case ModeType::HYBRID_LOOP_VISITS: TraverseTreeHybridLoopVisits(); break;
    default: TraverseTreeAuto();
    }
    exception_handler_.Rethrow();
  }
protected:
  uint current_step_tuning_ = 0;
  uint fastest_step_tuning_ = 0;

  double min_duration_tuning_ = std::numeric_limits<double>::max();
  std::vector<double> durations_tuning_;

  const uvec min_sizes_chunk_ = {8}; //, 4, 8, 16, 32};

  const std::vector<ModeType> choices_mode_auto_ = {
    ModeType::SINGLE_THREAD_LOOP_POSTORDER,
    ModeType::SINGLE_THREAD_LOOP_PRUNES,
    ModeType::SINGLE_THREAD_LOOP_VISITS,
    ModeType::MULTI_THREAD_LOOP_VISITS_THEN_LOOP_PRUNES,
    ModeType::MULTI_THREAD_LOOP_VISITS,
    ModeType::MULTI_THREAD_VISIT_QUEUE
  };

  const std::vector<ModeType> choices_hybrid_mode_auto_ = {
    ModeType::HYBRID_LOOP_PRUNES,
    ModeType::HYBRID_LOOP_VISITS,
    ModeType::HYBRID_LOOP_VISITS_THEN_LOOP_PRUNES
  };

public:
  bool IsTuning() const {
    return current_step_tuning_ < choices_mode_auto_.size() +
      min_sizes_chunk_.size() * choices_hybrid_mode_auto_.size();
  }


  std::string ModeAutoCurrent() const {
    std::ostringstream oss;
    oss<<ModeAuto();
    return oss.str();
  }

  std::string ModeAutoStep(uint step) const {
    std::ostringstream oss;
    oss<<ModeAuto(step);
    return oss.str();
  }

  ModeType ModeAuto() const {
    auto step = IsTuning()? current_step_tuning_ : fastest_step_tuning_;
    return ModeAuto(step);
  }

  ModeType ModeAuto(uint step) const {
    if( step < choices_mode_auto_.size() ) {
      return choices_mode_auto_[step];
    } else {
      uint k = choices_hybrid_mode_auto_.size();
      uint l = step - choices_mode_auto_.size();
      return choices_hybrid_mode_auto_[(l/k) % k];
    }

  }

  uint IndexMinSizeChunkVisit() const {
    auto step = IsTuning()? current_step_tuning_ : fastest_step_tuning_;
    //return (step / min_sizes_chunk_.size()) % min_sizes_chunk_.size();
    return step % min_sizes_chunk_.size();
  }

  uint IndexMinSizeChunkPrune() const {
    auto step = IsTuning()? current_step_tuning_ : fastest_step_tuning_;
    return step % min_sizes_chunk_.size();
  }

  uint min_size_chunk_visit() const {
    return min_sizes_chunk_[IndexMinSizeChunkVisit()];
  }

  uint min_size_chunk_prune() const {
    return min_sizes_chunk_[IndexMinSizeChunkPrune()];
  }

  uint fastest_step_tuning() const {
    return fastest_step_tuning_;
  }

  std::vector<double>  durations_tuning() const {
    return durations_tuning_;
  }

protected:
  void TraverseTreeAuto() {

    std::chrono::steady_clock::time_point start, end;
    double duration;

    ModeType mode = ModeAuto();

    if( IsTuning() ) {

      start = std::chrono::steady_clock::now();
      TraverseTree(mode);
      end = std::chrono::steady_clock::now();

      duration = std::chrono::duration<double, std::milli>(end - start).count();
      durations_tuning_.push_back(duration);
      if(duration < min_duration_tuning_) {
        min_duration_tuning_ = duration;
        fastest_step_tuning_ = current_step_tuning_;
      }
      current_step_tuning_++;

    } else {
      TraverseTree(mode);
    }
  }
  void TraverseTreeSingleThreadLoopPostorder() {
    _PRAGMA_OMP_SIMD
    for(uint i = 0; i < ParentType::ref_tree_.num_nodes(); i++) {
      ParentType::ref_spec_.InitNode(i);
    }

    for(uint i = 0; i < ParentType::ref_tree_.num_nodes() - 1; i++) {
      ParentType::ref_spec_.VisitNode(i);
      ParentType::ref_spec_.PruneNode(i, ParentType::ref_tree_.FindIdOfParent(i));
    }
  }

  void TraverseTreeSingleThreadLoopPrunes() {
    _PRAGMA_OMP_SIMD
    for(uint i = 0; i < ParentType::ref_tree_.num_nodes(); i++) {
      ParentType::ref_spec_.InitNode(i);
    }

    for(uint i_prune = 0;
        i_prune < ParentType::ref_tree_.num_parallel_ranges_prune();
        i_prune++) {
      auto range_prune = ParentType::ref_tree_.RangeIdPruneNode(i_prune);

    _PRAGMA_OMP_SIMD
      for(uint i = range_prune[0]; i <= range_prune[1]; i++) {
        ParentType::ref_spec_.VisitNode(i);
        ParentType::ref_spec_.PruneNode(i, ParentType::ref_tree_.FindIdOfParent(i));
      }
    }
  }

  void TraverseTreeSingleThreadLoopVisits() {
    _PRAGMA_OMP_SIMD
    for(uint i = 0; i < ParentType::ref_tree_.num_nodes(); i++) {
      ParentType::ref_spec_.InitNode(i);
    }

    for(uint i_level = 0; i_level < ParentType::ref_tree_.num_levels(); i_level++) {
      auto range_visit = ParentType::ref_tree_.RangeIdVisitNode(i_level);
      _PRAGMA_OMP_SIMD
      for(uint i = range_visit[0]; i <= range_visit[1]; i++) {
        if(i < ParentType::ref_tree_.num_tips()) {
          // i is a tip (only Visit)
          ParentType::ref_spec_.VisitNode(i);
        } else {
          // i is internal
          for(uint j: ParentType::ref_tree_.FindChildren(i)) {
            ParentType::ref_spec_.PruneNode(j, i);
          }
          ParentType::ref_spec_.VisitNode(i);
        }
      }
    }

    // VisitNode not called on the root node
    for(uint j: ParentType::ref_tree_.FindChildren(ParentType::ref_tree_.num_nodes() - 1)) {
      ParentType::ref_spec_.PruneNode(j, ParentType::ref_tree_.num_nodes() - 1);
    }
  }

  void TraverseTreeMultiThreadLoopVisitsThenLoopPrunes() {

#pragma omp parallel
{
  _PRAGMA_OMP_FOR_SIMD
  for(uint i = 0; i < ParentType::ref_tree_.num_nodes(); i++) {
    exception_handler_.Run([=]{
      ParentType::ref_spec_.InitNode(i);  
    });
  }

  uint i_prune = 0;
  for(uint i_level = 0; i_level < ParentType::ref_tree_.num_levels(); i_level++) {

#pragma omp barrier

    auto range_visit = ParentType::ref_tree_.RangeIdVisitNode(i_level);
    _PRAGMA_OMP_FOR_SIMD
      for(uint i = range_visit[0]; i <= range_visit[1]; i++) {
        exception_handler_.Run([=]{
          ParentType::ref_spec_.VisitNode(i);
        });
      }

      uint num_branches_done = 0;

    while(num_branches_done != range_visit[1] - range_visit[0] + 1) {
#pragma omp barrier
      auto range_prune = ParentType::ref_tree_.RangeIdPruneNode(i_prune);

      _PRAGMA_OMP_FOR_SIMD
        for(uint i = range_prune[0]; i <= range_prune[1]; i++) {
          exception_handler_.Run([=]{
            ParentType::ref_spec_.PruneNode(i, ParentType::ref_tree_.FindIdOfParent(i));
          });
        }

        num_branches_done +=  range_prune[1] - range_prune[0] + 1;
      ++i_prune;
    }
  }
}
  }

  void TraverseTreeMultiThreadLoopVisits() {
#pragma omp parallel
{
  uint tid;
#ifdef _OPENMP
  tid = omp_get_thread_num();
#else
  tid = 0;
#endif

  _PRAGMA_OMP_FOR_SIMD
    for(uint i = 0; i < ParentType::ref_tree_.num_nodes(); i++) {
      exception_handler_.Run([=]{
        ParentType::ref_spec_.InitNode(i);
      });
    }

    for(uint i_level = 0; i_level < ParentType::ref_tree_.num_levels(); i_level++) {
      auto range_visit = ParentType::ref_tree_.RangeIdVisitNode(i_level);
    _PRAGMA_OMP_FOR_SIMD
      for(uint i = range_visit[0]; i <= range_visit[1]; i++) {
        exception_handler_.Run([=]{
          if(i < ParentType::ref_tree_.num_tips()) {
            // i is a tip (only Visit)
            ParentType::ref_spec_.VisitNode(i);
          } else {
            // i is internal
            for(uint j: ParentType::ref_tree_.FindChildren(i)) {
              ParentType::ref_spec_.PruneNode(j, i);
            }
            ParentType::ref_spec_.VisitNode(i);
          }
        });
      }
    }
}
    // VisitNode not called on the root node
    for(uint j: ParentType::ref_tree_.FindChildren(ParentType::ref_tree_.num_nodes() - 1)) {
      ParentType::ref_spec_.PruneNode(j, ParentType::ref_tree_.num_nodes() - 1);
    }
  }

  void TraverseTreeMultiThreadVisitQueue() {
    ParentType::visit_queue_.Init(ParentType::num_children_);
#pragma omp parallel
{
  exception_handler_.Run([=]{
    while(true) {
      uint i = ParentType::visit_queue_.NextInQueue();
      if(i == NA_UINT) {
        continue;
      } else if(i == ParentType::ref_tree_.num_nodes()) {
        break;
      } else if(i < ParentType::ref_tree_.num_tips()) {
        // i is a tip (only Visit)
        ParentType::ref_spec_.InitNode(i);
        ParentType::ref_spec_.VisitNode(i);
        ParentType::visit_queue_.RemoveVisitedNode(i);
      } else if(i < ParentType::ref_tree_.num_nodes() - 1){
        // i is internal
        ParentType::ref_spec_.InitNode(i);
        uvec const& children = ParentType::ref_tree_.FindChildren(i);
        for(uint j: children) {
          ParentType::ref_spec_.PruneNode(j, i);
        }
        ParentType::ref_spec_.VisitNode(i);
        ParentType::visit_queue_.RemoveVisitedNode(i);
      } else {
        // i is the root
        ParentType::ref_spec_.InitNode(i);
        uvec const& children = ParentType::ref_tree_.FindChildren(i);
        for(uint j: children) {
          ParentType::ref_spec_.PruneNode(j, i);
        }
        // don't visit the root
      }
    }
  });
}
}

  void TraverseTreeMultiThreadLoopPrunes() {

#pragma omp parallel
{
  _PRAGMA_OMP_FOR_SIMD
  for(uint i = 0; i < ParentType::ref_tree_.num_nodes(); i++) {
    exception_handler_.Run([=]{
      ParentType::ref_spec_.InitNode(i);
    });
  }

  for(uint i_prune = 0; i_prune < ParentType::ref_tree_.num_parallel_ranges_prune(); i_prune++) {
    auto range_prune = ParentType::ref_tree_.RangeIdPruneNode(i_prune);

    _PRAGMA_OMP_FOR_SIMD
      for(uint i = range_prune[0]; i <= range_prune[1]; i++) {
        exception_handler_.Run([=]{
          ParentType::ref_spec_.VisitNode(i);
          ParentType::ref_spec_.PruneNode(i, ParentType::ref_tree_.FindIdOfParent(i));
        });
      }
  }
}
  }

  void TraverseTreeHybridLoopVisitsThenLoopPrunes() {
    uint min_size_chunk_visit = this->min_size_chunk_visit();
#pragma omp parallel
{
  uint tid;
#ifdef _OPENMP
  tid = omp_get_thread_num();
#else
  tid = 0;
#endif

  _PRAGMA_OMP_FOR_SIMD
    for(uint i = 0; i < ParentType::ref_tree_.num_nodes(); i++) {
      exception_handler_.Run([=]{
        ParentType::ref_spec_.InitNode(i);
      });
    }

    uint i_prune = 0;
  for(uint i_level = 0; i_level < ParentType::ref_tree_.num_levels(); i_level++) {
    auto range_visit = ParentType::ref_tree_.RangeIdVisitNode(i_level);
#pragma omp barrier
    if(range_visit[1] - range_visit[0] + 1 >
         ParentType::num_threads_ * min_size_chunk_visit) {
      _PRAGMA_OMP_FOR_SIMD
        for(uint i = range_visit[0]; i <= range_visit[1]; i++) {
          exception_handler_.Run([=]{
            ParentType::ref_spec_.VisitNode(i);
          });
        }
    } else if(tid == 0) {
      // only the master thread executes this
      _PRAGMA_OMP_SIMD
      for(uint i = range_visit[0]; i <= range_visit[1]; i++) {
        exception_handler_.Run([=]{
          ParentType::ref_spec_.VisitNode(i);
        });
      }
    }

    if (tid == 0) {
      // only one (master) thread executes this
      uint num_branches_done = 0;
      while(num_branches_done != range_visit[1] - range_visit[0] + 1) {
        auto range_prune = ParentType::ref_tree_.RangeIdPruneNode(i_prune);
        _PRAGMA_OMP_SIMD
          for(uint i = range_prune[0]; i <= range_prune[1]; i++) {
            exception_handler_.Run([=]{
              ParentType::ref_spec_.PruneNode(i, ParentType::ref_tree_.FindIdOfParent(i));
            });
          }

          num_branches_done +=  range_prune[1] - range_prune[0] + 1;
        ++i_prune;
      }
    }
  }
}
  }

  void TraverseTreeHybridLoopPrunes() {
    uint min_size_chunk_prune = this->min_size_chunk_prune();
#pragma omp parallel
{
  uint tid;
#ifdef _OPENMP
  tid = omp_get_thread_num();
#else
  tid = 0;
#endif

  _PRAGMA_OMP_FOR_SIMD
    for(uint i = 0; i < ParentType::ref_tree_.num_nodes(); i++) {
      exception_handler_.Run([=]{
        ParentType::ref_spec_.InitNode(i);
      });
    }


  for(uint i_prune = 0; i_prune < ParentType::ref_tree_.num_parallel_ranges_prune(); i_prune++) {
      auto range_prune = ParentType::ref_tree_.RangeIdPruneNode(i_prune);
#pragma omp barrier
      if (range_prune[1] - range_prune[0] + 1 >
            ParentType::num_threads_ * min_size_chunk_prune) {
        _PRAGMA_OMP_FOR_SIMD
        for(uint i = range_prune[0]; i <= range_prune[1]; i++) {
          exception_handler_.Run([=]{
            ParentType::ref_spec_.VisitNode(i);
            ParentType::ref_spec_.PruneNode(i, ParentType::ref_tree_.FindIdOfParent(i));
          });
        }
      } else if (tid == 0) {
        // only one (master) thread executes this
        _PRAGMA_OMP_SIMD
        for(uint i = range_prune[0]; i <= range_prune[1]; i++) {
          exception_handler_.Run([=]{
            ParentType::ref_spec_.VisitNode(i);
            ParentType::ref_spec_.PruneNode(i, ParentType::ref_tree_.FindIdOfParent(i));
          });
        }
      }
    }
}
  }

  void TraverseTreeHybridLoopVisits() {
    uint min_size_chunk_visit = this->min_size_chunk_visit();
#pragma omp parallel
{
  uint tid;
#ifdef _OPENMP
  tid = omp_get_thread_num();
#else
  tid = 0;
#endif

  _PRAGMA_OMP_FOR_SIMD
    for(uint i = 0; i < ParentType::ref_tree_.num_nodes(); i++) {
      exception_handler_.Run([=]{
        ParentType::ref_spec_.InitNode(i);
      });
    }

  for(uint i_level = 0; i_level < ParentType::ref_tree_.num_levels(); i_level++) {
    auto range_visit = ParentType::ref_tree_.RangeIdVisitNode(i_level);
#pragma omp barrier
    if(range_visit[1] - range_visit[0] + 1 >
         ParentType::num_threads_ * min_size_chunk_visit) {
      _PRAGMA_OMP_FOR_SIMD
      for(uint i = range_visit[0]; i <= range_visit[1]; i++) {
        exception_handler_.Run([=]{
          if(i < ParentType::ref_tree_.num_tips()) {
            // i is a tip (only Visit)
            ParentType::ref_spec_.VisitNode(i);
          } else if(i < ParentType::ref_tree_.num_nodes() - 1){
            // i is internal
            for(uint j: ParentType::ref_tree_.FindChildren(i)) {
              ParentType::ref_spec_.PruneNode(j, i);
            }
            ParentType::ref_spec_.VisitNode(i);
          }
        });
      }
    } else if(tid == 0) {
      // only the master thread executes this
      _PRAGMA_OMP_SIMD
      for(uint i = range_visit[0]; i <= range_visit[1]; i++) {
        exception_handler_.Run([=]{
          if(i < ParentType::ref_tree_.num_tips()) {
            // i is a tip (only Visit)
            ParentType::ref_spec_.VisitNode(i);
          } else if(i < ParentType::ref_tree_.num_nodes() - 1){
            // i is internal
            for(uint j: ParentType::ref_tree_.FindChildren(i)) {
              ParentType::ref_spec_.PruneNode(j, i);
            }
            ParentType::ref_spec_.VisitNode(i);
          }
        });
      }
    }
  }
}
    // VisitNode not called on the root
    for(uint j: ParentType::ref_tree_.FindChildren(ParentType::ref_tree_.num_nodes() - 1)) {
      exception_handler_.Run([=]{
        ParentType::ref_spec_.PruneNode(j, ParentType::ref_tree_.num_nodes() - 1);
      });
    }
  }
};


enum PreOrderMode {
  PREORDER_AUTO = 0,
  PREORDER_SINGLE_THREAD_LOOP_PREORDER = 10,
  PREORDER_SINGLE_THREAD_LOOP_VISITS = 12,
  PREORDER_MULTI_THREAD_LOOP_VISITS = 22
};

inline std::ostream& operator<< (std::ostream& os, PreOrderMode mode) {
  switch(mode) {
  case PreOrderMode::PREORDER_AUTO: os<<"PREORDER_AUTO"; break;
  case PreOrderMode::PREORDER_SINGLE_THREAD_LOOP_PREORDER: os<<"PREORDER_SINGLE_THREAD_LOOP_PREORDER"; break;
  case PreOrderMode::PREORDER_SINGLE_THREAD_LOOP_VISITS: os<<"PREORDER_SINGLE_THREAD_LOOP_VISITS"; break;
  case PreOrderMode::PREORDER_MULTI_THREAD_LOOP_VISITS: os<<"PREORDER_MULTI_THREAD_LOOP_VISITS"; break;
  };
  return os<< static_cast<int>(mode);
}

template<class TraversalSpecification>
class PreOrderTraversal: public TraversalAlgorithm<TraversalSpecification> {

  typedef TraversalAlgorithm<TraversalSpecification> ParentType;

public:
  typedef PreOrderMode ModeType;

  PreOrderTraversal(typename TraversalSpecification::TreeType const& tree,
                     TraversalSpecification& spec): ParentType(tree, spec) { }

  void TraverseTree(ModeType mode) {
    switch(mode) {
    case ModeType::PREORDER_SINGLE_THREAD_LOOP_PREORDER: TraverseTreeSingleThreadLoopPreorder(); break;
    case ModeType::PREORDER_SINGLE_THREAD_LOOP_VISITS: TraverseTreeSingleThreadLoopVisits(); break;
    case ModeType::PREORDER_MULTI_THREAD_LOOP_VISITS: TraverseTreeMultiThreadLoopVisits(); break;
    default: TraverseTreeAuto();
    }
  }
protected:
  uint current_step_tuning_ = 0;
  uint fastest_step_tuning_ = 0;

  double min_duration_tuning_ = std::numeric_limits<double>::max();
  std::vector<double> durations_tuning_;

  const uvec min_sizes_chunk_ = {8}; //, 4, 8, 16, 32};

  const std::vector<ModeType> choices_mode_auto_ = {
    ModeType::PREORDER_SINGLE_THREAD_LOOP_PREORDER,
    ModeType::PREORDER_SINGLE_THREAD_LOOP_VISITS,
    ModeType::PREORDER_MULTI_THREAD_LOOP_VISITS
  };

public:
  bool IsTuning() const {
    return current_step_tuning_ < choices_mode_auto_.size();
  }


  std::string ModeAutoCurrent() const {
    std::ostringstream oss;
    oss<<ModeAuto();
    return oss.str();
  }

  std::string ModeAutoStep(uint step) const {
    std::ostringstream oss;
    oss<<ModeAuto(step);
    return oss.str();
  }

  ModeType ModeAuto() const {
    auto step = IsTuning()? current_step_tuning_ : fastest_step_tuning_;
    return ModeAuto(step);
  }

  ModeType ModeAuto(uint step) const {
    return choices_mode_auto_[step%choices_mode_auto_.size()];
  }

  uint fastest_step_tuning() const {
    return fastest_step_tuning_;
  }

  std::vector<double>  durations_tuning() const {
    return durations_tuning_;
  }

protected:
  void TraverseTreeAuto() {

    std::chrono::steady_clock::time_point start, end;
    double duration;

    ModeType mode = ModeAuto();

    if( IsTuning() ) {

      start = std::chrono::steady_clock::now();
      TraverseTree(mode);
      end = std::chrono::steady_clock::now();

      duration = std::chrono::duration<double, std::milli>(end - start).count();
      durations_tuning_.push_back(duration);
      if(duration < min_duration_tuning_) {
        min_duration_tuning_ = duration;
        fastest_step_tuning_ = current_step_tuning_;
      }
      current_step_tuning_++;

    } else {
      TraverseTree(mode);
    }
  }

  void TraverseTreeSingleThreadLoopPreorder() {
    _PRAGMA_OMP_SIMD
    for(uint i = 0; i < ParentType::ref_tree_.num_nodes(); i++) {
      ParentType::ref_spec_.InitNode(i);
    }

    for(uint i = ParentType::ref_tree_.num_nodes() - 1; ; i--) {
      ParentType::ref_spec_.VisitNode(i);
      if(i == 0) {
        break;
      }
    }
  }

  void TraverseTreeSingleThreadLoopVisits() {
    _PRAGMA_OMP_SIMD
    for(uint i = 0; i < ParentType::ref_tree_.num_nodes(); i++) {
      ParentType::ref_spec_.InitNode(i);
    }

    ParentType::ref_spec_.VisitNode(ParentType::ref_tree_.num_nodes() - 1);

    for(uint i_level = ParentType::ref_tree_.num_levels(); i_level > 0; i_level--) {
      auto range_visit = ParentType::ref_tree_.RangeIdVisitNode(i_level - 1);
      _PRAGMA_OMP_SIMD
        for(uint i = range_visit[0]; i <= range_visit[1]; i++) {
          ParentType::ref_spec_.VisitNode(i);
        }
    }
  }

  void TraverseTreeMultiThreadLoopVisits() {
#pragma omp parallel
{
  uint tid;
#ifdef _OPENMP
  tid = omp_get_thread_num();
#else
  tid = 0;
#endif

  _PRAGMA_OMP_FOR_SIMD
    for(uint i = 0; i < ParentType::ref_tree_.num_nodes(); i++) {
      ParentType::ref_spec_.InitNode(i);
    }

    ParentType::ref_spec_.VisitNode(ParentType::ref_tree_.num_nodes() - 1);

    for(uint i_level = ParentType::ref_tree_.num_levels(); i_level > 0; i_level--) {
      auto range_visit = ParentType::ref_tree_.RangeIdVisitNode(i_level - 1);
      _PRAGMA_OMP_FOR_SIMD
        for(uint i = range_visit[0]; i <= range_visit[1]; i++) {
          ParentType::ref_spec_.VisitNode(i);
        }
    }
}
  }

};

template<class TraversalSpecification>
class TraversalTask {
public:
  typedef TraversalSpecification TraversalSpecificationType;
  typedef typename TraversalSpecification::TreeType TreeType;
  typedef typename TraversalSpecification::AlgorithmType AlgorithmType;
  typedef typename AlgorithmType::ModeType ModeType;
  typedef typename TreeType::NodeType NodeType;
  typedef typename TreeType::LengthType LengthType;
  typedef typename TraversalSpecificationType::DataType DataType;
  typedef typename TraversalSpecificationType::ParameterType ParameterType;
  typedef typename TraversalSpecificationType::StateType StateType;

  TraversalTask(
    std::vector<NodeType> const& branch_start_nodes,
    std::vector<NodeType> const& branch_end_nodes,
    std::vector<LengthType> const& branch_lengths,
    DataType const& data):
  tree_(branch_start_nodes, branch_end_nodes, branch_lengths),
  spec_(tree_, data),
  algorithm_(tree_, spec_) {}

  StateType TraverseTree(ParameterType const& par, uint mode) {
    spec_.SetParameter(par);
    algorithm_.TraverseTree(static_cast<ModeType>(mode));
    return spec_.StateAtRoot();
  }

  TreeType & tree() {
    return tree_;
  }
  TraversalSpecification & spec() {
    return spec_;
  }
  AlgorithmType & algorithm() {
    return algorithm_;
  }
  
protected:
  TreeType tree_;
  TraversalSpecification spec_;
  AlgorithmType algorithm_;
};

// A lighter TraversalTask class which gets a reference to an already 
// constructed tree.
template<class TraversalSpecification>
class TraversalTaskLightweight {
public:
  typedef TraversalSpecification TraversalSpecificationType;
  typedef typename TraversalSpecification::TreeType TreeType;
  typedef typename TraversalSpecification::AlgorithmType AlgorithmType;
  typedef typename AlgorithmType::ModeType ModeType;
  typedef typename TreeType::NodeType NodeType;
  typedef typename TreeType::LengthType LengthType;
  typedef typename TraversalSpecificationType::DataType DataType;
  typedef typename TraversalSpecificationType::ParameterType ParameterType;
  typedef typename TraversalSpecificationType::StateType StateType;
  
  TraversalTaskLightweight(
    TreeType const& tree,
    DataType const& data):
    tree_(tree),
    spec_(tree_, data),
    algorithm_(tree_, spec_) {}
  
  StateType TraverseTree(ParameterType const& par, uint mode) {
    spec_.SetParameter(par);
    algorithm_.TraverseTree(static_cast<ModeType>(mode));
    return spec_.StateAtRoot();
  }
  
  TraversalSpecification & spec() {
    return spec_;
  }
  AlgorithmType & algorithm() {
    return algorithm_;
  }
  
protected:
  TreeType const& tree_;
  TraversalSpecification spec_;
  AlgorithmType algorithm_;
};

// The following class defines the main interface of the SPLITT library.
// The user must provide a TraversalSpecificationImplementation class implementing
// this class' methods as described in the comments below. It is
// highly recommended to inherit from this class. However, this is not at all
// obligatory (it is not checked during compilation).
template<class Tree> class TraversalSpecification {
protected:
  // A reference to a Tree available for inheriting classes
  Tree const& ref_tree_;
  // A protected constructor that initializes the tree-reference. This constructor
  // must be called explicitly in the initalization list of inheriting class constructors.
  TraversalSpecification(Tree const& tree): ref_tree_(tree) {}
public:
  // public typedefs. These typedefs must be provided by an implementation class.
  // 1. typedef Tree TreeType;
  // 2. typedef PostOrderTraversal<ImlementationClass> AlgorithmType;
  // 3. typedef ImplementationSpecificParameterType ParameterType;
  // 4. typedef ImplementationSpecificDataType DataType;
  // 5. typedef ImplementationSpecificStateType StateType;


  // The following methods must be present any implementation
  // 6. constructor: will be called by a TraversalTask object; Here, it is
  // commented out, because the DataType is not known.
  // ImplementationClassName(TreeType & tree, DataType & input_data) :
  //   TraversalSpecification(tree) {
  //     implementation specific initialization using the tree and the input_data.
  // }


  // The following methods get called by the TraversalAlgorithm implementation:

  // 7. Setting the model parameters prior to starting the pruning procedure on the tree.
  // This method is called by the TraversalTask.TraverseTree(ParamterType const&, uint mode)
  // method. The method declaration is commented out because ParameterType is not known
  // and must be specified by the implementing class.
  // void SetParameter(ParameterType const& par);

  // 8. InitNode(i) is called on each node in the tree right after SetParameter and
  // before any of the VisitNode and PruneNode has been called. There is no predefined
  // order of the calls to InitNode and they may be executed in parallel. Therefore, only
  // node-specific data initialization, including the length of the branch
  // leading to node i, can take place in this method.
  void InitNode(uint i) {}


  // 9. VisitNode(i) is called on each tip or internal node (EXCLUDING THE ROOT),
  // in the tree after PruneNode has been called on each descendant of i.
  // The method is the perfect place to calculate the state of node i using the
  // pre-calculated states of its descendants. Although, it is guaranteed
  // that VisitNode(i) is called before VisitNode(i_parent), this method SHOULD NOT BE USED
  // FOR ALTERING THE STATE of i_parent, because this would conflict with
  // a concurrent execution of VisitNode on a sibling of i (see also PruneNode).
  void VisitNode(uint i) {}

  // 10. PruneNode(i, i_parent) is called on each tip or internal node (EXCLUDING THE ROOT)
  // after VisitNode(i) and in sync with PruneNode(k, i_parent), for any sibling k of i.
  // Thus, it is safe to use PruneNode to update the state of i_parent.
  void PruneNode(uint i, uint i_parent) {}

  // 11. StateType StateAtRoot() is called after PruneNode has been called on each
  // direct descendant of the root node. If necessary, VisitNode(i_root) can be called
  // here, in order to calculate the final state of the root. The value returned by this
  // function is also returned by the TraversalTask.TraverseTree(ParameterType const& par, uint mode)
  // method.
};

// 12. After the class TraversalSpecificationImplementation has been defined it is
// time to specify the TraversalTask template. This is not obligatory but can be very
// convinient for creating TraversalTask objects with the user specific implementation
// and to call their TraverseTree method.
// typedef TraversalTask<TraversalSpecificationImplementation> > MyTraversalTask;
}
#endif // SPLITT_SPLITT_H_
