#include <RcppArmadillo.h>
#include <R_ext/Rdynload.h>

#include "AbcPOUMM.h"

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

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

using namespace splittree;

ParallelPruningAbcPOUMM* CreateParallelPruningAbcPOUMM(
    Rcpp::List const& tree, vec const& z, vec const& se) {
  arma::umat branches = tree["edge"];
  uvec br_0 = arma::conv_to<uvec>::from(branches.col(0));
  uvec br_1 = arma::conv_to<uvec>::from(branches.col(1));
  vec t = Rcpp::as<vec>(tree["edge.length"]);
  uint num_tips = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size();
  uvec node_names = Seq(uint(1), num_tips);
  typename ParallelPruningAbcPOUMM::DataType data(node_names, z, se);
  return new ParallelPruningAbcPOUMM(br_0, br_1, t, data);
}

RCPP_EXPOSED_CLASS_NODECL(ParallelPruningAbcPOUMM::TreeType)
RCPP_EXPOSED_CLASS_NODECL(ParallelPruningAbcPOUMM::TraversalSpecificationType)
RCPP_EXPOSED_CLASS_NODECL(ParallelPruningAbcPOUMM::AlgorithmType)

RCPP_MODULE(POUMM_AbcPOUMM) {
  Rcpp::class_<ParallelPruningAbcPOUMM::TreeType::Tree> ( "POUMM_Tree" )
  .property("num_nodes", &ParallelPruningAbcPOUMM::TreeType::Tree::num_nodes )
  .property("num_tips", &ParallelPruningAbcPOUMM::TreeType::Tree::num_tips )
  .property("map_id_to_node", &ParallelPruningAbcPOUMM::TreeType::Tree::map_id_to_node )
  .method("LengthOfBranch", &ParallelPruningAbcPOUMM::TreeType::Tree::LengthOfBranch )
  .method("FindNodeWithId", &ParallelPruningAbcPOUMM::TreeType::Tree::FindNodeWithId )
  .method("FindIdOfNode", &ParallelPruningAbcPOUMM::TreeType::Tree::FindIdOfNode )
  .method("FindIdOfParent", &ParallelPruningAbcPOUMM::TreeType::Tree::FindIdOfParent )
  .method("OrderNodes", &ParallelPruningAbcPOUMM::TreeType::Tree::OrderNodes )
  ;
  Rcpp::class_<ParallelPruningAbcPOUMM::TreeType>( "POUMM_OrderedTree" )
    .derives<ParallelPruningAbcPOUMM::TreeType::Tree> ( "POUMM_Tree" )
    .method("RangeIdPruneNode", &ParallelPruningAbcPOUMM::TreeType::RangeIdPruneNode )
    .method("RangeIdVisitNode", &ParallelPruningAbcPOUMM::TreeType::RangeIdVisitNode )
    .property("num_levels", &ParallelPruningAbcPOUMM::TreeType::num_levels )
    .property("num_parallel_ranges_prune", &ParallelPruningAbcPOUMM::TreeType::num_parallel_ranges_prune )
    .property("ranges_id_visit", &ParallelPruningAbcPOUMM::TreeType::ranges_id_visit )
    .property("ranges_id_prune", &ParallelPruningAbcPOUMM::TreeType::ranges_id_prune )
  ;
  Rcpp::class_<ParallelPruningAbcPOUMM::AlgorithmType::ParentType>( "POUMM_TraversalAlgorithm" )
    .property( "VersionOPENMP", &ParallelPruningAbcPOUMM::AlgorithmType::ParentType::VersionOPENMP )
    .property( "num_threads", &ParallelPruningAbcPOUMM::AlgorithmType::num_threads )
  ;
  Rcpp::class_<ParallelPruningAbcPOUMM::AlgorithmType> ( "POUMM_ParallelPruning" )
    .derives<ParallelPruningAbcPOUMM::AlgorithmType::ParentType>( "POUMM_TraversalAlgorithm" )
    .method( "ModeAutoStep", &ParallelPruningAbcPOUMM::AlgorithmType::ModeAutoStep )
    .property( "ModeAutoCurrent", &ParallelPruningAbcPOUMM::AlgorithmType::ModeAutoCurrent )
    .property( "IsTuning", &ParallelPruningAbcPOUMM::AlgorithmType::IsTuning )
    .property( "min_size_chunk_visit", &ParallelPruningAbcPOUMM::AlgorithmType::min_size_chunk_visit )
    .property( "min_size_chunk_prune", &ParallelPruningAbcPOUMM::AlgorithmType::min_size_chunk_prune )
    .property( "durations_tuning", &ParallelPruningAbcPOUMM::AlgorithmType::durations_tuning )
    .property( "fastest_step_tuning", &ParallelPruningAbcPOUMM::AlgorithmType::fastest_step_tuning )
  ;
  Rcpp::class_<ParallelPruningAbcPOUMM::TraversalSpecificationType> ( "POUMM_PruningSpec" )
  .field( "a", &ParallelPruningAbcPOUMM::TraversalSpecificationType::a )
  .field( "b", &ParallelPruningAbcPOUMM::TraversalSpecificationType::b )
  .field( "c", &ParallelPruningAbcPOUMM::TraversalSpecificationType::c )
  .field( "se", &ParallelPruningAbcPOUMM::TraversalSpecificationType::se )
  .field( "z", &ParallelPruningAbcPOUMM::TraversalSpecificationType::z )
  .field( "alpha", &ParallelPruningAbcPOUMM::TraversalSpecificationType::alpha )
  .field( "sigma2", &ParallelPruningAbcPOUMM::TraversalSpecificationType::sigma2 )
  .field( "sigmae2", &ParallelPruningAbcPOUMM::TraversalSpecificationType::sigmae2 )
  .field( "theta", &ParallelPruningAbcPOUMM::TraversalSpecificationType::theta )
  ;
  Rcpp::class_<ParallelPruningAbcPOUMM>( "POUMM_AbcPOUMM" )
  .factory<Rcpp::List const&, vec const&, vec const&>(&CreateParallelPruningAbcPOUMM)
  .method( "DoPruning", &ParallelPruningAbcPOUMM::TraverseTree )
  .property( "tree", &ParallelPruningAbcPOUMM::tree )
  .property( "spec", &ParallelPruningAbcPOUMM::spec )
  .property( "algorithm", &ParallelPruningAbcPOUMM::algorithm )
  ;
}