% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/likPOUMM.R
\name{pruneTree}
\alias{pruneTree}
\title{Extract information for fast likelihood calculation using the
 breadth-first pruning algorithm.}
\usage{
pruneTree(tree, z, se = 0)
}
\arguments{
\item{tree}{a phylo object}

\item{z}{Numeric vector with length(tree$tip.label) values 
corresponding to tree$tip.label.}

\item{se}{Non-negative numerical or N-vector indicating known
standard measurement error.}
}
\value{
a list of objects used for likelihood evaluation
}
\description{
Extract information for fast likelihood calculation using the
 breadth-first pruning algorithm.
}
