% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/paramsPOUMM.R
\name{H2}
\alias{H2}
\title{Phylogenetic heritability estimated at time t}
\usage{
H2(alpha, sigma, sigmae, t = Inf, tm = 0)
}
\arguments{
\item{alpha, sigma}{numeric, parameters of the OU process acting on the
genetic contributions}

\item{sigmae}{numeric, environmental standard deviation}

\item{t}{time from the beginning of the process at which heritability should
be calculated, i.e. epidemiologic time}

\item{tm}{average time for within host evolution from getting infected until
getting measured or passing on the infection to another host}
}
\description{
Phylogenetic heritability estimated at time t
}
