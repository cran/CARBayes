\name{moran.permutation}
\alias{moran.permutation}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Moran's I permutation test
}

\description{
Implement a permutation test for spatial correlation based on Moran's I statistic.
}

\usage{
moran.permutation(x, W, n.simulation)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A vector of observations relating to n spatial units that you want to test for spatial correlation.
}
  \item{W}{
A binary n by n neighbourhood matrix (where n is the number of spatial units). The jkth element equals one if areas (j, k) are spatially close (e.g. share a common border) and is zero otherwise. 
}
  \item{n.simulation}{
The number of random permutations for the test.
}
}


\details{
Moran's test statistic is returned, as well as a p-value for a one-sided test whose null hypothesis is that there is no positive spatial correlation.
}


\value{
\item{statistic }{The Moran's I statistic.}
\item{rank }{The rank of the Moran's I statistic compared to the set of values from the permuted data sets.}
\item{value }{The p-value from the permutation test.}
}


\references{
Moran, P. (1950). Notes on continuous stochastic phenomena. Biometrika 37, 17-23.}

\author{
Duncan Lee}