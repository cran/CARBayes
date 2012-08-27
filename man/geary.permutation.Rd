\name{geary.permutation}
\alias{geary.permutation}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Geary's C permutation test
}

\description{
Implement a permutation test for spatial correlation based on Geary's C statistic.
}

\usage{
geary.permutation(x, W, n.simulation)
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
Geary's C test statistic is returned, as well as a p-value for a one-sided test whose null hypothesis is that there is no positive spatial correlation.
}


\value{
\item{statistic }{The Geary's C statistic.}
\item{rank }{The rank of the Geary's C statistic compared to the set of values from the permuted data sets.}
\item{value }{The p-value from the permutation test.}
}


\references{
Geary, R. C. (1954). The Contiguity Ratio and Statistical Mapping. The Incorporated Statistician 5, 115???145.
}

\author{
Duncan Lee
}


\examples{
###################################################
#### Apply Moran's I to simulated data on a lattice
###################################################
#### Set up a square lattice region
x.easting <- 1:10
x.northing <- 1:10
Grid <- expand.grid(x.easting, x.northing)
n <- nrow(Grid)

#### set up distance and neighbourhood (W, based on sharing a common border) matrices
distance <-array(0, c(n,n))
W <-array(0, c(n,n))
  for(i in 1:n)
  {
     	for(j in 1:n)
		{
		temp <- (Grid[i,1] - Grid[j,1])^2 + (Grid[i,2] - Grid[j,2])^2
		distance[i,j] <- sqrt(temp)
			if(temp==1)  W[i,j] <- 1 
		}	
	}
	
	
#### Generate the data
phi <- mvrnorm(n=1, mu=rep(0,n), Sigma=0.4 * exp(-0.1 * distance))
Y <- rnorm(n=n, mean=phi, sd=rep(0.1,n))


#### Compute Geary's C
geary.permutation(x=Y, W=W, n.simulation=500)
}