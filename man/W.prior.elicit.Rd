\name{W.prior.elicit}
\alias{W.prior.elicit}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Elicit a prior for the binary neighbourhood matrix W
}

\description{
This function elicits a prior for the binary neighbourhood matrix W. For pairs of areas that do not share a common border the corresponding element of W is fixed at zero. For pairs of areas that do share a common border the prior probability that the corresponding element of W is equal to one is elicited. The elicitation is based on either Geary's C or Moran's I ordinates, and details are given in the reference listed below. The output of this function is a matrix of prior probabilities that each element of W equals one, and is used as input for the functions binomial.priorWCAR, gaussian.priorWCAR and poisson.priorWCAR.
}

\usage{
W.prior.elicit(phistar, W, ordinate)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phistar}{
A vector of observations relating to the n spatial units that you want to elicit prior probabilities about the spatial structure.
}
  \item{W}{
A binary n by n neighbourhood matrix (where n is the number of spatial units). The jkth element equals one if areas (j, k) are spatially close (e.g. share a common border) and is zero otherwise.
}
  \item{ordinate}{
Either "geary" or "moran" depending on which ordinates you wish to use for the prior elicitation.
}
}

\details{
A matrix of prior probabilities is return for the binary neighbourhood matrix W.
}

\value{
\item{W.prior }{The matrix of prior probabilities for W.}
}

\references{
Lee, D. (2012). Using prior information to identify boundaries in disease risk maps. arXiv:1208.4944v1.
}

\author{
Duncan Lee
}

\examples{
#########################################################################
#### Apply the prior elicitation procedure to simulated data on a lattice
#########################################################################
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


#### Apply the prior elicitation procedure
W.prior <- W.prior.elicit(phistar=Y, W=W, ordinate="geary")
}