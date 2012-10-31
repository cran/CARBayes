\name{CARBayes-package}
\alias{CARBayes-package}
\alias{CARBayes}
\docType{package}


\title{
Spatial areal unit modelling
}

\description{
This package implements Bayesian hierarchical spatial areal unit models. In such models the spatial correlation is modelled by  a set of random effects, which are assigned a conditional autoregressive (CAR) prior distribution. Examples of the models included are the BYM model as well as recently developed localised spatial smoothing models. In addition, Moran's I and Geary's C permutation tests for spatial correlation can be implemented.

The creation of this package was supported by the Economic and Social Research Council (ESRC) grant RES-000-22-4256.
}

\details{
\tabular{ll}{
Package: \tab CARBayes\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2012-10-31\cr
License: \tab GPL (>= 2)\cr
}

The package implements Bayesian hierarchical models for spatial areal unit data. The models are an extension of generalised linear models, and the response variable can be binomial, Gaussian or Poisson. The mean function for the data comprises a covariate component and a set of random effects. This package implements a number of different random effects models, which represent different types of spatial correlation structure. These models are within the conditional autoregressive (CAR) class, and capture spatial correlation via a binary neighbourhood matrix W. This matrix has a jkth element equal to one and zero, depending on whether areas j and k are spatially close (typically chosen to mean sharing a common border). Inference in each case is via Markov Chain Monte Carlo (MCMC) simulation. In addition, permutation tests based on Geary's C and Moran's I statistics for spatial correlation is also available. Finally, there is a function that can elicit a matrix of prior distributions for the spatial structure of the random effects based n Geary's C or Moran's I ordinates. The following random effects models are implemented in this package.

Independent and identically distributed random effects.

Spatially correlated random effects modelled by the intrinsic autoregressive (IAR) model proposed by Besag et al (1991).

Spatially correlated random effects modelled by the Besag-York-Mollie (BYM) model proposed by Besag et al (1991).

Spatially correlated random effects modelled by the proper CAR model proposed by Stern and Cressie (1999).

Spatially correlated random effects modelled by the proper CAR model proposed by Leroux (1999).

Spatially correlated random effects modelled by the localised spatial smoothing approaches proposed by Lee and Mitchell (2012).

Spatially correlated random effects modelled by the localised spatial smoothing approaches proposed by Lee (2012).

Further details about the model and how to implement them are given in the vignette accompanying this package. More detailed descriptions of the methods are given by the references below.
}


\author{
Maintainer: Duncan Lee <Duncan.Lee@glasgow.ac.uk>
}

\references{
Besag, J., J. York, and A. Mollie (1991). Bayesian image restoration with two applications in spatial statistics. Annals of the Institute of Statistics and Mathematics 43, 1-59.

Stern, H and Cressie, N. (1999). Inference for extremes in disease mapping, Chapter
Disease mapping and Risk Assessment for Public Health. Lawson, A and Biggeri, D
and Boehning, E and Lesaffre, E and Viel, J and Bertollini, R (eds). Wiley.

Lee, D. (2011). A comparison of conditional autoregressive model used in Bayesian disease mapping. Spatial and Spatio-temporal Epidemiology 2, 79-89.

Lee, D. and R. Mitchell (2012). Boundary detection in disease mapping studies. Biostatistics, 13, 415-426.

Lee, D. (2012). Using prior information to identify boundaries in disease risk maps. arXiv:1208.4944v1.

Leroux, B., X. Lei, and N. Breslow (1999). Estimation of disease rates in small areas: A new mixed model for spatial dependence, Chapter Statistical Models in Epidemiology, the Environment and Clinical Trials, Halloran, M and Berry, D (eds), pp. 135-178. Springer-Verlag, New York.
}


\examples{
##################################################
#### Run the models on simulated data on a lattice
##################################################

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
	
	
#### Generate the covariates and response data
x1 <- rnorm(n)
x2 <- rnorm(n)
E <- rep(40,n)
theta <- rnorm(n, sd=0.05)
phi <- mvrnorm(n=1, mu=rep(0,n), Sigma=0.4 * exp(-0.1 * distance))
risk <- exp(-0.2 +  0.1 * x1 + 0.1*x2 + theta + phi)
fitted <- E * risk
Y <- rpois(n=n, lambda=fitted)



#### Run the BYM model for Poisson data
#### Let the function randomly generate starting values for the parameters
#### Use the default priors specified by the function (for details see the help files)
formula <- Y ~ x1 + x2 + offset(log(E))
\dontrun{model <- poisson.bymCAR(formula=formula, W=W, burnin=2000, n.sample=5000)}
\dontshow{model <- poisson.bymCAR(formula=formula, W=W, burnin=20, n.sample=50)}
}
