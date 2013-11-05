\name{poisson.properCAR}
\alias{poisson.properCAR}
%- Also NEED an '\alias' for EACH other topic documented here.


\title{
Fit the proper  conditional autoregressive (CAR) model to spatial Poisson data
}


\description{
The function fits a Poisson log-normal random effects models to spatial count data, where the random effects are modelled by the conditional autoregressive (CAR) model proposed by Stern and Cressie (1999). The model represents the natural log of the mean function for the set of Poisson responses  by a combination of covariates and a set of random effects. The latter are spatially correlated and come from the Stern and Cressie CAR model. A set of offsets can also be included on the linear predictor scale. Inference is based on Markov Chain Monte Carlo (MCMC) simulation, using a combination of Gibbs sampling and Metropolis steps.
}


\usage{
poisson.properCAR(formula, data = NULL, beta = NULL, phi = NULL, tau2 = NULL, rho = NULL,
 W, burnin = 0, n.sample = 1000, thin=1, blocksize.beta = 5, blocksize.phi = 10, 
prior.mean.beta = NULL, prior.var.beta = NULL, prior.max.tau2 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{formula}{
A formula for the covariate part of the model, using the same notation as for the lm() function. The offsets should also be included here using the offset() function.
}
  \item{data}{
A data.frame containing the  variables in the formula.
}
  \item{beta}{
A vector of starting values for the regression parameters (including the intercept term). If this argument is not specified the function will randomly generate starting values.
}
  \item{phi}{
A vector of starting values for the random effects. If this argument is not specified the function will randomly generate starting values.
}
  \item{tau2}{
A starting value for the variance parameter of the random effects. If this argument is not specified the function will randomly generate a starting value.
}
  \item{rho}{
A starting value for the spatial correlation parameter rho. If this argument is not specified the function will randomly generate a starting value.
}
  \item{W}{
A binary n by n neighbourhood matrix (where n is the number of spatial units). The jkth element equals one if areas (j, k) are spatially close (e.g. share a common border) and is zero otherwise. 
}
  \item{burnin}{
The number of MCMC samples to discard as the burnin period. Defaults to 0.
}
  \item{n.sample}{
The number of MCMC samples to generate. Defaults to 1,000.
}
  \item{thin}{
The level of thinning to apply to the MCMC samples to reduce their temporal autocorrelation. Defaults to 1.
}
  \item{blocksize.beta}{
The size of the blocks in which to update the regression parameters in the MCMC algorithm. Defaults to 5.
}
  \item{blocksize.phi}{
The size of the blocks in which to update the random effects in the MCMC algorithm. Defaults to 10.
}
  \item{prior.mean.beta}{
A vector of prior means for the regression parameters beta (Gaussian priors are assumed). Defaults to a vector of zeros.
}
  \item{prior.var.beta}{
A vector of prior variances for the regression parameters beta (Gaussian priors are assumed). Defaults to a vector with values 1000.
}
  \item{prior.max.tau2}{
The maximum allowable value for the random effects variance tau2 (a Uniform(0,M) prior is assumed). Defaults to M=1000.
}
}

\details{
For further details about how to apply the function see the examples below and in the main CARBayes helpfile.
}


\value{
\item{formula }{The formula for the covariate and offset part of the model.}
\item{samples }{A list containing the MCMC samples from the model.}
\item{fitted.values }{A summary matrix of the posterior distributions of the fitted values for each area. The summaries include: Mean, Sd, Median, and credible interval.}
\item{random.effects }{A summary matrix of the posterior distributions of the random effects for each area. The summaries include: Mean, Sd, Median, and credible interval.}
\item{residuals }{A summary matrix of the posterior distributions of the residuals for each area. The summaries include: Mean, Sd, Median, and credible interval.}
\item{W.summary }{The neighbourhood matrix W from the model.}
\item{DIC }{The Deviance Information Criterion.}
\item{p.d }{The effective number of parameters in the model.}
\item{summary.results }{A summary table of the parameters.}
\item{model }{A text string describing the model fit.}
\item{accept }{The acceptance probabilities for the parameters.}
}


\references{
Stern, H. and N. Cressie (1999). Inference for extremes in disease mapping, Chapter
Disease mapping and Risk Assessment for Public Health. Lawson, A and Biggeri, D
and Boehning, E and Lesaffre, E and Viel, J and Bertollini, R (eds). Wiley.
}


\author{
Duncan Lee
}


\examples{
##################################################
#### Run the model on simulated data on a lattice
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



#### Run the Proper model
#### Let the function randomly generate starting values for the parameters
#### Use the default priors specified by the function (for details see the help files)
formula <- Y ~ x1 + x2 + offset(log(E))
\dontrun{model <- poisson.properCAR(formula=formula, W=W, burnin=5000, n.sample=10000)}
\dontshow{model <- poisson.properCAR(formula=formula, W=W, burnin=20, n.sample=50)}
}
