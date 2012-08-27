\name{binomial.priorWCAR}
\alias{binomial.priorWCAR}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Fit a localised conditional autoregressive (CAR) model to spatial binomial data where the neighbourhood matrix is treated as a collection of binary random quantities
}

\description{
The function fits a binomial logistic random effects models to spatial data, where the random effects are modelled by a localised conditional autoregressive (CAR) model, which is a generalisation of that proposed by Lee (2012). The random effects in neighbouring areas (e.g. those that share a common border) are modelled as correlated or conditionally independent, with this choice being determined by the data via a random neighbourhood matrix. The model represents the logit of the mean function for the set of binomial responses  by a combination of covariates and a set of random effects. Inference is based on Markov Chain Monte Carlo (MCMC) simulation, using a combination of Gibbs sampling and Metropolis steps.
}

\usage{
binomial.priorWCAR(formula, beta = NULL, phi = NULL, tau2 = NULL, 
rho = NULL, fix.rho = FALSE, trials, W, blocksize.W = 2, prior.W = NULL, 
burnin = 0, n.sample = 1000, blocksize.beta = 5, blocksize.phi = 10, 
prior.mean.beta = NULL, prior.var.beta = NULL, prior.max.tau2 = NULL)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
A formula for the covariate part of the model, using the same notation as for the lm() function. The offsets should also be included here using the offset() function.
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
A starting value for the global spatial correlation parameter rho. If this argument is not specified the function will randomly generate a starting value unless fix.rho=TRUE, in which case it is required.
}
  \item{fix.rho}{
Logical, whether the global correlation parameter is fixed or estimated in the MCMC algorithm.
}
 \item{trials}{
A vector the same length as the response containing the total number of trials for each area.
}
  \item{W}{
A binary n by n neighbourhood matrix (where n is the number of spatial units). The jkth element equals one if areas (j, k) are spatially close (e.g. share a common border) and is zero otherwise. 
}
  \item{blocksize.W}{
The number of elements of W to update in one block at each MCMC iteration.
}
  \item{prior.W}{
An n by n matrix specifying the prior probabilities that each element of W equals one. Note, for non-neighbouring areas the elements of this matrix should be zero.
}
  \item{burnin}{
The number of MCMC samples to discard as the burnin period. Defaults to 0.
}
  \item{n.sample}{
The number of MCMC samples to generate. Defaults to 1,000.
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
For further details about how to apply the function see the vignette to the CARBayes package.
}


\value{
\item{formula }{The formula for the covariate and offset part of the model.}
\item{samples.beta }{A matrix of MCMC samples for the regression parameters beta.}
\item{samples.phi }{A matrix of MCMC samples for the random effects phi.}
\item{samples.tau2 }{A matrix of MCMC samples for the random effects variance tau2.}
\item{samples.rho }{A matrix of MCMC samples for the global spatial correlation parameter rho. Only included if fix.rho=FALSE.}
\item{samples.W }{A matrix of MCMC samples for the elements of the neighbourhood matrix W.}
\item{fitted.values }{A summary matrix of the posterior distributions of the fitted values for each area. The summaries include: Mean, Sd, Median, and credible interval.}
\item{random.effects }{A summary matrix of the posterior distributions of the random effects for each area. The summaries include: Mean, Sd, Median, and credible interval.}
\item{W.posterior }{The posterior median for the n by n neighbourhood matrix W, which shows the locations of any boundaries in the risk surface.}
\item{W.border.prob }{An n by n matrix showing the posterior probability that there is a boundary in the risk surface between each pair of neighbours in the study region.}
\item{residuals }{A vector of raw residuals from the model.}
\item{DIC }{The Deviance Information Criterion.}
\item{p.d }{The effective number of parameters in the model.}
\item{summary.results }{A summary table of the parameters.}
\item{accept.W }{The acceptance probability for the elements of W.}
}


\references{
Lee, D. (2012). Using prior information to identify boundaries in disease risk maps. arXiv:1208.4944v1.
}

\author{
Duncan Lee
}


\examples{
###########################################################
#### Run the model on simulated data - localised CAR model
###########################################################
#### Set up a square lattice region
x.easting <- 1:10
x.northing <- 1:10
Grid <- expand.grid(x.easting, x.northing)
n <- nrow(Grid)

#### Split the area into two groups between which there will be a boundary.
groups <-rep(1, n) 
groups[Grid$Var1>5] <- 2

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
		
#### Generate the response data
phi <- mvrnorm(n=1, mu=groups, Sigma=0.2 * exp(-0.1 * distance))
logit <- phi
prob <- exp(logit) / (1 + exp(logit))
trials <- rep(50,n)
Y <- rbinom(n=n, size=trials, prob=prob)

#### Run the localised smoothing model
#### Let the function randomly generate starting values for the parameters
#### Use the default priors specified by the function (for details see the help files)
#### Specify a weakly informative prior for W that P(wkj=1)=0.5.
formula <- Y ~ 1
\dontrun{model <- binomial.priorWCAR(formula=formula, rho=0.99, fix.rho=TRUE, 
trials=trials, W=W, blocksize.W=2, prior.W=W/2, burnin=5000, n.sample=10000)}
\dontshow{model <- binomial.priorWCAR(formula=formula, rho=0.99, fix.rho=TRUE, 
trials=trials, W=W, blocksize.W=2, prior.W=W/2, burnin=20, n.sample=50)}
}