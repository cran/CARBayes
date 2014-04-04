\name{bymCAR.re}
\alias{bymCAR.re}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Fit the BYM  conditional autoregressive (CAR) model to spatial data
}

\description{
The function fits a Bayesian hierarchical model with spatially correlated random effects to the data, whre the data likelihood can be binomial, Gaussian or Poisson. The random effects are modelled by the BYM conditional autoregressive (CAR) model (Besag et. al. 1991). The model represents the linear predictor of the data  by a combination of covariates and two sets of random effects. For the latter, the first set are independent, while the second are spatially correlated and come from the IAR model. A set of offsets can also be included on the linear predictor scale. Inference is based on Markov Chain Monte Carlo (MCMC) simulation, using a combination of Gibbs sampling and Metropolis steps.
}

\usage{
bymCAR.re(formula, family=NULL, data=NULL, trials=NULL, W, burnin=0, n.sample=1000, 
thin=1, blocksize.beta=5, prior.mean.beta=NULL, prior.var.beta=NULL, prior.nu2=NULL, 
prior.tau2=NULL, prior.sigma2=NULL, verbose=TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
A formula for the covariate part of the model, using the same notation as for the lm() function. The offsets should also be included here using the offset() function.
}
  \item{family}{
One of either 'binomial', 'gaussian' or 'poisson', which respectively specify a binomial likelihood model with a logistic link function, a Gaussian likelihood model with an identity link function, or a Poisson likelihood model with a log link function. 
}
  \item{data}{
A data.frame containing the  variables in the formula.
}
  \item{trials}{
Only used if family='binomial'. A vector the same length as the response containing the total number of trials for each area.
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
  \item{prior.mean.beta}{
A vector of prior means for the regression parameters beta (Gaussian priors are assumed). Defaults to a vector of zeros.
}
  \item{prior.var.beta}{
A vector of prior variances for the regression parameters beta (Gaussian priors are assumed). Defaults to a vector with values 1000.
}
  \item{prior.nu2}{
Only used for the Gaussian model. The prior shape and scale in the form of c(shape, scale) for an Inverse-Gamma(shape, scale) prior for nu2. Defaults to c(0.001, 0.001).  
}
     \item{prior.tau2}{
The prior shape and scale in the form of c(shape, scale) for an Inverse-Gamma(shape, scale) prior for tau2. Defaults to c(0.001, 0.001).  
}
  \item{prior.sigma2}{
The prior shape and scale in the form of c(shape, scale) for an Inverse-Gamma(shape, scale) prior for sigma2. Defaults to c(0.001, 0.001).  
}
  \item{verbose}{
Logical, should the function update the user on its progress.  
}
}



\details{
For further details about how to apply the function see the examples below and in the vignette.
}


\value{
\item{formula }{The formula for the covariate and offset part of the model.}
\item{samples }{A list containing the MCMC samples from the model.}
\item{fitted.values }{A summary matrix of the posterior distributions of the fitted values for each area. The summaries include: Mean, Sd, Median, and credible interval.}
\item{residuals }{A summary matrix of the posterior distributions of the residuals for each area. The summaries include: Mean, Sd, Median, and credible interval.}
\item{W.summary }{The neighbourhood matrix W from the model.}
\item{modelfit }{Model fit criteria including the Deviance Information Criterion (DIC), the effective number of parameters in the model(p.d), DIC3 and the Marginal Predictive Likelihood (MPL).}
\item{summary.results }{A summary table of the parameters.}
\item{model }{A text string describing the model fit.}
\item{accept }{The acceptance probabilities for the parameters.}
}

\references{
Besag, J., J. York, and A. Mollie (1991). Bayesian image restoration with two applications in spatial statistics. Annals of the Institute of Statistics and Mathematics 43, 1-59.
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
theta <- rnorm(n, sd=0.05)
phi <- mvrnorm(n=1, mu=rep(0,n), Sigma=0.4 * exp(-0.1 * distance))
logit <- x1 + x2 + theta + phi
prob <- exp(logit) / (1 + exp(logit))
trials <- rep(50,n)
Y <- rbinom(n=n, size=trials, prob=prob)


#### Run the BYM model
formula <- Y ~ x1 + x2
\dontrun{model <- bymCAR.re(formula=formula, family="binomial, trials=trials,
W=W, burnin=5000, n.sample=10000)}
}