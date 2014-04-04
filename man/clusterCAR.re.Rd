\name{clusterCAR.re}
\alias{clusterCAR.re}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Fit a cluster model with spatially correlated random effects to spatial data.
}

\description{
The function fits a Bayesian hierarchical model with spatially correlated random effects and a cluster component to the data, whre the data likelihood can be binomial, Gaussian or Poisson. The random effects are modelled by the conditional autoregressive (CAR) model proposed by Leroux et. al. (1999). The model represents the linear predictor of the data  by a cluster component and a set of random effects. The latter are spatially correlated and come from the Leroux CAR model. Inference is based on Markov Chain Monte Carlo (MCMC) simulation, using a combination of Gibbs sampling and Metropolis steps.
}

\usage{
clusterCAR.re(Y, q, family=NULL, expected=NULL, trials=NULL, W, burnin=0, 
n.sample=1000, thin=1, prior.nu2=NULL, prior.tau2=NULL, prior.rho=NULL, 
verbose=TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
A vector of data representing the response variable to be modelled.
}
  \item{q}{
The number of clusters to fit to the data. Must be an interger with minimum value q=2.
}
  \item{family}{
One of either 'binomial', 'gaussian' or 'poisson', which respectively specify a binomial likelihood model with a logistic link function, a Gaussian likelihood model with an identity link function, or a Poisson likelihood model with a log link function. 
}
  \item{expected}{
Only used if family='poisson'. A vector of expected counts used as an offset.
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
  \item{prior.nu2}{
The prior shape and scale in the form of c(shape, scale) for an Inverse-Gamma(shape, scale) prior for nu2. Only for the Gaussian model. Defaults to c(0.001, 0.001).  
}
     \item{prior.tau2}{
The prior shape and scale in the form of c(shape, scale) for an Inverse-Gamma(shape, scale) prior for tau2. Defaults to c(0.001, 0.001).  
}
  \item{prior.rho}{
A discrete prior is assigned for rho. The prior is the set of possible values in the interval [0,1). Defaults to (0, 0.01,...,0.98, 0.99). 
}
  \item{verbose}{
Logical, should the function update the user on its progress.  
}
}



\details{
For further details about how to apply the function see the examples below and in the vignette.
}


\value{
\item{formula }{A text variable stating a cluster model was fitted.}
\item{samples }{A list containing the MCMC samples from the model.}
\item{fitted.values }{A summary matrix of the posterior distributions of the fitted values for each area. The summaries include: Mean, Sd, Median, and credible interval.}
\item{residuals }{A summary matrix of the posterior distributions of the residuals for each area. The summaries include: Mean, Sd, Median, and credible interval.}
\item{W.summary }{The neighbourhood matrix W from the model.}
\item{modelfit }{Model fit criteria including the Deviance Information Criterion (DIC), the effective number of parameters in the model(p.d), DIC3 and the Marginal Predictive Likelihood (MPL). Additionally for this cluster model two criteria are given for choosing between the number of clusters, including  modified versions of the ratio of the within to between cluster sum of squares and the diagnostic proposed by Calinski, T. and J. Harabasz (1974).}
\item{summary.results }{A summary table of the parameters.}
\item{model }{A text string describing the model fit.}
\item{accept }{The acceptance probabilities for the parameters.}
}

\references{
Calinski, T. and J. Harabasz (1974). Reference -  A dendrite method for cluster analysis. Communications in Statistics 3, 1-27.
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
theta <- c(rep(1, 0.5*n), rep(-1, 0.5*n))
phi <- mvrnorm(n=1, mu=rep(0,n), Sigma=0.4 * exp(-0.1 * distance))
logit <- x1 + x2 + theta + phi
prob <- exp(logit) / (1 + exp(logit))
trials <- rep(50,n)
Y <- rbinom(n=n, size=trials, prob=prob)


#### Run the cluster model
formula <- Y ~ x1 + x2
\dontrun{model <- clusterCAR.re(Y=Y, q=2, family="binomial, trials=trials, 
W=W,burnin=5000, n.sample=10000)}
}