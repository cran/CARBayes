\name{binomial.independent}
\alias{binomial.independent}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Fit the independent random effects model to spatial binomial data
}

\description{
The function fits a binomial logistic random effects models to spatial data, where the random effects are modelled as independent. The model represents the logit transform of the mean function for the set of binomial responses  by a combination of covariates and a set of independent random effects. A set of offsets can also be included on the linear predictor scale. Inference is based on Markov Chain Monte Carlo (MCMC) simulation, using a combination of Gibbs sampling and Metropolis steps.
}


\usage{
binomial.independent(formula, data = NULL, beta = NULL, theta = NULL, sigma2 = NULL, 
trials, burnin = 0, n.sample = 1000, thin=1, blocksize.beta = 5, 
blocksize.theta = 10, prior.mean.beta = NULL, prior.var.beta = NULL, 
prior.max.sigma2 = NULL)
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
  \item{theta}{
A vector of starting values for the random effects. If this argument is not specified the function will randomly generate starting values.
}
  \item{sigma2}{
A starting value for the variance parameter of the random effects. If this argument is not specified the function will randomly generate a starting value.
}  
  \item{trials}{
A vector the same length as the response containing the total number of trials for each area.
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
  \item{blocksize.theta}{
The size of the blocks in which to update the random effects in the MCMC algorithm. Defaults to 10.
}
  \item{prior.mean.beta}{
A vector of prior means for the regression parameters beta (Gaussian priors are assumed). Defaults to a vector of zeros.
}
  \item{prior.var.beta}{
A vector of prior variances for the regression parameters beta (Gaussian priors are assumed). Defaults to a vector with values 1000.
}
  \item{prior.max.sigma2}{
The maximum allowable value for the random effects variance sigma2 (a Uniform(0,M) prior is assumed). Defaults to M=1000.
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
\item{W.summary }{NULL.}
\item{DIC }{The Deviance Information Criterion.}
\item{p.d }{The effective number of parameters in the model.}
\item{summary.results }{A summary table of the parameters.}
\item{model }{A text string describing the model fit.}
\item{accept }{The acceptance probabilities for the parameters.}

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
	
	
#### Generate the covariates and response data
x1 <- rnorm(n)
x2 <- rnorm(n)
theta <- rnorm(n, sd=0.05)
logit <- x1 + x2 + theta
prob <- exp(logit) / (1 + exp(logit))
trials <- rep(50,n)
Y <- rbinom(n=n, size=trials, prob=prob)


#### Run the independent model
#### Let the function randomly generate starting values for the parameters
#### Use the default priors specified by the function (for details see the help files)
formula <- Y ~ x1 + x2
\dontrun{model <- binomial.independent(formula=formula, trials=trials, burnin=5000, 
n.sample=10000)}
\dontshow{model <- binomial.independent(formula=formula, trials=trials, burnin=20, 
n.sample=50)}
}
