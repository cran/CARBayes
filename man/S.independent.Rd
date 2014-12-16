\name{S.independent}
\alias{S.independent}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Fit a generalised linear model with a set of independent random effect to spatial data.}

\description{
The function fits a Bayesian hierarchical model with independent random effects to the data, whre the data likelihood can be binomial, Gaussian or Poisson. The model represents the linear predictor by a combination of covariates and a set of independent random effects. A set of offsets can also be included on the linear predictor scale. Inference is based on Markov chain Monte Carlo (McMC) simulation, using a combination of Gibbs sampling and Metropolis steps.
}


\usage{
S.independent(formula, family, data=NULL,  trials=NULL, burnin=0, 
n.sample=1000, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, 
prior.nu2=NULL, prior.sigma2=NULL, verbose=TRUE)
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
  \item{prior.mean.beta}{
A vector of prior means for the regression parameters beta (Gaussian priors are assumed). Defaults to a vector of zeros.
}
  \item{prior.var.beta}{
A vector of prior variances for the regression parameters beta (Gaussian priors are assumed). Defaults to a vector with values 1000.
}
     \item{prior.nu2}{
Only used for the Gaussian model. The prior shape and scale in the form of c(shape, scale) for an Inverse-Gamma(shape, scale) prior for nu2. Defaults to c(0.001, 0.001).  
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
\item{summary.results }{A summary table of the parameters.}
\item{samples }{A list containing the McMC samples from the model.}
\item{fitted.values }{A vector of fitted values for each area.}
\item{residuals }{A vector of residuals for each area.}
\item{modelfit }{Model fit criteria including the Deviance Information Criterion (DIC), the effective number of parameters in the model(p.d), and the Log Marginal Predictive Likelihood (LMPL).}
\item{accept }{The acceptance probabilities for the parameters.}
\item{localised.structure }{NULL, for compatability with the other models.}
\item{formula }{The formula for the covariate and offset part of the model.}
\item{model }{A text string describing the model fit.}
\item{X }{The design matrix of covariates.}
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
formula <- Y ~ x1 + x2
\dontrun{model <- S.independent(formula=formula, family="binomial", trials=trials, 
burnin=5000, n.sample=10000)}
}
