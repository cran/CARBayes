\name{S.independent}
\alias{S.independent}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Fit a spatial generalised linear mixed model to data, where the random effects 
are independent.
}

\description{
Fit a spatial generalised linear mixed model to areal unit data, where the response
variable can be binomial or Poisson. Note, a Gaussian likelihood is not allowed 
because of a lack of identifiability among the parameters. The linear predictor is 
modelled by known covariates and a vector of random effects. The latter are modelled as 
independent and identically distributed, and further  details are given in the 
vignette accompanying this package. Inference is conducted in a Bayesian setting
using Markov chain Monte Carlo (McMC) simulation. Missing (NA) values are allowed in 
the response. Missing (NA) values are allowed in the response, and posterior predictive 
distributions are created for the missing values for predictive purposes. These are s
aved in the `samples' argument in the output of the function and are denoted by `Y'.
}


\usage{
S.independent(formula, family, data=NULL,  trials=NULL, burnin, n.sample, thin=1, 
prior.mean.beta=NULL, prior.var.beta=NULL, prior.sigma2=NULL, verbose=TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
A formula for the covariate part of the model using the syntax of the
lm() function. Offsets can be included here using the offset() function. The 
response can contain missing (NA) values.
}
  \item{family}{
One of either `binomial' or `poisson', which respectively specify a 
binomial likelihood model with a logistic link function, or a Poisson likelihood 
model with a log link function. 
}
     \item{data}{
An optional data.frame containing the  variables in the formula.
}
  \item{trials}{
A vector the same length as the response containing the total number of trials 
for each area. Only used if family=`binomial'. 
}
  \item{burnin}{
The number of McMC samples to discard as the burnin period.
}
  \item{n.sample}{
The number of McMC samples to generate.
}
  \item{thin}{
The level of thinning to apply to the McMC samples to reduce their temporal 
autocorrelation. Defaults to 1.
}
  \item{prior.mean.beta}{
A vector of prior means for the regression parameters beta (Gaussian priors are 
assumed). Defaults to a vector of zeros.
}
  \item{prior.var.beta}{
A vector of prior variances for the regression parameters beta (Gaussian priors 
are assumed). Defaults to a vector with values 1000.
}  
  \item{prior.sigma2}{
The prior shape and scale in the form of c(shape, scale) for an Inverse-Gamma(shape, scale)
prior for sigma2. Defaults to c(0.001, 0.001).  
} 
     \item{verbose}{
Logical, should the function update the user on its progress.  
}
}



\value{
\item{summary.results }{A summary table of the parameters.}
\item{samples }{A list containing the McMC samples from the model.}
\item{fitted.values }{A vector of fitted values for each area.}
\item{residuals }{A vector of residuals for each area.}
\item{modelfit }{Model fit criteria including the Deviance Information Criterion 
(DIC), the effective number of parameters in the model (p.d), and the Log Marginal 
Predictive Likelihood (LMPL).}
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
K <- nrow(Grid)
	
	
#### Generate the covariates and response data
x1 <- rnorm(K)
x2 <- rnorm(K)
theta <- rnorm(K, sd=0.05)
logit <- x1 + x2 + theta
prob <- exp(logit) / (1 + exp(logit))
trials <- rep(50,K)
Y <- rbinom(n=K, size=trials, prob=prob)


#### Run the independent model
formula <- Y ~ x1 + x2
\dontrun{model <- S.independent(formula=formula, family="binomial", trials=trials, 
burnin=20000, n.sample=100000)}
}
