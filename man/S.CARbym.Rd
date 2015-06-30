\name{S.CARbym}
\alias{S.CARbym}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Fit a spatial generalised linear mixed model to data, where the random effects 
have a BYM conditional autoregressive prior.
}

\description{
Fit a spatial generalised linear mixed model to areal unit data, where the response
variable can be binomial, or Poisson. Note, a Gaussian likelihood is not allowed 
because of a lack of identifiability among the parameters. The linear predictor is 
modelled byknown covariates and a vector of random effects. The latter are modelled 
by the BYM conditional autoregressive prior proposed by Besag et al. (1991), and 
further details are given in the vignette accompanying this package. Inference is 
conducted in a Bayesian setting using Markov chain Monte Carlo (McMC) simulation. 
}


\usage{
S.CARbym(formula, family, data=NULL, trials=NULL, W, burnin, n.sample, 
thin=1, prior.mean.beta=NULL, prior.var.beta=NULL,  prior.tau2=NULL, 
prior.sigma2=NULL, verbose=TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
A formula for the covariate part of the model using the syntax of the
lm() function. Offsets can be included here using the offset() function.
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
   \item{W}{
A K by K neighbourhood matrix (where K is the number of spatial units). Typically 
a binary specification is used, where the jkth element equals one if areas (j, k) 
are spatially close (e.g. share a common border) and is zero otherwise. The matrix 
can be non-binary, but each row must contain at least one non-zero entry.
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
     \item{prior.tau2}{
The prior shape and scale in the form of c(shape, scale) for an Inverse-Gamma(shape, scale) 
prior for tau2. Defaults to c(0.001, 0.001).  
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




\references{
Besag, J., J. York, and A. Mollie (1991). Bayesian image restoration with two 
applications in spatial statistics. Annals of the Institute of Statistics and 
Mathematics 43, 1-59.
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

#### set up distance and neighbourhood (W, based on sharing a common border) matrices
distance <-array(0, c(K,K))
W <-array(0, c(K,K))
  for(i in 1:K)
	{
		for(j in 1:K)
		{
		temp <- (Grid[i,1] - Grid[j,1])^2 + (Grid[i,2] - Grid[j,2])^2
		distance[i,j] <- sqrt(temp)
			if(temp==1)  W[i,j] <- 1 
		}	
	}
	
	
#### Generate the covariates and response data
x1 <- rnorm(K)
x2 <- rnorm(K)
theta <- rnorm(K, sd=0.05)
phi <- mvrnorm(n=1, mu=rep(0,K), Sigma=0.4 * exp(-0.1 * distance))
logit <- x1 + x2 + theta + phi
prob <- exp(logit) / (1 + exp(logit))
trials <- rep(50,K)
Y <- rbinom(n=K, size=trials, prob=prob)


#### Run the BYM model
formula <- Y ~ x1 + x2
\dontrun{model <- S.CARbym(formula=formula, family="binomial", trials=trials,
W=W, burnin=20000, n.sample=100000)}
}