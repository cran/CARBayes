\name{S.CARcluster}
\alias{S.CARcluster}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Fit a spatial generalised linear mixed model to data, where the random effects 
have a localised clustered conditional autoregressive prior.
}

\description{
Fit a spatial generalised linear mixed model to areal unit data, where the response
variable can be binomial, Gaussian or Poisson. The linear predictor is modelled by
known covariates and a vector of random effects. The latter are modelled by the 
localised clustered conditional autoregressive prior proposed by Lee and Sarran
(2014), and further details are given in the vignette accompanying this package. 
Inference is conducted in a Bayesian setting using Markov chain Monte Carlo (McMC)
simulation. 
}

\usage{
S.CARcluster(formula, exposure=NULL, family, data=NULL, G, trials=NULL, W,
burnin, n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL,
prior.mean.alpha=NULL, prior.var.alpha=NULL, prior.nu2=NULL, prior.tau2=NULL,
prior.delta=NULL, verbose=TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
A formula for the covariate part of the model using the syntax of the
lm() function. Offsets can be included here using the offset() function.
}
  \item{exposure}{
This optional argument allows a single covariate to have variation in its values 
within an areal unit if the data are Poisson counts. The argument is a list object,
and each element of the list corresponds to a single areal unit. Each element in
the list is a 2 column matrix, where the first column is the vector of values of 
the covariate while the second column is a vector of weights for each value. 
}
  \item{family}{
One of either `binomial', `gaussian' or `poisson', which respectively specify a 
binomial likelihood model with a logistic link function, a Gaussian likelihood 
model with an identity link function, or a Poisson likelihood model with a 
log link function. 
}
  \item{data}{
An optional data.frame containing the  variables in the formula.
}

\item{G}{
The maximum number of distinct intercept terms (clusters) to allow in the model.
}
  \item{trials}{
A vector the same length as the response containing the total number of trials 
for each area. Only used if family=`binomial'. 
}
   \item{W}{
A binary n by n neighbourhood matrix (where n is the number of spatial units). 
The jkth element equals one if areas (j, k) are spatially close (e.g. share a common
border) and is zero otherwise. 
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
  \item{prior.mean.alpha}{
The prior mean for the regression parameter alpha (Gaussian priors are assumed) 
multiplying the covariate with within area variation in its values. Defaults to zero.
}
  \item{prior.var.alpha}{
The prior variance for the regression parameter alpha (Gaussian priors are assumed) 
multiplying the covariate with within area variation in its values. Defaults to 1000.
}  
  \item{prior.nu2}{
The prior shape and scale in the form of c(shape, scale) for an Inverse-Gamma(shape, scale) 
prior for nu2. Defaults to c(0.001, 0.001) and only used if family=`Gaussian'.   
}
     \item{prior.tau2}{
The prior shape and scale in the form of c(shape, scale) for an Inverse-Gamma(shape, scale) 
prior for tau2. Defaults to c(0.001, 0.001).  
}
     \item{prior.delta}{
The prior maximum for the cluster smoothing parameter delta. Defaults to 10.  
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
\item{localised.structure }{A vector giving the posterior median of which cluster 
(group) each area is in.}
\item{formula }{The formula for the covariate and offset part of the model.}
\item{model }{A text string describing the model fit.}
\item{X }{The design matrix of covariates.}
}

\references{
Lee, D and Sarran, C  (2014). Controlling for unmeasured confounding and spatial 
misalignment in long-term air pollution and health studies. arXiv:1412.4479.
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
formula <- Y ~ 1
\dontrun{model <- S.CARcluster(formula=formula, family="binomial", G=2,  
W=W,burnin=20000, n.sample=100000)}
}