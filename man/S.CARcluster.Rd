\name{S.CARcluster}
\alias{S.CARcluster}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Fit a generalised linear model with a set of spatially autocorrelated random effects following a conditional autoregressive (CAR) prior and a piecewise constant jump component proposed by Lee and Sarran (2014) to spatial data.
}

\description{
The function fits a Bayesian hierarchical model with spatially correlated random effects and a cluster component to the data, whre the data likelihood can be binomial, Gaussian or Poisson. The random effects are modelled by the conditional autoregressive (CAR) model proposed by Leroux et. al. (1999). The model represents the linear predictor of the data  by a cluster component and a set of random effects. The latter are spatially correlated and come from the Leroux CAR model. Inference is based on Markov chain Monte Carlo (McMC) simulation, using a combination of Gibbs sampling and Metropolis steps.
}

\usage{
S.CARcluster(formula, exposure=NULL, family, data=NULL, G, trials=NULL, W,
burnin=0, n.sample=1000, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL,
prior.mean.alpha=NULL, prior.var.alpha=NULL, prior.nu2=NULL, prior.tau2=NULL,
prior.delta=NULL, verbose=TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
A formula for the covariate part of the model, using the same notation as for the lm() function. The offsets should also be included here using the offset() function.
}
  \item{exposure}{
This argument allows a single covariate to have variation in its values within an areal unit if the data are Poisson counts. This argument is not required and is only applicable if family="poisson". If such a covariate exists the argument is a list object, in which each element in the list corresponds to each areal unit. Each element in the list is a matrix where the first column is the vector of values of the covariate and the second column is the vector of weights for each value. See Lee and Sarran (2014) for more details.
}
  \item{family}{
One of either 'binomial', 'gaussian' or 'poisson', which respectively specify a binomial likelihood model with a logistic link function, a Gaussian likelihood model with an identity link function, or a Poisson likelihood model with a log link function. 
}
  \item{data}{
A data.frame containing the  variables in the formula.
}

\item{G}{
The number of clusters to fit to the data. Must be an integer with minimum value G=2.
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
  \item{prior.mean.beta}{
A vector of prior means for the regression parameters beta (Gaussian priors are assumed). Defaults to a vector of zeros.
}
  \item{prior.var.beta}{
A vector of prior variances for the regression parameters beta (Gaussian priors are assumed). Defaults to a vector with values 1000.
}  
  \item{prior.mean.alpha}{
The prior mean for the regression parameter alpha (Gaussian priors are assumed) multiplying the covariate with within area variation in values. Defaults to zero.
}
  \item{prior.var.alpha}{
The prior variance for the regression parameter alpha (Gaussian priors are assumed) multiplying the covariate with within area variation in values. Defaults to 1000.
}  
  \item{prior.nu2}{
The prior shape and scale in the form of c(shape, scale) for an Inverse-Gamma(shape, scale) prior for nu2. Only for the Gaussian model. Defaults to c(0.001, 0.001).  
}
     \item{prior.tau2}{
The prior shape and scale in the form of c(shape, scale) for an Inverse-Gamma(shape, scale) prior for tau2. Defaults to c(0.001, 0.001).  
}
     \item{prior.delta}{
The prior maximum for the cluster smoothing parameter delta. Defaults to 10.  
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
\item{localised.structure }{A vector giving the posterior median of which cluster group each area is in.}
\item{formula }{The formula for the covariate and offset part of the model.}
\item{model }{A text string describing the model fit.}
\item{X }{The design matrix of covariates.}
}

\references{
Leroux, B., X. Lei, and N. Breslow (1999). Estimation of disease rates in small areas: A new mixed model for spatial dependence, Chapter Statistical Models in Epidemiology, the Environment and Clinical Trials, Halloran, M and Berry, D (eds), pp. 135-178. Springer-Verlag, New York.

Lee, D and Sarran, C  (2014). Controlling for unmeasured confounding and spatial misalignment in long-term air pollution and health studies. arXiv:1412.4479.
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
\dontrun{model <- S.CARcluster(formula=formula, family="poisson", G=2,  
W=W,burnin=5000, n.sample=10000)}
}