\name{poisson.dissimilarityCAR}
\alias{poisson.dissimilarityCAR}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Fit a localised conditional autoregressive (CAR) model based on dissimilarity metrics with binary neighbourhood relations to spatial Poisson data
}

\description{
The function fits a Poisson log-normal random effects models to spatial count data, where the random effects are modelled by the localised conditional autoregressive (CAR) model proposed by Lee and Mitchell (2012). The random effects in neighbouring areas (e.g. those that share a common border) are modelled as correlated or conditionally independent, depending on whether the populations living in the two areas are similar (correlated random effects) or very different (conditionally independent). The model represents the natural log of the mean function for the set of Poisson responses  by a combination of covariates and a set of random effects. Inference is based on Markov Chain Monte Carlo (MCMC) simulation, using a combination of Gibbs sampling and Metropolis steps.
}

\usage{
poisson.dissimilarityCAR(formula, data=NULL, beta=NULL, phi=NULL, tau2=NULL, 
rho=0.99, alpha=NULL, W, Z, burnin=0, n.sample=1000, thin=1, blocksize.beta=5, 
prior.mean.beta=NULL, prior.var.beta=NULL, prior.tau2=NULL)
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
The fixed value for the global spatial correlation parameter rho. The default is 0.99 which should be retained.
}
  \item{alpha}{
A vector of starting values for the regression parameters that control the identification of boundaries in the random effects surface. If this argument is not specified the function will randomly generate starting values.
}
  \item{W}{
A binary n by n neighbourhood matrix (where n is the number of spatial units). The jkth element equals one if areas (j, k) are spatially close (e.g. share a common border) and is zero otherwise. 
}
  \item{Z}{
A list, where each element is an n by n matrix of non-negative dissimilarity metrics. For details of their construction see the main help page for CARBayes
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
  \item{prior.tau2}{
The prior shape and scale in the form of c(shape, scale) for an Inverse-Gamma(shape, scale) prior for tau2. Defaults to c(0.001, 0.001).
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
\item{W.summary }{A list containing two matrices: W.posterior contains posterior medians for each element wkj of the n by n neighbourhood matrix W; W.border.prob contains posterior probabilities that each wkj element of the n by n neighbourhood matrix W equals zero. This corresponds to the posterior probability of a boundary in the random effects surface. In both cases elements which correspond to two non-neighbouring areas have NA values. }
\item{DIC }{The Deviance Information Criterion.}
\item{p.d }{The effective number of parameters in the model.}
\item{MPL }{The Marginal Predictive Likelihood of the model.}
\item{summary.results }{A summary table of the parameters.}
\item{model }{A text string describing the model fit.}
\item{accept }{The acceptance probabilities for the parameters.}
}


\references{
Lee, D. and R. Mitchell (2012). Boundary detection in disease mapping studies. Biostatistics, 13, 415-426.
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
E <- rep(40,n)
phi <- mvrnorm(n=1, mu=groups, Sigma=0.2 * exp(-0.1 * distance))
risk <- exp(-1 + phi)
fitted <- E * risk
Y <- rpois(n=n, lambda=fitted)


#### Generate a dissimilarity metric
dissimilarity <- cbind(groups) + rnorm(n, sd=0.2)
dissimilarity.matrix <- as.matrix(dist(cbind(dissimilarity, dissimilarity), 
method="manhattan", diag=TRUE, upper=TRUE)) * W/2

Z <- list(dissimilarity.matrix=dissimilarity.matrix)

#### Run the localised smoothing model
#### Let the function randomly generate starting values for the parameters
#### Use the default priors specified by the function (for details see the help files)
formula <- Y ~ offset(log(E))
\dontrun{model <- poisson.dissimilarityCAR(formula=formula, rho=0.99,
W=W, Z=Z, burnin=5000, n.sample=10000)}
\dontshow{model <- poisson.dissimilarityCAR(formula=formula, rho=0.99,
W=W, Z=Z, burnin=20, n.sample=50)}
}