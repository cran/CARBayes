\name{S.RAB}
\alias{S.RAB}
%- Also NEED an '\alias' for EACH other topic documented here.


\title{
Fit a spatial generalised linear model with anisotropic basis functions to data
for computationally efficient localised spatial smoothing, where the parameters 
are estimated by penalised maximum likelihood estimation with a ridge regression 
penalty.
}

\description{
Fit a spatial generalised linear model to areal unit data, where the response
variable can be binomial, Gaussian or Poisson. The linear predictor is modelled 
by known covariates and a set of K anisotropic spatial basis functions. The basis
functions are constructed from the set of geodesic distances between all pairs of 
areal units and a vector of ancillary data V, and the latter should have a similar
spatial pattern to the residual (after covariate adjustment) spatial structure in 
the data on the linear predictor scale. Parameter estimtion is carried out via 
penalised maximum likelihood methods, and the basis function coefficients are
constrained by a ridge regression penalty to prevent overfitting. The glmnet 
package is used for parameter estimation. Missing (NA) values are allowed in the 
response, and predictions are made for these values. This model implements localised
spatial smoothing and allows for boundaries in the data surface using a 
computationally efficient approach. 
}


\usage{
S.RAB(formula, family, data=NULL, trials=NULL, W, V, nlambda=100, verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
A formula for the covariate part of the model using the syntax of the
lm() function. Offsets can be included here using the offset() function. The 
response, offset and each covariate is a vector of length K*1. The response can 
contain missing (NA) values.
}
  \item{family}{
One of either "binomial", "gaussian", or "poisson", which respectively 
specify a binomial likelihood model with a logistic link function, 
a Gaussian likelihood model with an identity link function, or a Poisson likelihood 
model with a log link function. 
}
  \item{data}{
An optional data.frame containing the  variables in the formula.
}
 \item{trials}{
A vector the same length as the response containing the total number of trials 
for each area. Only used if family="binomial". 
}
   \item{W}{A non-negative K by K neighbourhood matrix (where K is the number of 
   spatial units). Typically a binary specification is used, where the jkth 
   element equals one if areas (j, k) are spatially close (e.g. share a common 
   border) and is zero otherwise. The matrix can be non-binary, but each row must 
   contain at least one non-zero entry.
}
  \item{V}{
A vector of ancillary data of length K, which should have a similar spatial 
pattern to the residual (after covariate adjustment) spatial structure in 
the data on the linear predictor scale.
}
  \item{nlambda}{
The number of possible values to use for the penalty parameter lambda in the glmnet()
estimation function. Defaults to 100.
}
  \item{verbose}{
Logical, should the function update the user on its progress.  
}
}



\value{
\item{beta.summary }{A summary table of the estimated regression parameters.}
\item{beta.covariance }{The estimated covariance matrix for the regression parameters.}
\item{sigma2.hat }{The estimated error variance in the Gaussian data likelihood 
model. If a Gaussian model is not specified it is NA.}
\item{lambda.hat }{The estimated ridge regression penalty parameter.}
\item{I }{The level of residual spatial autocorrelation as measured by
Moran's I statistic.}
\item{fitted.values }{The fitted values from the model.}
    \item{residuals }{A matrix with 2 columns where each column is a type of 
        residual and each row relates to an area. The types are "response" (raw), 
        and "pearson".}
\item{formula }{The formula (as a text string) for the response, covariate and 
offset parts of the model.}
\item{model.string }{A text string describing the model fit.}
\item{X }{The design matrix of covariates and spatial basis functions.}
\item{model }{The fitted model object from the glmnet() function.}
}

\references{
Leroux B, Lei X, Breslow N (2000). "Estimation of Disease Rates in SmallAreas: A 
New Mixed Model for Spatial Dependence." In M Halloran, D Berry (eds.), 
\emph{Statistical Models in Epidemiology, the Environment and Clinical Trials},
pp. 179-191. Springer-Verlag, New York.
}

\author{
Duncan Lee
}




\examples{
#################################################
#### Run the model on simulated data on a lattice
#################################################
#### Load other libraries required
library(MASS)

#### Set up a square lattice region
x.easting <- 1:10
x.northing <- 1:10
Grid <- expand.grid(x.easting, x.northing)
K <- nrow(Grid)

#### set up distance and neighbourhood (W, based on sharing a common border) matrices
distance <- as.matrix(dist(Grid))
W <-array(0, c(K,K))
W[distance==1] <-1 	
	
#### Generate the spatial covariance structure
dists <- as.numeric(distance[upper.tri(distance)])
dists.quant <- quantile(dists, 0.05)
rho <- log(0.75) / -dists.quant
Sigma <- exp(-rho * distance)

#### Generate the boundaries
groups <-rep(0, K) 
groups[Grid$Var2>5] <- 1

#### Generate the covariates and response data
x1 <- rnorm(K)
x2 <- rnorm(K)
phi <- mvrnorm(n=1, mu=rep(0,K), Sigma=0.1 * exp(-rho * distance))
logit <- x1 +  x2 + phi + 0.4 * groups
prob <- exp(logit) / (1 + exp(logit))
trials <- rep(100,K)
Y <- rbinom(n=K, size=trials, prob=prob)

#### Generate the ancillary data
V <- rnorm(n=K, mean=phi + 0.4*groups , sd=rep(0.05,K))

#### Run the RAB model
mod <- S.RAB(formula=Y~x1+x2, family="binomial", data=NULL, trials=trials, W=W, 
        V=V, nlambda=50, verbose=TRUE)
}
