\name{MVS.CARleroux}
\alias{MVS.CARleroux}
%- Also NEED an '\alias' for EACH other topic documented here.


\title{
    Fit a multivariate spatial generalised linear mixed model to data, where the 
    random effects are modelled by a multivariate conditional autoregressive model.
}

\description{
    Fit a multivariate spatial generalised linear mixed model to areal unit data, 
    where the response variable can be binomial or Poisson. The linear predictor is 
    modelled by known covariates and a vector of random effects. The latter account
    for both spatial and between variable correlation, via a Kronecker product
    formulation. Spatial correlation is captured by the conditional autoregressive 
    (CAR) prior proposed by Leroux et al. (1999), and between variable correlation is 
    captured by a between variable covariance matrix with no fixed structure. This is
    a type of multivariate conditional autoregressive (MCAR).Further details are given 
    in the vignette accompanying this package. Independent (over space) random effects 
    can be obtained by setting (fix.rho=TRUE, rho=0), in which case the neighbourhood 
    matrix W is not part of the model. In this case enter a fake W matrix that is a 
    K by K matrix of zeros, where K is the number of spatial units. Similarly, the 
    intrinsic MCAR model can be obtained by setting (fix.rho=TRUE, rho=1). Inference 
    is conducted in a Bayesian setting using Markov chain Monte Carlo (McMC) simulation. 
    Missing (NA) values are allowed in the response, and posterior predictive 
    distributions are created for the missing values for predictive purposes. These are 
    saved in the`samples' argument in the output of the function and are denoted by `Y'. 
}


\usage{
MVS.CARleroux(formula, family, data=NULL,  trials=NULL, W, burnin, 
n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL,
prior.Sigma.df=NULL, prior.Sigma.scale=NULL, fix.rho=FALSE, rho=NULL, verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{formula}{
A formula for the covariate part of the model using the syntax of the
lm() function. Offsets can be included here using the offset() function. The 
response and each covariate should be vectors of length (KJ)*1, where
K is the number of spatial units and J is the number of different variables. Each
vector should be ordered so that the first J data points are the set of all J
variables for spatial location 1, the next J are the set of variables for spatial 
location 2 and so on. Missing (NA) values are allowed in the response.
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
    The number of McMC samples to discard as the burn-in period.
}
\item{n.sample}{
    The number of McMC samples to generate.
}
\item{thin}{
    The level of thinning to apply to the McMC samples to reduce their temporal 
    autocorrelation. Defaults to 1 (no thinning).
}
\item{prior.mean.beta}{
    A vector of prior means for the regression parameters beta (Gaussian priors are 
    assumed). Defaults to a vector of zeros.
}
\item{prior.var.beta}{
    A vector of prior variances for the regression parameters beta (Gaussian priors 
     are assumed). Defaults to a vector with values 1000.
}  
\item{prior.Sigma.df}{
    The prior degrees of freedom for the Inverse-Wishart prior for Sigma. 
    Defaults to J+1.   
}
\item{prior.Sigma.scale}{
    The prior J times J scale matrix for the Inverse-Wishart prior for Sigma. 
    Defaults to the identity matrix. 
}
\item{fix.rho}{
Logical, should the spatial dependence parameter rho be fixed rather than being estimated
in the model.
}
\item{rho}{
The value that the spatial dependence parameter is set to if fix.rho=TRUE. Must be in the
interval [0, 1], where 0 corresponds to independence while 1 corresponds to fitting the
intrinsic CAR prior. Note, if rho=0 then the W matrix is not included in the model
so enter a K by K matrix of zeros instead.
}
\item{verbose}{
Logical, should the function update the user on its progress.  
}
}




\value{
\item{summary.results }{A summary table of the parameters.}
\item{samples }{A list containing the McMC samples from the model.}
\item{fitted.values }{A vector of fitted values for each area and response variable
in the same order as the data vectors.}
\item{residuals }{A vector of fitted values for each area and response variable
in the same order as the data vectors.}
\item{modelfit }{Model fit criteria including the Deviance Information Criterion 
(DIC) and its corresponding estimated effective number of parameters (p.d), the Log 
Marginal Predictive Likelihood (LMPL), and the Watanabe-Akaike Information Criterion 
(WAIC) and its corresponding estimated number of effective parameters (p.w).}
\item{accept }{The acceptance probabilities for the parameters.}
\item{localised.structure }{NULL, for compatability with the other models.}
\item{formula }{The formula for the covariate and offset part of the model.}
\item{model }{A text string describing the model fit.}
\item{X }{The design matrix of covariates.}
}

\references{
Gelfand, A and Vounatsou, P (2003). Proper multivariate conditional autoregressive
models for spatial data analysis, Biostatistics, 4, 11-25.

Kavanagh, L., D. Lee, and G. Pryce (2016). Is Poverty Decentralising? Quantifying 
Uncertainty in the Decentralisation of Urban Poverty, Annals of the American 
Association of Geographers, to appear.

Leroux, B., X. Lei, and N. Breslow (1999). Estimation of disease rates in small 
areas: A new mixed model for spatial dependence, Chapter Statistical Models in 
Epidemiology, the Environment and Clinical Trials, Halloran, M and Berry, D (eds),
pp. 135-178. Springer-Verlag, New York.
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
K <- nrow(W)

#### Generate the correlation structures
Q.W <- 0.99 * (diag(apply(W, 2, sum)) - W) + 0.01 * diag(rep(1,K))
Q.W.inv <- solve(Q.W)

Sigma <- matrix(c(1,0.5,0, 0.5,1,0.3, 0, 0.3, 1), nrow=3)
Sigma.inv <- solve(Sigma)
J <- nrow(Sigma)

precision.phi <- kronecker(Q.W, Sigma.inv)
var.phi <- solve(precision.phi)

#### Generate data
N.all <- K * J
x1 <- rnorm(N.all)
x2 <- rnorm(N.all)
phi <- mvrnorm(n=1, mu=rep(0,N.all), Sigma=var.phi)

lp <- 0.1 * x1 - 0.1 * x2 + phi
prob <- exp(lp) / (1 + exp(lp))
trials <- rep(100,N.all)
Y <- rbinom(n=N.all, size=trials, prob=prob)

#### Run the Leroux model
formula <- Y ~ x1 + x2
\dontrun{model <- MVS.CARleroux(formula=formula, family="binomial", 
trials=trials, W=W, burnin=20000, n.sample=100000)}
}
