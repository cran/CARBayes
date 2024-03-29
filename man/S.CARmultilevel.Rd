\name{S.CARmultilevel}
\alias{S.CARmultilevel}
%- Also NEED an '\alias' for EACH other topic documented here.


\title{
    Fit a spatial generalised linear mixed model  to multi-level areal unit data, 
    where the spatial random effects have a Leroux conditional autoregressive 
    prior.
}

\description{
    Fit a spatial generalised linear mixed model to multi-level areal unit data, where
    the response variable can be binomial, Gaussian or Poisson. The data are 
    structured with individuals within areal units, and different numbers of individuals
    are allowed within each areal unit. The linear predictor is modelled by
    known covariates (either individual or areal level) and a vector of areal level
    random effects that are modelled by the conditional autoregressive prior 
    proposed by Leroux et al. (2000). Independent random effects can be obtained 
    by setting rho=0, while the intrinsic CAR model can be obtained by setting 
    rho=1. Inference is conducted in a Bayesian setting using Markov chain Monte 
    Carlo (MCMC) simulation. Missing (NA) values are allowed in the response, and 
    posterior predictive  distributions are created for the missing values using 
    data augmentation. These are  saved in the "samples" argument in the output 
    of the function and are denoted by "Y". For a full model specification see 
    the vignette accompanying this package. 
}


\usage{
    S.CARmultilevel(formula, family, data=NULL,  trials=NULL, W, ind.area, 
    burnin, n.sample, thin=1, n.chains=1, n.cores=1, prior.mean.beta=NULL, 
    prior.var.beta=NULL, prior.nu2=NULL, prior.tau2=NULL, rho=NULL, MALA=TRUE, 
    verbose=TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{formula}{
        A formula for the covariate part of the model using the syntax of the
        lm() function. Offsets can be included here using the offset() function.
        The response, offset and each covariate are vectors with length equal to
        the number of individuals. The response can contain missing (NA) values.
    }
    \item{family}{
        One of either "binomial", "gaussian", or "poisson", which respectively 
        specify a binomial likelihood model with a logistic link function, a 
        Gaussian likelihood model with an identity link function, or a Poisson 
        likelihood model with a log link function. 
    }
    \item{data}{
        An optional data.frame containing the  variables in the formula.
    }
    \item{trials}{
        A vector the same length as the response containing the total number of trials 
        for each individual. Only used if family="binomial". 
    }
    \item{W}{A non-negative K by K neighbourhood matrix (where K is the number of 
   spatial units). Typically a binary specification is used, where the jkth 
   element equals one if areas (j, k) are spatially close (e.g. share a common 
   border) and is zero otherwise. The matrix can be non-binary, but each row must 
   contain at least one non-zero entry.
}
        \item{ind.area}{
        A vector of integers the same length as the number of data points (individuals) 
        giving which spatial unit (nunmbered from 1 to K to align with the rows of
        the W matrix) each individual belongs to.
        }
  \item{burnin}{
The number of MCMC samples to discard as the burn-in period in each chain.
}
  \item{n.sample}{
The overall number of MCMC samples to generate in each chain.
}
  \item{thin}{
The level of thinning to apply to the MCMC samples in each chain to reduce their 
autocorrelation. Defaults to 1 (no thinning).
}
  \item{n.chains}{
The number of MCMC chains to run when fitting the model. Defaults to 1.
}
  \item{n.cores}{
The number of computer cores to run the MCMC chains on. Must be less than or 
equal to n.chains. Defaults to 1.
}
        \item{prior.mean.beta}{
        A vector of prior means for the regression parameters beta (Gaussian priors are 
        assumed). Defaults to a vector of zeros.
        }
        \item{prior.var.beta}{
        A vector of prior variances for the regression parameters beta (Gaussian priors 
        are assumed). Defaults to a vector with values 100,000.
        }  
        \item{prior.nu2}{
        The prior shape and scale in the form of c(shape, scale) for an Inverse-Gamma(shape, scale) 
        prior for nu2. Defaults to c(1, 0.01) and only used if family="Gaussian".   
        }
        \item{prior.tau2}{
        The prior shape and scale in the form of c(shape, scale) for an Inverse-Gamma(shape, scale) 
        prior for tau2. Defaults to c(1, 0.01).  
        }
        \item{rho}{
        The value in the interval [0, 1] that the spatial dependence parameter rho is 
        fixed at if it should not be estimated. If this arugment is NULL then rho is
        estimated in the model.
        }
        \item{MALA}{
        Logical, should the function use Metropolis adjusted Langevin algorithm 
        (MALA) updates (TRUE, default) or simple random walk updates (FALSE) for 
        the regression parameters. Not applicable if family="gaussian".   
}
        \item{verbose}{
        Logical, should the function update the user on its progress.  
        }
        }




\value{
    \item{summary.results }{A summary table of the parameters.}
    \item{samples }{A list containing the MCMC samples from the model.}
    \item{fitted.values }{The fitted values based on posterior means from the model.}
    \item{residuals }{A matrix with 2 columns where each column is a type of 
        residual and each row relates to an area. The types are "response" (raw), 
        and "pearson".}
\item{modelfit }{Model fit criteria including the Deviance Information Criterion 
(DIC) and its corresponding estimated effective number of parameters (p.d), the Log 
Marginal Predictive Likelihood (LMPL), the Watanabe-Akaike Information Criterion 
(WAIC) and its corresponding estimated number of effective parameters (p.w), and
the loglikelihood.}
    \item{localised.structure }{NULL, for compatability with other models.}
\item{formula }{The formula (as a text string) for the response, covariate and 
offset parts of the model}
    \item{model }{A text string describing the model fit.}
    \item{mcmc.info }{A vector giving details of the numbers of MCMC samples generated.}
    \item{X }{The design matrix of covariates.}
}



\author{
    Duncan Lee
}




\examples{
#################################################
#### Run the model on simulated data on a lattice
#################################################
    
#### Set up a square lattice region
x.easting <- 1:10
x.northing <- 1:10
Grid <- expand.grid(x.easting, x.northing)
K <- nrow(Grid)
    
#### set up distance and neighbourhood (W, based on sharing a common border) matrices
distance <- as.matrix(dist(Grid))
W <-array(0, c(K,K))
W[distance==1] <-1 	

#### Generate the number of individuals per area and which individuals to which areas
n <- sample(5:30, K, replace=TRUE)
n.total <- sum(n)
ind.area.temp <- rep(1:K, n)
ind.area <- sample(ind.area.temp, n.total, replace=FALSE)

#### Generate the covariates and response data
x1 <- rnorm(n.total)
x2 <- rnorm(n.total)
phi <- mvrnorm(n=1, mu=rep(0,K), Sigma=0.4 * exp(-0.1 * distance))
phi.extend <- phi[ind.area]
logit <- x1 + x2 + phi.extend
prob <- exp(logit) / (1 + exp(logit))
trials <- rep(50,n.total)
Y <- rbinom(n=n.total, size=trials, prob=prob)
    
    
#### Run the model
formula <- Y ~ x1 + x2
\dontrun{model <- S.CARmultilevel(formula=formula, family="binomial", ind.area=ind.area,
                trials=trials, W=W, burnin=20000, n.sample=100000)}
                
#### Toy example for checking
model <- S.CARmultilevel(formula=formula, family="binomial", ind.area=ind.area,
                trials=trials, W=W, burnin=10, n.sample=50)
}
