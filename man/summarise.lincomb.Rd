\name{summarise.lincomb}
\alias{summarise.lincomb}


\title{
Compute the posterior distribution and quantiles of a linear combination of the covariate component of the linear predictor.  
}

\description{
This function takes in a CARBayes model object and computes the posterior distribution and quantiles of a linear combination of the covariate component of the linear predictor. For example, if in the model the user specified a quadratic effect of a covariate on the response, then this function allows you to compute the posterior distribution of the entire quadratic relationship.
}
\usage{
summarise.lincomb(model, columns=NULL, quantiles=0.5, distribution=FALSE)
}


\arguments{
  \item{model}{
A CARBayes model object from fitting one of the models in this package.
}
  \item{columns}{
A vector of numbers stating which columns in the design matrix of covariates you wish to compute the posterior distribution for.
}
 \item{quantiles}{
A vector of posterior quantiles required for each column of McMC samples.
}
     \item{distribution}{
A logical value stating whether the entire posterior distribution should be returned or just the specified quantiles.
}
}

\details{
For further details about how to apply the function see the example below.
}

\value{
\item{quantiles }{A 2 dimensional array containing the requied posterior quantiles. Each row relates to a data value, and each column to a different requested quantile. 
}
\item{posterior }{A 2 dimensional array containing the requied posterior distribution. Each column relates to a different data value. Is NULL if distribution=FALSE in the function arguments.
}
}



\author{
Duncan Lee
}


\examples{
## See Vignette for an example of its use.
}

