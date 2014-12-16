\name{summarise.samples}
\alias{summarise.samples}

\title{
Summarise a matrix of Markov chain Monte Carlo samples.
}

\description{
This function takes in a matrix of Markov chain Monte Carlo (McMC) samples from a CARBayes model object, such as from a set of parameters or fitted values, and calculates posterior quantiles and exceeedence probabilities. The latter are probabilities of the form P(quantity > c|data), where c is a threshold chosen by the user.
}

\usage{
summarise.samples(samples, columns=NULL, quantiles=0.5, exceedences=NULL)
}

\arguments{
  \item{samples}{
A matrix of McMC samples resulting from fitting one of the models. The object must be of class mcmc from the coda package.
}
  \item{columns}{
A vector of numbers stating which columns in the matrix of McMC samples you wish to summarise. Defaults to all columns.
}
 \item{quantiles}{
A vector of posterior quantiles required for each column of McMC samples.
}
     \item{exceedences}{
A vector of threshold levels, c, that you wish to calculate exceedence probabilities for, for each column of McMC samples.
}
}


\details{
For further details about how to apply the function see the example below.
}

\value{
\item{quantiles }{A 2 dimensional array containing the requied posterior quantiles. Each row relates to a parameter, and each column to a different requested quantile. If the argument quantiles is missing this object is NULL.
}
\item{exceedences }{A 2 dimensional array containing the requied exceedence probabilities. Each row relates to a parameter, and each column to a different requested exceedence probability. If the argument exceedences is missing this object is NULL.
}
}


\author{
Duncan Lee
}




\examples{
## See Vignette for an example of its use.
}

