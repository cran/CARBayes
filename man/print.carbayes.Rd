\name{print.carbayes}
\alias{print.carbayes}


\title{
Print a summary of a fitted carbayes model to the screen.
}
\description{
This function takes a carbayes object and returns a summary of the fitted model.
The summary includes posterior medians and 95 percent credible intervals, the
number of samples, the acceptance rate, the effective number of independent samples 
and the Geweke convergence diagnostic in the form of a Z-score.
}

\usage{
\method{print}{carbayes}(x, \dots)
}



\arguments{
  \item{x}{
A `carbayes' fitted model object
}
  \item{\dots}{
Ignored
}
}
\author{
Duncan Lee
}
