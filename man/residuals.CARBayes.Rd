\name{residuals.CARBayes}
\alias{residuals.CARBayes}


\title{
    Extract the residuals from a model.
}
\description{
    This function takes a CARBayes object and returns a set of residuals. 
    The allowable types of residual are "response" (raw), "pearson", and "deviance". 
}

\usage{
    \method{residuals}{CARBayes}(object, type, \dots)
}



\arguments{
    \item{object}{
        A CARBayes fitted model object.
    }
    \item{type}{
        A text string and one of c("response", "pearson", "deviance").
    }
    \item{\dots}{
        Ignored.
    }
}
\author{
    Duncan Lee
}
