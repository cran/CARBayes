\name{respdata}
\alias{respdata}
\docType{data}

\title{
Respiratory hosptial admissions data for the 271 Intermediate Geographies (IG) in the Greater Glasgow and Clyde health board.
}

\description{
Respiratory disease data for the 271 Intermediate Geographies (IG) in the Greater Glasgow and Clyde health board.These data  are used in the worked example in the main help page of CARBayes.
}
\usage{data(respdata)}
\format{
  A data frame with 271 observations on the following 3 variables.
  \describe{
    \item{\code{observed2010}}{The observed numbers of cases of hospitalisation due to respiratory disease in each IG in 2010.}
    \item{\code{expected2010}}{The expected numbers of cases of hospitalisation due to respiratory disease in each IG in 2010 computed using external standardisation from Scotland wide disease rates.}
    \item{\code{incomedep2010}}{The percentage of people in each IG who are defined to be income deprived.}
  }
}

\source{
These data were provided by the Scottish Government, and can also be downloaded from http://www.sns.gov.uk/. Further details about these data are also available at the same location.
}


\examples{
data(respdata)
}

\keyword{datasets}
