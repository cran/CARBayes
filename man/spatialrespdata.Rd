\name{spatialrespdata}
\alias{spatialrespdata}
\docType{data}

\title{
A SpatialPolygonsDataFrame object (from the sp package) containing respiratory hosptial admissions data for the 134 Intermediate Geographies (IG) to the north of the river Clyde in the Greater Glasgow and Clyde health board.
}

\description{
A SpatialPolygonsDataFrame object containing data on respiratory disease data for the 134 Intermediate Geographies (IG) to the north of the river Clyde in the Greater Glasgow and Clyde health board. These data  are used in the CARBayes vignette.
}
\usage{data(spatialrespdata)}
\format{
  A SpatialPolygonsDataFrame object where the dataframe contains 134 observations on the following 3 variables.
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
data(spatialrespdata)
}

\keyword{datasets}
