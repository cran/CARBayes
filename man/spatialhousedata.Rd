\name{spatialhousedata}
\alias{spatialhousedata}

\docType{data}

\title{
A SpatialPolygonsDataFrame object (from the sp package) containing property price data for the 271 Intermediate Geographies (IG) in the Greater Glasgow and Clyde health board.
}

\description{
A SpatialPolygonsDataFrame object containing data on median property prices in 2008 and possible covariates for the 270 Intermediate Geographies (IG) in the Greater Glasgow and Clyde health board. These data  are used in the worked example in the main help page of CARBayes.
}
\usage{data(spatialhousedata)}
\format{
  A SpatialPolygonsDataFrame object where the dataframe contains 271 observations on the following 6 variables.
  \describe{
    \item{\code{price}}{Median property price in each IG.}
    \item{\code{crime}}{The crime rate (number of crimes per 10,000 people) in each IG.}
    \item{\code{rooms}}{The median number of rooms in a property in each IG.}
    \item{\code{sales}}{The percentage of properties that sold in each IG in a year.}
    \item{\code{driveshop}}{The average time taken to drive to a shopping centre in minutes.}
    \item{\code{type}}{The predominant property type in each IG with levels: \code{detached}, \code{flat}, \code{semi}, \code{terrace}}
  }
}

\source{
These data were provided by the Scottish Government, and can also be downloaded from http://www.sns.gov.uk/. Further details about these data are also available at the same location.
}

\examples{
data(spatialhousedata)
}

\keyword{datasets}
