\name{combine.data.shapefile}
\alias{combine.data.shapefile}

\title{
Combine a data frame with a shapefile to create a SpatialPolygonsDataFrame object.
}

\description{
This function combines a data frame with a shapefile and creates a SpatialPolygonsDataFrame object from the 'sp' package. The creation of this object allows the variables in the data frame to be plotted as a spatial map using the 'spplot()' function. For the function to work the row names in the data frame must be a subset of the names in the first column of the dbf file. See the main CARBayes help file for an example.
}

\usage{
combine.data.shapefile(data, shp, dbf)
}

\arguments{
\item{data}{
A data frame containing the variables you wish to map. Note, the row names of this data frame must appear in the .dbf file (see argument 3 below).
}
  \item{shp}{
The .shp part of a shapefile containing the polygons of the areas that the data relate to.
}
  \item{dbf}{
The .dbf part of the shapefile containing a lookup table, whose first column includes the row names of the data frame. 
}
}


\value{
\item{SpatialPolygonsDataFrame }{A SpatialPolygonsDataFrame object from the 'sp' package, which has combined the data frame with the shapefile and allows the variables in the data frame to be plotted as a spatial map using the spplot() function.}
}

\author{
Duncan Lee
}

\examples{
## See the examples in the main help file for CARBayes.
}
