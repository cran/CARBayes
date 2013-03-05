\name{highlight.borders}
\alias{highlight.borders}

\title{
Create a SpatialPoints object identifying a subset of borders between neighbouring areas, which allows them to be overlaid on a spatial map.
}


\description{
Create a SpatialPoints object identifying a subset of borders between neighbouring areas, which allows them to be overlaid on a spatial map. The borders that are identified could relate to neighbouring areas that exhibit very different values. See the main CARBayes help file for an example.
}

\usage{
highlight.borders(border.locations, ID, shp, dbf)
}

\arguments{
  \item{border.locations}{
An n by n matrix where n is the number of areas. The matrix contains 3 distinct values: NA for non-neighbouring areas, 0 for borders between neighbouring areas that are to be highlighted on the map, and 1 for borders between neighbouring areas that are not to be highlighted.
}
  \item{ID}{
A vector containing the row names of the data frame which are unique identifiers for each area.
}
  \item{shp}{
The .shp part of a shapefile containing the polygons of the areas that the data relate to.
}
  \item{dbf}{
The .dbf part of the shapefile containing a lookup table, whose first column includes the row names of the data frame. 
}
}



\value{
\item{SpatialPoints }{A SpatialPoints object from the 'sp' package, which contains the vertices of all the borders to be highlighted on the map. The mapping can be done using the spplot() function, see the main help file in CARBayes for an example.
}
}

\author{
Duncan Lee
}


\examples{
## See the examples in the main help file for CARBayes.
}