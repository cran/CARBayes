combine.data.shapefile <-
function(data, shp, dbf)
{
#### This function merges a data frame with a shapefile to enable spaital plotting
#### The data frame must have rownames which match the shapefile id in the dbf file
#### The data frame can be a subset of the areas in the shapefile
     

###############################################################################
#### Subset the shapefile to only contain areas corresponding to the data frame
#### Then turn each portion of the shapefile into a polygon object
###############################################################################
n <- nrow(data)
polygons <- as.list(rep(NA,n))
names(polygons) <- rownames(data)

     for(i in 1:n)
     {
     #### Select the appropriate shapefile component
	index <- which(dbf$dbf[,1]==names(polygons)[i])
          if(length(index)==0) stop("the identifier a row in the dataframe does not appear in the shapfile.", call.=FALSE)
	shapefile <- shp$shp[[index]]

          
     #### Turn the shapefile into a polygon object
          if(shapefile$num.parts==1)
          {
          temp <- Polygon(shapefile$points)
          polygons[[i]] <- Polygons(list(temp), names(polygons)[i])
          }else
          {
          # If the polygon has multiple areas or holes
          n.parts <- shapefile$num.parts
          results.part <- as.list(rep(0,n.parts))
          breakpoints <- c(shapefile$parts, shapefile$num.points)

               for(k in 1:n.parts)
               {
               start <- breakpoints[k]+1
               end <- breakpoints[(k+1)]
               results.part[[k]] <- Polygon(shapefile$points[start:end, ])
               }
               
          polygons[[i]] <- Polygons(results.part, names(polygons)[i])
          }    
 	}

     
# Merge the dataframe and polygons object
poly <- SpatialPolygons(polygons)
combined.data <- SpatialPolygonsDataFrame(poly, data)
return(combined.data)
}
