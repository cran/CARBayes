highlight.borders <- function(border.locations, sfdata)
{
#### This function takes in an n by n matrix where values of zero represent borders to be highlighted
######################################
#### Identify the borders to highlight
######################################
border.temp <- border.locations
border.temp[upper.tri(border.temp)] <- NA
boundary.list <- which(border.temp==0, arr.ind=TRUE)
boundary.dat <- data.frame(area1=boundary.list[ ,1], area2=boundary.list[ ,2])
M <- nrow(boundary.dat)



################################################
#### Add the geometry to the boundary data frame
################################################
boundary.dat$geometry <- rep(NA, M)
   for(j in 1:M)
   {
   intersect.all <- st_intersection(sfdata$geometry[boundary.list[j, ]])
   intersect.type <- sapply(intersect.all, class)
   intersect.final <- intersect.all[intersect.type[2, ] %in% c("LINESTRING",  "MULTILINESTRING", "GEOMETRYCOLLECTION")]
      if(length(intersect.final)>0)   boundary.dat$geometry[j] <- intersect.final
   }
#boundary.dat2 <- boundary.dat[which(!is.na(boundary.dat$geometry)), ]
boundary.final <- st_as_sf(x=boundary.dat)



############################
#### Return the final object
############################
return(boundary.final) 
}

