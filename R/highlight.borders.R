highlight.borders <- function(border.locations, sfdata)
{
#### This function takes in an n by n matrix where values of one represent non-borders and
#### values equal to zero represent borders. The function then links this matrix with an sf
#### object to identify the borders. The result is a multilinestring (output=lines) or 
#### multipoint (output=points) object which can be plotted.

#############################################
#### Create the object to save the results in
#############################################
n <- nrow(border.locations)
border.locations[is.na(border.locations)] <- 2
n.boundaries <- sum(border.locations==0) / 2
boundary.all <- array(c(NA, NA), c(1,2))   



#####################################
#### Identify and save the boundaries
#####################################
     for(i in 1:n)
     {
          for(j in 1:n)
          {
                if(border.locations[i,j]==0 & i>j)
                {
                #### Obtain the intersection points
                b1 <- st_boundary(sfdata$geometry[i])
                b2 <- st_boundary(sfdata$geometry[j])
                intersect.points <- st_intersection(x=b1, y=b2)
                intersect.points2 <- st_cast(x=intersect.points, to="MULTIPOINT") 
                intersect.points3 <- st_coordinates(intersect.points2)[ ,1:2]
                
                #### Save them to the boundary object
                boundary.all <- rbind(boundary.all, intersect.points3)
                }else
                {
                }
           }
     }



###############################
#### Create a MULTIPOINT object
###############################
boundary.all <- boundary.all[-1, ]
boundary.all2 <- boundary.all[!duplicated(boundary.all), ]
boundary.final <- st_multipoint(x=boundary.all2)
return(boundary.final) 
}

