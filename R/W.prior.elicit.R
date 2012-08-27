W.prior.elicit <-
function(phistar, W, ordinate)
{
     ################################################     
     #### Check the elements take on allowable values	
     ################################################
     ## response variable phistar
     if(missing(phistar)) stop("phistar has not been specified.", call.=FALSE)
     if(sum(is.na(phistar))>0) stop("phistar has missing 'NA' values.", call.=FALSE)
     if(!is.numeric(phistar)) stop("phistar has non-numeric values.", call.=FALSE)	
     n <- length(phistar)	
     
     
     #### W matrix
     if(missing(W)) stop("W has not been specified.", call.=FALSE)
     if(!is.matrix(W)) stop("W is not a matrix.", call.=FALSE)
     if(nrow(W)!= n) stop("W has the wrong number of rows.", call.=FALSE)
     if(ncol(W)!= n) stop("W has the wrong number of columns.", call.=FALSE)
     if(sum(is.na(W))>0) stop("W has missing 'NA' values.", call.=FALSE)
     if(!is.numeric(W)) stop("W has non-numeric values.", call.=FALSE)
     if(!sum(names(table(W))==c(0,1))==2) stop("W has non-binary (zero and one) values.", call.=FALSE)
     
     
     ## ordinate	
     if(missing(ordinate)) stop("ordinate has not been specified.", call.=FALSE)
     if(!(ordinate=="geary" | ordinate=="moran")) stop("ordinate must be either 'geary' or 'moran'.", call.=FALSE)
     
     
     
     ####################################
     #### Compute the prior probabilities
     ####################################
     #### Create a lookup table for W
     n.W <- sum(W)/2
     W.lookup <- array(NA, c(n.W, 3))
     colnames(W.lookup) <- c("number", "row", "col")
     W.lookup[ ,1] <- 1:n.W
     current <- 1
     
     for(i in 1:n)
     {		
          neighbours.temp <- which(W[i, ]==1)
          neighbours.final <- neighbours.temp[neighbours.temp>i]	
          num.neighbours <- length(neighbours.final)
          if(num.neighbours==0)
          {
          }else
          {
               W.lookup[current:(current+num.neighbours-1) , 2] <- i
               W.lookup[current:(current+num.neighbours-1) , 3] <- neighbours.final
               current <- current + num.neighbours	
          }
     }
     
     
     #### Compute the prior probabilities
     W.prior <- array(0, c(n,n))
     W.upper <- upper.tri(W, diag=FALSE)
     
     
     if(ordinate=="geary")
     {
          Diff.geary.mat <- as.matrix(dist(cbind(phistar,phistar), method="maximum", diag=TRUE, upper=TRUE))^2
          Diff.geary.all <- Diff.geary.mat[W.upper]
          
          for(i in 1:n.W)
          {
               row <- W.lookup[i, 2]
               col <- W.lookup[i, 3]	
               prob.geary <-  length(which(Diff.geary.mat[row, col] < Diff.geary.all)) / length(Diff.geary.all)	
               W.prior[row, col] <- prob.geary
               W.prior[col, row] <- prob.geary
          }
     }else
     {
          Diff.moran.mat <- (phistar-mean(phistar)) %*% t(phistar-mean(phistar))	
          Diff.moran.all <- Diff.moran.mat[W.upper]
          
          for(i in 1:n.W)
          {
               row <- W.lookup[i, 2]
               col <- W.lookup[i, 3]	
               prob.moran <-  length(which(Diff.moran.mat[row, col] > Diff.moran.all)) / length(Diff.moran.all)		
               W.prior[row, col] <- prob.moran
               W.prior[col, row] <- prob.moran
          }		
     }
     
     
     #######################
     #### Return the results
     #######################
     return(W.prior)	
}
