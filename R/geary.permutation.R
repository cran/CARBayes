geary.permutation <-
function(x, W, n.simulation)
{
     #### Implement some simple checks     
     n <- length(x)
     
     if(n<10) stop("x contains less than 10 elements and Geary's C should not be calculated.", call.=FALSE)
     if(sum(is.na(x))>0) stop("x contains missing 'NA' values.", call.=FALSE)
     if(!is.numeric(x)) stop("x has non-numeric values.", call.=FALSE)
     
     if(!is.matrix(W)) stop("W is not a matrix.", call.=FALSE)
     if(nrow(W)!= n) stop("W has the wrong number of rows.", call.=FALSE)
     if(ncol(W)!= n) stop("W has the wrong number of columns.", call.=FALSE)
     if(sum(is.na(W))>0) stop("W has missing 'NA' values.", call.=FALSE)
     if(!is.numeric(W)) stop("W has non-numeric values.", call.=FALSE)
     if(!sum(names(table(W))==c(0,1))==2) stop("W has non-binary (zero and one) values.", call.=FALSE)
     
     if(!is.numeric(n.simulation)) stop("n.simulation is not a number", call.=FALSE)    
     if(n.simulation < 100) stop("n.sample is less than 100, the p-value may not be very accurate.", call.=FALSE)
     
     
     
     #### Calculate Geary's C statistic
     diff <- x - mean(x)
     denominator <- 2 * sum(diff^2) * sum(W)
     numerator <- (n-1) * sum(as.matrix(dist(cbind(x,x), method="maximum", diag=TRUE, upper=TRUE))^2 * W)
     G <- numerator / denominator

     
     #### Compute Geary's C for random permutations
     G.perm <- rep(NA, n.simulation)
     
     for(i in 1:n.simulation)
     {
          x.perm <- sample(x=x, size=n, replace=FALSE)
          G.perm[i] <- (n-1) * sum(as.matrix(dist(cbind(x.perm, x.perm), method="maximum", diag=TRUE, upper=TRUE))^2 * W) / denominator
     }
     
     
     
     #### Compute the p-value	 
     rank <- rank(c(G, G.perm))[1]
     pval <- 1 - punif((n.simulation-rank + 2)/(n.simulation + 1))
     
     
     
     #### Return the results
     results <- list(statistic=G, rank=rank, pvalue=pval)
     return(results)
}
