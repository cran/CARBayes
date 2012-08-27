moran.permutation <-
function(x, W, n.simulation)
{
#### Implement some simple checks	
n <- length(x)

    if(n<10) stop("x contains less than 10 elements and Moran's I should not be calculated.", call.=FALSE)
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
	
	
	
#### Calculate Moran's I statistic
diff <- x - mean(x)
denominator <- sum(diff^2) * sum(W)
numerator <- n * sum(diff %*% t(diff) * W)
I <- numerator / denominator



#### Compute Moran's I for random permutations
I.perm <- rep(NA, n.simulation)

	for(i in 1:n.simulation)
	{
	x.perm <- sample(x=x, size=n, replace=FALSE)
	diff.perm <- x.perm - mean(x.perm)
	I.perm[i] <- n * sum(diff.perm %*% t(diff.perm) * W) / denominator
	}
	
	
	
#### Compute the p-value	 
rank <- rank(c(I, I.perm))[1]
pval <- punif((n.simulation-rank + 2)/(n.simulation + 1))



#### Return the results
results <- list(statistic=I, rank=rank, pvalue=pval)
return(results)
}
