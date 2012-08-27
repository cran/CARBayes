poisson.priorWCAR <-
function(formula, beta=NULL, phi=NULL, tau2=NULL, rho=NULL, fix.rho=FALSE, W, blocksize.W=2, prior.W=NULL, burnin=0, n.sample=1000, blocksize.beta=5, blocksize.phi=10, prior.mean.beta=NULL, prior.var.beta=NULL, prior.max.tau2=NULL)
{
##############################################
#### Format the arguments and check for errors
##############################################
#### Overall formula object
frame <- try(suppressWarnings(model.frame(formula, na.action=na.pass)), silent=TRUE)
if(class(frame)=="try-error") stop("the formula inputted contains an error, e.g the variables may be different lengths.", call.=FALSE)



#### Design matrix
## Create the matrix
X <- try(suppressWarnings(model.matrix(object=attr(frame, "terms"), data=frame)), silent=TRUE)
    if(class(X)=="try-error") stop("the covariate matrix contains inappropriate values.", call.=FALSE)
    if(sum(is.na(X))>0) stop("the covariate matrix contains missing 'NA' values.", call.=FALSE)

n <- nrow(X)
p <- ncol(X)

## Check for linearly related columns
cor.X <- suppressWarnings(cor(X))
diag(cor.X) <- 0

    if(max(cor.X, na.rm=TRUE)==1) stop("the covariate matrix has two exactly linearly related columns.", call.=FALSE)
    if(min(cor.X, na.rm=TRUE)==-1) stop("the covariate matrix has two exactly linearly related columns.", call.=FALSE)
 
 	 if(p>1)
	 {
    	if(sort(apply(X, 2, sd))[2]==0) stop("the covariate matrix has two intercept terms.", call.=FALSE)
	 }else
	 {
	 }
 
## Standardise the matrix
X.standardised <- X
X.sd <- apply(X, 2, sd)
X.mean <- apply(X, 2, mean)
X.indicator <- rep(NA, p)       # To determine which parameter estimates to transform back

    for(j in 1:p)
    {
        if(length(table(X[ ,j]))>2)
        {
        X.indicator[j] <- 1
        X.standardised[ ,j] <- (X[ ,j] - mean(X[ ,j])) / sd(X[ ,j])
        }else if(length(table(X[ ,j]))==1)
        {
        X.indicator[j] <- 2
        }else
        {
        X.indicator[j] <- 0
        }
    }



#### Response variable
## Create the response
Y <- model.response(frame)
    
## Check for errors
    if(sum(is.na(Y))>0) stop("the response has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(Y)) stop("the response variable has non-numeric values.", call.=FALSE)
int.check <- n-sum(ceiling(Y)==floor(Y))
    if(int.check > 0) stop("the respons variable has non-integer values.", call.=FALSE)
    if(min(Y)<0) stop("the response variable has negative values.", call.=FALSE)



#### Offset variable
## Create the offset
offset <- try(model.offset(frame), silent=TRUE)

## Check for errors
    if(class(offset)=="try-error")   stop("the offset is not numeric.", call.=FALSE)
    if(is.null(offset))  offset <- rep(0,n)
    if(sum(is.na(offset))>0) stop("the offset has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(offset)) stop("the offset variable has non-numeric values.", call.=FALSE)



#### Initial parameter values
## Regression parameters beta
	 if(is.null(beta)) beta <- glm(Y~X.standardised-1, offset=offset, family=poisson)$coefficients
    if(length(beta)!= p) stop("beta is the wrong length.", call.=FALSE)
    if(sum(is.na(beta))>0) stop("beta has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(beta)) stop("beta has non-numeric values.", call.=FALSE)

## Random effects phi
    if(is.null(phi)) phi <- rnorm(n=n, mean=rep(0,n), sd=rep(0.1, n))
    if(length(phi)!= n) stop("phi is the wrong length.", call.=FALSE)
    if(sum(is.na(phi))>0) stop("phi has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(phi)) stop("phi has non-numeric values.", call.=FALSE)

## Random effects variance tau2
    if(is.null(tau2)) tau2 <- runif(1)
    if(length(tau2)!= 1) stop("tau2 is the wrong length.", call.=FALSE)
    if(sum(is.na(tau2))>0) stop("tau2 has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(tau2)) stop("tau2 has non-numeric values.", call.=FALSE)
    if(tau2 <= 0) stop("tau2 is negative or zero.", call.=FALSE)    

## Global correlation parameter rho
    if(is.null(rho) & fix.rho==TRUE) stop("rho is fixed yet a value has not been specified.", call.=FALSE)
    if(is.null(rho)) rho <- runif(1)
    if(length(rho)!= 1) stop("rho is the wrong length.", call.=FALSE)
    if(sum(is.na(rho))>0) stop("rho has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(rho)) stop("rho has non-numeric values.", call.=FALSE)
    if(rho < 0 | rho >=1) stop("rho is outside the interval [0,1).", call.=FALSE)
 


#### MCMC quantities
## Checks
    if(!is.numeric(burnin)) stop("burn-in is not a number", call.=FALSE)
    if(!is.numeric(n.sample)) stop("n.sample is not a number", call.=FALSE)    
    if(n.sample <= 0) stop("n.sample is less than or equal to zero.", call.=FALSE)
    if(burnin < 0) stop("burn-in is less than zero.", call.=FALSE)
    if(n.sample <= burnin)  stop("Burn-in is greater than n.sample.", call.=FALSE)

    if(!is.numeric(blocksize.beta)) stop("blocksize.beta is not a number", call.=FALSE)
    if(blocksize.beta <= 0) stop("blocksize.beta is less than or equal to zero", call.=FALSE)
    if(!(floor(blocksize.beta)==ceiling(blocksize.beta))) stop("blocksize.beta has non-integer values.", call.=FALSE)
    if(!is.numeric(blocksize.phi)) stop("blocksize.phi is not a number", call.=FALSE)
    if(blocksize.phi <= 0) stop("blocksize.phi is less than or equal to zero", call.=FALSE)
    if(!(floor(blocksize.phi)==ceiling(blocksize.phi))) stop("blocksize.phi has non-integer values.", call.=FALSE)

    if(!is.numeric(blocksize.W)) stop("blocksize.W is not a number", call.=FALSE)
    if(blocksize.W < 1) stop("blocksize.W is less than 1.", call.=FALSE)
    if(ceiling(blocksize.W)!=floor(blocksize.W))  stop("blocksize.W is not an integer.", call.=FALSE)



## Matrices to store samples
samples.beta <- array(NA, c(n.sample, p))
samples.phi <- array(NA, c(n.sample, n))
samples.tau2 <- array(NA, c(n.sample, 1))
samples.deviance <- array(NA, c(n.sample, 1))

## Metropolis quantities
	if(fix.rho)
	{
	accept <- rep(0,6)
	accept.all <- accept
	}else
	{
	samples.rho <- array(NA, c(n.sample, 1))
	accept <- rep(0,8)
	accept.all <- accept
	proposal.sd.rho <- 0.05
	}

proposal.sd.beta <- 0.01
proposal.sd.phi <- 0.1
proposal.corr.beta <- solve(t(X.standardised) %*% X.standardised)


#### Priors
## Put in default priors
## N(0, 100) for beta 
## U(0, 10) for tau2
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(1000, p)
    if(is.null(prior.max.tau2)) prior.max.tau2 <- 1000

    
## Checks    
    if(length(prior.mean.beta)!=p) stop("the vector of prior means for beta is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.mean.beta)) stop("the vector of prior means for beta is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.mean.beta))!=0) stop("the vector of prior means for beta has missing values.", call.=FALSE)    
 
    if(length(prior.var.beta)!=p) stop("the vector of prior variances for beta is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.var.beta)) stop("the vector of prior variances for beta is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.var.beta))!=0) stop("the vector of prior variances for beta has missing values.", call.=FALSE)    
    if(min(prior.var.beta) <=0) stop("the vector of prior variances has elements less than zero", call.=FALSE)
 
    if(length(prior.max.tau2)!=1) stop("the maximum prior value for tau2 is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.max.tau2)) stop("the maximum prior value for tau2 is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.max.tau2))!=0) stop("the maximum prior value for tau2 has missing values.", call.=FALSE)    
    if(min(prior.max.tau2) <=0) stop("the maximum prior value for tau2 is less than zero", call.=FALSE)



#### Checks for the original W matrix
    if(!is.matrix(W)) stop("W is not a matrix.", call.=FALSE)
    if(nrow(W)!= n) stop("W has the wrong number of rows.", call.=FALSE)
    if(ncol(W)!= n) stop("W has the wrong number of columns.", call.=FALSE)
    if(sum(is.na(W))>0) stop("W has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(W)) stop("W has non-numeric values.", call.=FALSE)
    if(!sum(names(table(W))==c(0,1))==2) stop("W has non-binary (zero and one) values.", call.=FALSE)

#### Checks for the prior.W matrix
	if(is.null(prior.W)) prior.W <- W/2
	if(!is.matrix(prior.W)) stop("prior.W is not a matrix.", call.=FALSE)
    if(nrow(prior.W)!= n) stop("prior.W has the wrong number of rows.", call.=FALSE)
    if(ncol(prior.W)!= n) stop("prior.W has the wrong number of columns.", call.=FALSE)
    if(sum(is.na(prior.W))>0) stop("prior.W has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(prior.W)) stop("prior.W has non-numeric values.", call.=FALSE)
    if(min(prior.W)<0 | max(prior.W)>1) stop("prior.W contains probabilities outside [0,1].", call.=FALSE)
W.variable <- length(which(as.numeric(prior.W)>0 & as.numeric(prior.W)<1))
	if(W.variable==0) stop("prior.W has all zero or one values, it does not need to be updated.", call.=FALSE)


#### Specify a look-up table for saving W
n.W <- sum(W)/2
W.lookup <- array(NA, c(n.W, 4))
colnames(W.lookup) <- c("number", "row", "col", "prior")
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

W.current <- W	
	for(i in 1:n.W)
	{
	W.lookup[i, 4] <- prior.W[W.lookup[i, 2], W.lookup[i,3]]		
		
		if(W.lookup[i,4]==0) 
		{
		W.current[W.lookup[i,2], W.lookup[1,3]] <- 0
		W.current[W.lookup[i,3], W.lookup[1,2]] <- 0
		}else
		{
		}
	}

indices.allowed <- which(W.lookup[ , 4]>0 & W.lookup[ , 4] <1)		
samples.W <- array(NA, c((n.sample+1), n.W))
samples.W[1, ] <- as.numeric(W.lookup[ ,4]!=0)



#### Specify the precision matrix
I.n <- diag(1,n)
W.star <- -W.current
diag(W.star) <- apply(W.current, 1, sum)
Q <- rho * W.star + (1-rho) * I.n
det.Q <- as.numeric(determinant(Q, logarithm=TRUE)$modulus)



###########################
#### Run the Bayesian model
###########################
	if(fix.rho)
	{	
    	for(j in 1:n.sample)
    	{
    	####################
    	## Sample from beta
    	####################
    	#### Create the blocking structure
    		if(blocksize.beta >= p)
    		{
    		n.block <- 1
    		beg <- 1
    		fin <- p
    		}else
    		{
    		init <- sample(1:blocksize.beta,  1)
    		n.standard <- floor((p-init) / blocksize.beta)
    		remainder <- p - (init + n.standard * blocksize.beta)
        
        		if(n.standard==0)
        		{
        		beg <- c(1,(init+1))
        		fin <- c(init,p)
        		}else if(remainder==0)
        		{
        		beg <- c(1,seq((init+1), p, blocksize.beta))
        		fin <- c(init, seq((init+blocksize.beta), p, blocksize.beta))
        		}else
        		{
        		beg <- c(1, seq((init+1), p, blocksize.beta))
        		fin <- c(init, seq((init+blocksize.beta), p, blocksize.beta), p)
        		}
    		n.block <- length(beg)
    		}


    	#### Update the parameters in blocks
    	proposal.beta <- beta    
        
        for(r in 1:n.block)
        {
        ## Propose a value
        n.current <- length(beg[r]:fin[r])
        proposal.beta[beg[r]:fin[r]] <- mvrnorm(n=1, mu=beta[beg[r]:fin[r]], Sigma=(proposal.sd.beta * proposal.corr.beta[beg[r]:fin[r], beg[r]:fin[r]]))
        fitted.proposal <- exp(as.numeric(X.standardised %*% proposal.beta) + phi + offset)
        fitted.current <- exp(as.numeric(X.standardised %*% beta) + phi + offset)

        ## Calculate the acceptance probability
        prob1 <- sum(Y * (log(fitted.proposal) - log(fitted.current)) + fitted.current - fitted.proposal)
        prob2 <- sum(((beta[beg[r]:fin[r]] - prior.mean.beta[beg[r]:fin[r]])^2 - (proposal.beta[beg[r]:fin[r]] - prior.mean.beta[beg[r]:fin[r]])^2) / (2 * prior.var.beta[beg[r]:fin[r]]))
        prob <- exp(prob1 + prob2)

        ## Accept or reject the value
            if(prob > runif(1))
            {
            beta <- proposal.beta
            accept[1] <- accept[1] + 1  
            accept[2] <- accept[2] + 1 
            }else
            {
            proposal.beta <- beta
            accept[2] <- accept[2] + 1 
            }
        }



    	####################
    	## Sample from phi
    	####################
    	#### Create the blocking structure
    		if(blocksize.phi >= n)
    		{
    		n.block <- 1
    		beg <- 1
    		fin <- n
    		}else
    		{
    		init <- sample(1:blocksize.phi,  1)
    		n.standard <- floor((n-init) / blocksize.phi)
    		remainder <- n - (init + n.standard * blocksize.phi)
        
        		if(n.standard==0)
        		{
        		beg <- c(1,(init+1))
        		fin <- c(init,n)
        		}else if(remainder==0)
        		{
        		beg <- c(1,seq((init+1), n, blocksize.phi))
        		fin <- c(init, seq((init+blocksize.phi), n, blocksize.phi))
        		}else
        		{
        		beg <- c(1, seq((init+1), n, blocksize.phi))
        		fin <- c(init, seq((init+blocksize.phi), n, blocksize.phi), n)
        		}
    		n.block <- length(beg)
    		}
    

    	#### Update the parameters in blocks
    	Q.temp <- Q / tau2
    	beta.offset <- as.numeric(X.standardised %*% beta) + offset
    	proposal.phi <- phi    
        
        for(r in 1:n.block)
        {
        ## Propose a value
        Q.current <- Q.temp[beg[r]:fin[r], beg[r]:fin[r]]
        block.var <- chol2inv(chol(Q.current))
        block.mean <- - block.var %*% Q.temp[beg[r]:fin[r], -(beg[r]:fin[r])] %*% phi[-(beg[r]:fin[r])]
        proposal.phi[beg[r]:fin[r]] <- mvrnorm(n=1, mu=phi[beg[r]:fin[r]], Sigma=(proposal.sd.phi * block.var))
        fitted.proposal <- exp(beta.offset[beg[r]:fin[r]] + proposal.phi[beg[r]:fin[r]])
        fitted.current <- exp(beta.offset[beg[r]:fin[r]] + phi[beg[r]:fin[r]])
        
        ## Calculate the acceptance probability
        prob1 <- sum(Y[beg[r]:fin[r]] * (log(fitted.proposal) - log(fitted.current)) + fitted.current - fitted.proposal)
        prob2 <- t(phi[beg[r]:fin[r]] - block.mean) %*% Q.current %*% (phi[beg[r]:fin[r]] - block.mean) - t(proposal.phi[beg[r]:fin[r]] - block.mean) %*% Q.current %*% (proposal.phi[beg[r]:fin[r]] - block.mean)
        prob <- exp(prob1 + 0.5 * prob2)

        ## Accept or reject the value
            if(prob > runif(1))
            {
            phi[beg[r]:fin[r]] <- proposal.phi[beg[r]:fin[r]]
            accept[3] <- accept[3] + 1  
            accept[4] <- accept[4] + 1 
            }else
            {
            proposal.phi[beg[r]:fin[r]] <- phi[beg[r]:fin[r]]
            accept[4] <- accept[4] + 1 
            }
        }
        
    	phi <- phi - mean(phi)
    
    

    	##################
    	## Sample from tau2
    	##################
    	tau2.posterior.scale <- 0.5 * t(phi) %*% Q %*% phi
    	tau2 <- rinvgamma(n=1, shape=(0.5*n-1), scale=tau2.posterior.scale)
            while(tau2 > prior.max.tau2)
            {
            tau2 <- rinvgamma(n=1, shape=(0.5*n-1), scale=tau2.posterior.scale)
            }   
   
   
   
 		##################################################
		## Update blocksize.W elements of W in each iteration      
 		##################################################
		#### Choose blocksize.W elements at random
		indices.change <- sample(W.lookup[indices.allowed, 1], blocksize.W, replace=FALSE)	
	
		#### Compute the proposed W matrix	
		elements.change <- W.lookup[indices.change, 2:3]
		W.proposal <- W.current
		diag(W.proposal[elements.change[ ,1], elements.change[ ,2]]) <- 1 - diag(W.proposal[elements.change[ ,1], elements.change[ ,2]])
		diag(W.proposal[elements.change[ ,2], elements.change[ ,1]]) <- 1 - diag(W.proposal[elements.change[ ,2], elements.change[ ,1]])
					
	 	#### Calculate Q.proposal
	 	W.star.proposal <- -W.proposal
		diag(W.star.proposal) <- apply(W.proposal, 1, sum)
	 	proposal.Q <- rho * W.star.proposal + (1-rho) * I.n
	 	proposal.det.Q <- as.numeric(determinant(proposal.Q, logarithm=TRUE)$modulus)

    	#### Calculate the acceptance probability
    	elements.current <- diag(W.current[elements.change[ ,1], elements.change[ ,2]])
    	prior.current <- rep(NA, blocksize.W)
    	prior.current[elements.current==1] <- W.lookup[indices.change[elements.current==1], 4]
    	prior.current[elements.current==0] <- 1 - W.lookup[indices.change[elements.current==0], 4]
    	prior.proposal <- 1 - prior.current
    	
 		logprob.current <- 0.5 * det.Q - 0.5 * t(phi) %*% Q %*% phi / tau2 + sum(log(prior.current))
    	logprob.proposal <- 0.5 * proposal.det.Q - 0.5 * t(phi) %*% proposal.Q %*% phi / tau2 + sum(log(prior.proposal))

    	prob <- exp(logprob.proposal - logprob.current)
    
    	#### Accept or reject the proposal and save the result
    	samples.W[(j+1), ] <- samples.W[j, ]
    	
            if(prob > runif(1))
            {
            W.current <- W.proposal	
            W.star <- W.star.proposal	
            Q <- proposal.Q
            det.Q <- proposal.det.Q
            accept[5] <- accept[5] + 1
            accept[6] <- accept[6] + 1
        	samples.W[(j+1), indices.change] <- 1 - samples.W[j, indices.change]
            }else
            {
            accept[6] <- accept[6] + 1
            }

 
 
     	#########################
    	## Calculate the deviance
    	#########################
    	fitted <- exp(as.numeric(X.standardised %*% beta) + phi + offset)
    	deviance <- -2 * sum(Y * log(fitted) -  fitted - lfactorial(Y))



    	###################
    	## Save the results
    	###################
        samples.beta[j, ] <- beta
        samples.phi[j, ] <- phi
        samples.tau2[j, ] <- tau2
        samples.deviance[j, ] <- deviance


 
    	########################################
    	## Self tune the acceptance probabilties
    	########################################
    	k <- j/100
        if(ceiling(k)==floor(k))
        {
        #### Determine the acceptance probabilities
        accept.beta <- 100 * accept[1] / accept[2]
        accept.phi <- 100 * accept[3] / accept[4]
        accept.all <- accept.all + accept
        accept <- c(0,0,0,0,0,0)
            
        #### beta tuning parameter
            if(accept.beta > 40)
            {
            proposal.sd.beta <- 2 * proposal.sd.beta
            }else if(accept.beta < 30)              
            {
            proposal.sd.beta <- 0.5 * proposal.sd.beta
            }else
            {
            }
            
        #### phi tuning parameter
            if(accept.phi > 40)
            {
            proposal.sd.phi <- 2 * proposal.sd.phi
            }else if(accept.phi < 30)              
            {
            proposal.sd.phi <- 0.5 * proposal.sd.phi
            }else
            {
            }
        }else
        {   
        }

    
    
    	#######################################
    	#### Print out the number of iterations
    	#######################################
    	k <- j/1000
        if(ceiling(k)==floor(k))
        {
        cat("Completed ",j, " samples\n")
        flush.console()
        }else
        {
        }
		}
	}else
	{
    	for(j in 1:n.sample)
    	{
    	####################
    	## Sample from beta
    	####################
    	#### Create the blocking structure
    		if(blocksize.beta >= p)
    		{
    		n.block <- 1
    		beg <- 1
    		fin <- p
    		}else
    		{
    		init <- sample(1:blocksize.beta,  1)
    		n.standard <- floor((p-init) / blocksize.beta)
    		remainder <- p - (init + n.standard * blocksize.beta)
        
        		if(n.standard==0)
        		{
        		beg <- c(1,(init+1))
        		fin <- c(init,p)
        		}else if(remainder==0)
        		{
        		beg <- c(1,seq((init+1), p, blocksize.beta))
        		fin <- c(init, seq((init+blocksize.beta), p, blocksize.beta))
        		}else
        		{
        		beg <- c(1, seq((init+1), p, blocksize.beta))
        		fin <- c(init, seq((init+blocksize.beta), p, blocksize.beta), p)
        		}
    		n.block <- length(beg)
    		}


    	#### Update the parameters in blocks
    	proposal.beta <- beta    
        
        for(r in 1:n.block)
        {
        ## Propose a value
        n.current <- length(beg[r]:fin[r])
        proposal.beta[beg[r]:fin[r]] <- mvrnorm(n=1, mu=beta[beg[r]:fin[r]], Sigma=(proposal.sd.beta * proposal.corr.beta[beg[r]:fin[r], beg[r]:fin[r]]))
        fitted.proposal <- exp(as.numeric(X.standardised %*% proposal.beta) + phi + offset)
        fitted.current <- exp(as.numeric(X.standardised %*% beta) + phi + offset)

        ## Calculate the acceptance probability
        prob1 <- sum(Y * (log(fitted.proposal) - log(fitted.current)) + fitted.current - fitted.proposal)
        prob2 <- sum(((beta[beg[r]:fin[r]] - prior.mean.beta[beg[r]:fin[r]])^2 - (proposal.beta[beg[r]:fin[r]] - prior.mean.beta[beg[r]:fin[r]])^2) / (2 * prior.var.beta[beg[r]:fin[r]]))
        prob <- exp(prob1 + prob2)

        ## Accept or reject the value
            if(prob > runif(1))
            {
            beta <- proposal.beta
            accept[1] <- accept[1] + 1  
            accept[2] <- accept[2] + 1 
            }else
            {
            proposal.beta <- beta
            accept[2] <- accept[2] + 1 
            }
        }



    	####################
    	## Sample from phi
    	####################
    	#### Create the blocking structure
    		if(blocksize.phi >= n)
    		{
    		n.block <- 1
    		beg <- 1
    		fin <- n
    		}else
    		{
    		init <- sample(1:blocksize.phi,  1)
    		n.standard <- floor((n-init) / blocksize.phi)
    		remainder <- n - (init + n.standard * blocksize.phi)
        
        		if(n.standard==0)
        		{
        		beg <- c(1,(init+1))
        		fin <- c(init,n)
        		}else if(remainder==0)
        		{
        		beg <- c(1,seq((init+1), n, blocksize.phi))
        		fin <- c(init, seq((init+blocksize.phi), n, blocksize.phi))
        		}else
        		{
        		beg <- c(1, seq((init+1), n, blocksize.phi))
        		fin <- c(init, seq((init+blocksize.phi), n, blocksize.phi), n)
        		}
    		n.block <- length(beg)
    		}
    

    	#### Update the parameters in blocks
    	Q.temp <- Q / tau2
    	beta.offset <- as.numeric(X.standardised %*% beta) + offset
    	proposal.phi <- phi    
        
        for(r in 1:n.block)
        {
        ## Propose a value
        Q.current <- Q.temp[beg[r]:fin[r], beg[r]:fin[r]]
        block.var <- chol2inv(chol(Q.current))
        block.mean <- - block.var %*% Q.temp[beg[r]:fin[r], -(beg[r]:fin[r])] %*% phi[-(beg[r]:fin[r])]
        proposal.phi[beg[r]:fin[r]] <- mvrnorm(n=1, mu=phi[beg[r]:fin[r]], Sigma=(proposal.sd.phi * block.var))
        fitted.proposal <- exp(beta.offset[beg[r]:fin[r]] + proposal.phi[beg[r]:fin[r]])
        fitted.current <- exp(beta.offset[beg[r]:fin[r]] + phi[beg[r]:fin[r]])
        
        ## Calculate the acceptance probability
        prob1 <- sum(Y[beg[r]:fin[r]] * (log(fitted.proposal) - log(fitted.current)) + fitted.current - fitted.proposal)
        prob2 <- t(phi[beg[r]:fin[r]] - block.mean) %*% Q.current %*% (phi[beg[r]:fin[r]] - block.mean) - t(proposal.phi[beg[r]:fin[r]] - block.mean) %*% Q.current %*% (proposal.phi[beg[r]:fin[r]] - block.mean)
        prob <- exp(prob1 + 0.5 * prob2)

        ## Accept or reject the value
            if(prob > runif(1))
            {
            phi[beg[r]:fin[r]] <- proposal.phi[beg[r]:fin[r]]
            accept[3] <- accept[3] + 1  
            accept[4] <- accept[4] + 1 
            }else
            {
            proposal.phi[beg[r]:fin[r]] <- phi[beg[r]:fin[r]]
            accept[4] <- accept[4] + 1 
            }
        }
        
    	phi <- phi - mean(phi)
    
    

    	##################
    	## Sample from tau2
    	##################
    	tau2.posterior.scale <- 0.5 * t(phi) %*% Q %*% phi
    	tau2 <- rinvgamma(n=1, shape=(0.5*n-1), scale=tau2.posterior.scale)
            while(tau2 > prior.max.tau2)
            {
            tau2 <- rinvgamma(n=1, shape=(0.5*n-1), scale=tau2.posterior.scale)
            }


    
    	##################
    	## Sample from rho
    	##################
    	#### Propose a value
    	proposal.rho <- rnorm(n=1, mean=rho, sd=proposal.sd.rho)
            while(proposal.rho >= 1 | proposal.rho < 0)
            {
            proposal.rho <- rnorm(n=1, mean=rho, sd=proposal.sd.rho)
            }      
   
	 	#### Calculate Q.proposal
	 	proposal.Q <- proposal.rho * W.star + (1-proposal.rho) * I.n
	 	proposal.det.Q <- as.numeric(determinant(proposal.Q, logarithm=TRUE)$modulus)


    	#### Calculate the acceptance probability
    	logprob.current <- 0.5 * det.Q - tau2.posterior.scale / tau2
    	logprob.proposal <- 0.5 * proposal.det.Q - 0.5 * t(phi) %*% proposal.Q %*% phi / tau2
    	prob <- exp(logprob.proposal - logprob.current)
    
    	#### Accept or reject the proposal
            if(prob > runif(1))
            {
            rho <- proposal.rho
            Q <- proposal.Q
            det.Q <- proposal.det.Q        
            accept[5] <- accept[5] + 1  
            accept[6] <- accept[6] + 1 
            }else
            {
            accept[6] <- accept[6] + 1 
            }

          
 
  		##################################################
		## Update blocksize.W elements of W in each iteration      
 		##################################################
		#### Choose blocksize.W elements at random
		indices.change <- sample(W.lookup[indices.allowed, 1], blocksize.W, replace=FALSE)	
	
		#### Compute the proposed W matrix	
		elements.change <- W.lookup[indices.change, 2:3]
		W.proposal <- W.current
		diag(W.proposal[elements.change[ ,1], elements.change[ ,2]]) <- 1 - diag(W.proposal[elements.change[ ,1], elements.change[ ,2]])
		diag(W.proposal[elements.change[ ,2], elements.change[ ,1]]) <- 1 - diag(W.proposal[elements.change[ ,2], elements.change[ ,1]])
					
	 	#### Calculate Q.proposal
	 	W.star.proposal <- -W.proposal
		diag(W.star.proposal) <- apply(W.proposal, 1, sum)
	 	proposal.Q <- rho * W.star.proposal + (1-rho) * I.n
	 	proposal.det.Q <- as.numeric(determinant(proposal.Q, logarithm=TRUE)$modulus)

    	#### Calculate the acceptance probability
    	elements.current <- diag(W.current[elements.change[ ,1], elements.change[ ,2]])
    	prior.current <- rep(NA, blocksize.W)
    	prior.current[elements.current==1] <- W.lookup[indices.change[elements.current==1], 4]
    	prior.current[elements.current==0] <- 1 - W.lookup[indices.change[elements.current==0], 4]
    	prior.proposal <- 1 - prior.current
    	
 		logprob.current <- 0.5 * det.Q - 0.5 * t(phi) %*% Q %*% phi / tau2 + sum(log(prior.current))
    	logprob.proposal <- 0.5 * proposal.det.Q - 0.5 * t(phi) %*% proposal.Q %*% phi / tau2 + sum(log(prior.proposal))
 		prob <- exp(logprob.proposal - logprob.current)
    
    	#### Accept or reject the proposal and save the result
    	samples.W[(j+1), ] <- samples.W[j, ]
    	
            if(prob > runif(1))
            {
            W.current <- W.proposal	
            W.star <- W.star.proposal	
            Q <- proposal.Q
            det.Q <- proposal.det.Q
            accept[7] <- accept[7] + 1
            accept[8] <- accept[8] + 1
        	samples.W[(j+1), indices.change] <- 1 - samples.W[j, indices.change]
            }else
            {
            accept[8] <- accept[8] + 1
            }



    	#########################
    	## Calculate the deviance
    	#########################
    	fitted <- exp(as.numeric(X.standardised %*% beta) + phi + offset)
    	deviance <- -2 * sum(Y * log(fitted) -  fitted - lfactorial(Y))



    	###################
    	## Save the results
    	###################
    	samples.beta[j, ] <- beta
        samples.phi[j, ] <- phi
        samples.tau2[j, ] <- tau2
        samples.rho[j, ] <- rho
        samples.deviance[j, ] <- deviance


 
    	########################################
    	## Self tune the acceptance probabilties
    	########################################
    	k <- j/100
        if(ceiling(k)==floor(k))
        {
        #### Determine the acceptance probabilities
        accept.beta <- 100 * accept[1] / accept[2]
        accept.phi <- 100 * accept[3] / accept[4]
        accept.all <- accept.all + accept
        accept <- c(0,0,0,0,0,0, 0, 0)
            
        #### beta tuning parameter
            if(accept.beta > 40)
            {
            proposal.sd.beta <- 2 * proposal.sd.beta
            }else if(accept.beta < 30)              
            {
            proposal.sd.beta <- 0.5 * proposal.sd.beta
            }else
            {
            }
            
        #### phi tuning parameter
            if(accept.phi > 40)
            {
            proposal.sd.phi <- 2 * proposal.sd.phi
            }else if(accept.phi < 30)              
            {
            proposal.sd.phi <- 0.5 * proposal.sd.phi
            }else
            {
            }
        }else
        {   
        }

    
    
    	#######################################
    	#### Print out the number of iterations
    	#######################################
    	k <- j/1000
        if(ceiling(k)==floor(k))
        {
        cat("Completed ",j, " samples\n")
        flush.console()
        }else
        {
        }
		}
	}



##############################
#### Remove the burnin samples
##############################
	if(p==1)
	{
	samples.beta <- matrix(samples.beta[(burnin+1):n.sample, ], ncol=1)
	}else
	{
	samples.beta <- samples.beta[(burnin+1):n.sample, ]		
	}

samples.phi <- samples.phi[(burnin+1):n.sample, ]
samples.tau2 <- matrix(samples.tau2[(burnin+1):n.sample, ], ncol=1)
samples.deviance <- matrix(samples.deviance[(burnin+1):n.sample, ], ncol=1)
samples.W <- samples.W[-1, ]
samples.W <- samples.W[(burnin+1):n.sample, ]
	if(fix.rho==FALSE) samples.rho <- matrix(samples.rho[(burnin+1):n.sample, ], ncol=1)




###################################
#### Summarise and save the results 
###################################
## Deviance information criterion (DIC)
median.beta <- apply(samples.beta, 2, median)
median.phi <- apply(samples.phi, 2, median)
fitted.median <- exp(X.standardised %*% median.beta + median.phi + offset)
deviance.fitted <- -2 * sum(Y * log(fitted.median) -  fitted.median - lfactorial(Y))
p.d <- mean(samples.deviance) - deviance.fitted
DIC <- 2 * mean(samples.deviance) - deviance.fitted
residuals <- Y - fitted.median



#### transform the parameters back to the origianl covariate scale.
samples.beta.orig <- samples.beta
    for(r in 1:p)
    {
        if(X.indicator[r]==1)
        {
        samples.beta.orig[ ,r] <- samples.beta[ ,r] / X.sd[r]
        }else if(X.indicator[r]==2 & p>1)
        {
        X.transformed <- which(X.indicator==1)
        samples.temp <- as.matrix(samples.beta[ ,X.transformed])
            for(s in 1:length(X.transformed))
            {
            samples.temp[ ,s] <- samples.temp[ ,s] * X.mean[X.transformed[s]]  / X.sd[X.transformed[s]]
            }
        intercept.adjustment <- apply(samples.temp, 1,sum) 
        samples.beta.orig[ ,r] <- samples.beta[ ,r] - intercept.adjustment
        }else
        {
        }
    }



#### Create a summary object
samples.beta.orig <- mcmc(samples.beta.orig)
summary.beta <- t(apply(samples.beta.orig, 2, quantile, c(0.5, 0.025, 0.975))) 
summary.beta <- cbind(summary.beta, rep((n.sample-burnin), p), as.numeric(100 * (1-rejectionRate(samples.beta.orig))))
rownames(summary.beta) <- colnames(X)
colnames(summary.beta) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept")

	if(fix.rho)
	{
	summary.hyper <- array(NA, c(1 ,5))
	summary.hyper[1, 1:3] <- quantile(samples.tau2, c(0.5, 0.025, 0.975))
	summary.hyper[1, 4:5] <- c((n.sample-burnin), as.numeric(100 * (1-rejectionRate(mcmc(samples.tau2)))))

	summary.results <- rbind(summary.beta, summary.hyper)
	rownames(summary.results)[(p+1)] <- c("tau2")
	summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
	summary.results[ , 4:5] <- round(summary.results[ , 4:5], 1)
	}else
	{
	summary.hyper <- array(NA, c(2 ,5))
	summary.hyper[1, 1:3] <- quantile(samples.tau2, c(0.5, 0.025, 0.975))
	summary.hyper[1, 4:5] <- c((n.sample-burnin), as.numeric(100 * (1-rejectionRate(mcmc(samples.tau2)))))
	summary.hyper[2, 1:3] <- quantile(samples.rho, c(0.5, 0.025, 0.975))
	summary.hyper[2, 4:5] <- c((n.sample-burnin), as.numeric(100 * (1-rejectionRate(mcmc(samples.rho)))))

	summary.results <- rbind(summary.beta, summary.hyper)
	rownames(summary.results)[(p+1):(p+2)] <- c("tau2", "rho")
	summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
	summary.results[ , 4:5] <- round(summary.results[ , 4:5], 1)
	}



#### Create the random effects summary
random.effects <- array(NA, c(n, 5))
colnames(random.effects) <- c("Mean", "Sd", "Median", "2.5%", "97.5%")
random.effects[ ,1] <- apply(samples.phi, 2, mean)
random.effects[ ,2] <- apply(samples.phi, 2, sd)
random.effects[ ,3:5] <- t(apply(samples.phi, 2, quantile, c(0.5, 0.025, 0.975)))
random.effects <- round(random.effects, 4)



#### Create the Fitted values
fitted.values <- array(NA, c(n, 5))
colnames(fitted.values) <- c("Mean", "Sd", "Median", "2.5%", "97.5%")
fitted.temp <- array(NA, c(nrow(samples.beta), n))
    for(i in 1:nrow(samples.beta))
    {
    fitted.temp[i, ] <- exp(X.standardised %*% samples.beta[i, ] + samples.phi[i, ] + offset)
    }
fitted.values[ ,1] <- apply(fitted.temp, 2, mean)
fitted.values[ ,2] <- apply(fitted.temp, 2, sd)
fitted.values[ ,3:5] <- t(apply(fitted.temp, 2, quantile, c(0.5, 0.025, 0.975)))
fitted.values <- round(fitted.values, 4)



#### Create the posterior median and probability of a border for W
W.posterior <- array(0, c(n,n))
W.border.prob <- array(NA, c(n,n))

	for(i in 1:n.W)
	{
    row <- W.lookup[i, 2]
    col <- W.lookup[i, 3]	
	W.single <- samples.W[ ,i]
	W.posterior[row, col] <- ceiling(median(W.single))
	W.posterior[col, row] <- ceiling(median(W.single))
	W.border.prob[row, col] <- 100 * (1 - sum(W.single) / length(W.single))
	W.border.prob[col, row] <- 100 * (1 - sum(W.single) / length(W.single))
	}
n.borders <- sum(abs(W - W.posterior)) / 2

	if(fix.rho)
	{
	accept.W <- 100 * accept.all[5] / accept.all[6]	
	}else
	{
	accept.W <- 100 * accept.all[7] / accept.all[8]			
	}



#### Print a summary of the results to the screen
	if(fix.rho)
	{
	cat("\n#################\n")
	cat("#### Model fitted\n")
	cat("#################\n\n")
	cat("Likelihood model - Poisson (log link function) \n")
	cat("Random effects model - Localised CAR\n")
	cat("Regression equation - ")
	print(formula)

	cat("\n\n############\n")
	cat("#### Results\n")
	cat("############\n\n")

	cat("Posterior quantiles and acceptance rates\n\n")
	print(summary.results)
	cat("\n\n")
	cat("The global spatial correlation parameter rho is fixed at ", rho,"\n\n", sep="")
	cat("Acceptance rate for the random effects is ", round(100 * accept.all[3] / accept.all[4],1), "%","\n\n", sep="")
	cat("DIC = ", DIC, "     ", "p.d = ", p.d, "\n")
	cat("\n")
	cat("There are ", n.borders, " boundaries in the study region based on a 0.5 threshold","\n\n", sep="")
	}else
	{
	cat("\n#################\n")
	cat("#### Model fitted\n")
	cat("#################\n\n")
	cat("Likelihood model - Poisson (log link function) \n")
	cat("Random effects model - Localised CAR\n")
	cat("Regression equation - ")
	print(formula)

	cat("\n\n############\n")
	cat("#### Results\n")
	cat("############\n\n")

	cat("Posterior quantiles and acceptance rates\n\n")
	print(summary.results)
	cat("\n\n")
	cat("Acceptance rate for the random effects is ", round(100 * accept.all[3] / accept.all[4],1), "%","\n\n", sep="")
	cat("DIC = ", DIC, "     ", "p.d = ", p.d, "\n")
	cat("\n")
	cat("There are ", n.borders, " boundaries in the study region based on a 0.5 threshold","\n\n", sep="")
	}


## Compile and return the results
	if(fix.rho)
	{
	results <- list(formula=formula, samples.beta=samples.beta.orig, samples.phi=mcmc(samples.phi), samples.tau2=mcmc(samples.tau2), samples.W=mcmc(samples.W), fitted.values=fitted.values, random.effects=random.effects, W.posterior=W.posterior, W.border.prob=W.border.prob, residuals=residuals, DIC=DIC, p.d=p.d, summary.results=summary.results, accept.W=accept.W)
	}else
	{
	results <- list(formula=formula, samples.beta=samples.beta.orig, samples.phi=mcmc(samples.phi), samples.tau2=mcmc(samples.tau2), samples.rho=mcmc(samples.rho), samples.W=mcmc(samples.W), fitted.values=fitted.values, random.effects=random.effects, W.posterior=W.posterior, W.border.prob=W.border.prob, residuals=residuals, DIC=DIC, p.d=p.d, summary.results=summary.results, accept.W=accept.W)
	}

return(results)
}
