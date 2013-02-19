gaussian.dissimilarityCAR <-
function(formula, beta=NULL, phi=NULL, nu2=NULL, tau2=NULL, rho=NULL, fix.rho=FALSE, alpha=NULL, W, Z, burnin=0, n.sample=1000, thin=1, blocksize.phi=10, prior.mean.beta=NULL, prior.var.beta=NULL, prior.max.nu2=NULL, prior.max.tau2=NULL, prior.max.alpha=NULL)
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



#### Dissimilarity metric matrix
## Check the list of dissimilarity metrics is appropriate
    if(class(Z)!="list") stop("Z is not a list object.", call.=FALSE)
    if(sum(is.na(as.numeric(lapply(Z, sum, na.rm=FALSE))))>0) stop("Z contains missing 'NA' values.", call.=FALSE)

q <- length(Z)

	if(sum(as.character(lapply(Z,class))=="matrix")<q) stop("Z contains non-matrix values.", call.=FALSE)
	if(sum(as.numeric(lapply(Z,nrow))==n) <q) stop("Z contains matrices of the wrong size.", call.=FALSE)
	if(sum(as.numeric(lapply(Z,ncol))==n) <q) stop("Z contains matrices of the wrong size.", call.=FALSE)
	if(min(as.numeric(lapply(Z,min)))<0) stop("Z contains negative values.", call.=FALSE)


## Determine the default values for the maximums for alpha and the threshold values to be significant
alpha.max <- rep(NA,q)
alpha.threshold <- rep(NA,q)
	for(k in 1:q)
	{
	Z.crit <- quantile(as.numeric(Z[[k]])[as.numeric(Z[[k]])!=0], 0.5)
	alpha.max[k] <- -log(0.5) / Z.crit
	alpha.threshold[k] <- -log(0.5) / max(Z[[k]])
	}



#### Response variable
## Create the response
Y <- model.response(frame)
    
## Check for errors
    if(sum(is.na(Y))>0) stop("the response has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(Y)) stop("the response variable has non-numeric values.", call.=FALSE)



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
    if(is.null(beta)) beta <- lm(Y~X.standardised-1, offset=offset)$coefficients
    if(length(beta)!= p) stop("beta is the wrong length.", call.=FALSE)
    if(sum(is.na(beta))>0) stop("beta has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(beta)) stop("beta has non-numeric values.", call.=FALSE)

## Data variance nu2
    if(is.null(nu2)) nu2 <- runif(1)
    if(length(nu2)!= 1) stop("nu2 is the wrong length.", call.=FALSE)
    if(sum(is.na(nu2))>0) stop("nu2 has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(nu2)) stop("nu2 has non-numeric values.", call.=FALSE)
    if(nu2 <= 0) stop("nu2 is negative or zero.", call.=FALSE)  

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

## Covariance parameters alpha
    if(is.null(alpha)) alpha <- runif(n=q, min=rep(0,q), max=(alpha.max/(2+q)))
    if(length(alpha)!= q) stop("alpha is the wrong length.", call.=FALSE)
    if(sum(is.na(alpha))>0) stop("alpha has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(alpha)) stop("alpha has non-numeric values.", call.=FALSE)



#### MCMC quantities
## Checks
    if(!is.numeric(burnin)) stop("burn-in is not a number", call.=FALSE)
    if(!is.numeric(n.sample)) stop("n.sample is not a number", call.=FALSE)    
    if(n.sample <= 0) stop("n.sample is less than or equal to zero.", call.=FALSE)
    if(burnin < 0) stop("burn-in is less than zero.", call.=FALSE)
    if(n.sample <= burnin)  stop("Burn-in is greater than n.sample.", call.=FALSE)
    if(!is.numeric(thin)) stop("thin is not a number", call.=FALSE)
    if(thin <= 0) stop("thin is less than or equal to zero.", call.=FALSE)

    if(!is.numeric(blocksize.phi)) stop("blocksize.phi is not a number", call.=FALSE)
    if(blocksize.phi <= 0) stop("blocksize.phi is less than or equal to zero", call.=FALSE)
    if(!(floor(blocksize.phi)==ceiling(blocksize.phi))) stop("blocksize.phi has non-integer values.", call.=FALSE)


## Compute the blocking structure for phi
     if(blocksize.phi >= n)
     {
     n.phi.block <- 1
     phi.beg <- 1
     phi.fin <- n  
     }else
     {
     n.standard <- 1 + floor((n-blocksize.phi) / blocksize.phi)
     remainder <- n - (n.standard * blocksize.phi)
     
          if(remainder==0)
          {
          phi.beg <- c(1,seq((blocksize.phi+1), n, blocksize.phi))
          phi.fin <- c(blocksize.phi, seq((blocksize.phi+blocksize.phi), n, blocksize.phi))
          n.phi.block <- length(phi.beg)
          }else if(remainder==1)
          {
          phi.beg <- c(1, seq((blocksize.phi), n, blocksize.phi))
          phi.fin <- c(blocksize.phi-1, seq((blocksize.phi+blocksize.phi-1), n, blocksize.phi), n)
          n.phi.block <- length(phi.beg)    
          }else
          {
          phi.beg <- c(1, seq((blocksize.phi+1), n, blocksize.phi))
          phi.fin <- c(blocksize.phi, seq((blocksize.phi+blocksize.phi), n, blocksize.phi), n)
          n.phi.block <- length(phi.beg)
          }
     }


## Matrices to store samples
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, p))
samples.phi <- array(NA, c(n.keep, n))
samples.nu2 <- array(NA, c(n.keep, 1))
samples.tau2 <- array(NA, c(n.keep, 1))
samples.alpha <- array(NA, c(n.keep, q))
samples.deviance <- array(NA, c(n.keep, 1))

	if(fix.rho)
	{
     accept <- c(0,0)
	}else
	{
     accept <- c(0,0,0,0)
     accept.all <- accept
	samples.rho <- array(NA, c(n.keep, 1))
	proposal.sd.rho <- 0.02
	}

nu2.posterior.shape <- 0.5*n-1
tau2.posterior.shape <- 0.5*n-1

#### Priors
## Put in default priors
## N(0, 100) for beta 
## U(0, 10) for tau2
## U(0, M_i) for alpha
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(1000, p)
    if(is.null(prior.max.tau2)) prior.max.tau2 <- 1000
    if(is.null(prior.max.nu2)) prior.max.nu2 <- 1000
    if(is.null(prior.max.alpha)) prior.max.alpha <- alpha.max

    
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
    
    if(length(prior.max.nu2)!=1) stop("the maximum prior value for nu2 is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.max.nu2)) stop("the maximum prior value for nu2 is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.max.nu2))!=0) stop("the maximum prior value for nu2 has missing values.", call.=FALSE)    
    if(min(prior.max.nu2) <=0) stop("the maximum prior value for nu2 is less than zero", call.=FALSE)

    if(length(prior.max.alpha)!=q) stop("the vector of prior maximums for alpha is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.max.alpha)) stop("the vector of prior maximums for alpha is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.max.alpha))!=0) stop("the vector of prior maximums for alpha has missing values.", call.=FALSE)    


#### Specify the proposal sd for alpha
proposal.sd.alpha <- 0.02 * prior.max.alpha


#### Checks for the original W matrix
    if(!is.matrix(W)) stop("W is not a matrix.", call.=FALSE)
    if(nrow(W)!= n) stop("W has the wrong number of rows.", call.=FALSE)
    if(ncol(W)!= n) stop("W has the wrong number of columns.", call.=FALSE)
    if(sum(is.na(W))>0) stop("W has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(W)) stop("W has non-numeric values.", call.=FALSE)
    if(!sum(names(table(W))==c(0,1))==2) stop("W has non-binary (zero and one) values.", call.=FALSE)


#### Specify the precision matrix
I.n <- diag(1,n)
Z.combined <- array(0, c(n,n))

	for(r in 1:q)
	{
	Z.combined <- Z.combined + alpha[r] * Z[[r]]
	}
	
W.temp <- array(as.numeric(exp(-Z.combined)>=0.5), c(n,n)) * W
W.star <- -W.temp
diag(W.star) <- apply(W.temp, 1, sum)
Q <- rho * W.star + (1-rho) * I.n
spam.Q <- as.spam(Q)
spam.Q.proposal <- spam.Q
det.Q <- sum(log(diag(chol.spam(spam.Q))))
indices.Q <- which(as.numeric(Q)!=0)


#### Beta update quantities
data.precision.beta <- t(X.standardised) %*% X.standardised
data.var.beta <- solve(data.precision.beta)
data.temp.beta <- data.var.beta %*% t(X.standardised)
	if(length(prior.var.beta)==1)
	{
	prior.precision.beta <- 1 / prior.var.beta
	}else
	{
	prior.precision.beta <- solve(diag(prior.var.beta))
	}



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
		data.mean.beta <- data.temp.beta %*% (Y - phi - offset)
		U <- chol((prior.precision.beta + data.precision.beta / nu2))
		Uinv <- backsolve(U, diag(rep(1,p)))
		fc.mean.beta <- Uinv %*% (t(Uinv) %*% (prior.precision.beta %*% prior.mean.beta + (data.precision.beta / nu2) %*% data.mean.beta))
		beta <- fc.mean.beta + Uinv %*% rnorm(p)
		     


		##################
		## Sample from nu2
		##################
		fitted.current <-  as.numeric(X.standardised %*% beta) + phi + offset
		nu2.posterior.scale <- 0.5 * sum((Y - fitted.current)^2)
		nu2 <- 1/rtrunc(n=1, spec="gamma", a=(1/prior.max.nu2), b=Inf,  shape=nu2.posterior.shape, scale=(1/nu2.posterior.scale))
		


		####################
		## Sample from phi
		####################
		Q.temp <- Q / tau2
		data.mean.phi <- Y - as.numeric(X.standardised %*% beta) - offset    
		data.precision.phi <- diag(rep((1/nu2), blocksize.phi))
		b <- rnorm(n)
		
		     for(r in 1:n.phi.block)
		     {
		     ## Create the prior  mean and variance
		     n.current <- length(phi.beg[r]:phi.fin[r])     
		     Q.current <- Q.temp[phi.beg[r]:phi.fin[r], phi.beg[r]:phi.fin[r]]
		     U <- chol(Q.current)
		     Uinv <- backsolve(U, diag(rep(1,n.current)))
		     block.mean <- - Uinv %*% (t(Uinv) %*% Q.temp[phi.beg[r]:phi.fin[r], -(phi.beg[r]:phi.fin[r])] %*% phi[-(phi.beg[r]:phi.fin[r])])
		     
		     ## Create the full conditional
		     U2 <- chol((Q.current + data.precision.phi[1:n.current, 1:n.current])) 
		     U2inv <- backsolve(U2, diag(rep(1,n.current)))
		     fc.mean.phi <- U2inv %*% t(U2inv) %*% (Q.current %*% block.mean + data.precision.phi[1:n.current, 1:n.current] %*% data.mean.phi[phi.beg[r]:phi.fin[r]])
		     
		     ## Update phi
		     phi[phi.beg[r]:phi.fin[r]] <- fc.mean.phi + U2inv %*% b[phi.beg[r]:phi.fin[r]]
		     }   
		phi <- phi - mean(phi)
		
    
    
		##################
		## Sample from tau2
		##################
		tau2.posterior.scale <- 0.5 * sum(phi * (Q %*% phi))
		tau2 <- 1/rtrunc(n=1, spec="gamma", a=(1/prior.max.tau2), b=Inf,  shape=tau2.posterior.shape, scale=(1/tau2.posterior.scale))
		
    
 
		######################
		#### Sample from alpha
		######################
		proposal.alpha <- alpha
		     for(r in 1:q)
		     {
		     proposal.alpha[r] <- rtrunc(n=1, spec="norm", a=0, b=prior.max.alpha[r],  mean=alpha[r], sd=proposal.sd.alpha[r])
		     }
		
		#### Calculate Q.proposal
		Z.combined.proposal <- array(0, c(n,n))
		     for(r in 1:q)
		     {
		     Z.combined.proposal <- Z.combined.proposal + proposal.alpha[r] * Z[[r]]
		     }
		
		W.temp.proposal <- array(as.numeric(exp(-Z.combined.proposal)>=0.5), c(n,n)) * W
		W.star.proposal <- -W.temp.proposal
		diag(W.star.proposal) <- apply(W.temp.proposal, 1, sum)
		proposal.Q <- rho * W.star.proposal + (1 - rho) * I.n
		spam.Q.proposal@entries <- as.numeric(proposal.Q)[indices.Q]   
		proposal.det.Q <- sum(log(diag(chol.spam(spam.Q.proposal))))
		
		#### Calculate the acceptance probability
		logprob.current <- det.Q - 0.5 * sum(phi * (spam.Q %*% phi)) / tau2
		logprob.proposal <- proposal.det.Q - 0.5 * sum(phi * (spam.Q.proposal %*% phi)) / tau2
		prob <- exp(logprob.proposal - logprob.current)
		
		#### Accept or reject the proposed value
		     if(prob > runif(1))
		     {
		     alpha <- proposal.alpha
		     Q <- proposal.Q
		     det.Q <- proposal.det.Q 
		     spam.Q <- spam.Q.proposal
		     W.star <- W.star.proposal 
		     accept[1] <- accept[1] + 1   
		     }else
		     {
		     }
    accept[2] <- accept[2] + 1   
		
      
             
    	     #########################
    	     ## Calculate the deviance
    	     #########################
    	     fitted <- as.numeric(X.standardised %*% beta) + phi + offset
    	     deviance <- -2 * sum(dnorm(Y, mean = fitted, sd = rep(sqrt(nu2),n), log = TRUE))



    	     ###################
    	     ## Save the results
    	     ###################
               if(j > burnin & (j-burnin)%%thin==0)
               {
               ele <- (j - burnin) / thin
               samples.beta[ele, ] <- beta
               samples.phi[ele, ] <- phi
               samples.nu2[ele, ] <- nu2
               samples.tau2[ele, ] <- tau2
               samples.alpha[ele, ] <- alpha
               samples.deviance[ele, ] <- deviance
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
		data.mean.beta <- data.temp.beta %*% (Y - phi - offset)
		U <- chol((prior.precision.beta + data.precision.beta / nu2))
		Uinv <- backsolve(U, diag(rep(1,p)))
		fc.mean.beta <- Uinv %*% (t(Uinv) %*% (prior.precision.beta %*% prior.mean.beta + (data.precision.beta / nu2) %*% data.mean.beta))
		beta <- fc.mean.beta + Uinv %*% rnorm(p)
		     


		##################
		## Sample from nu2
		##################
		fitted.current <-  as.numeric(X.standardised %*% beta) + phi + offset
		nu2.posterior.scale <- 0.5 * sum((Y - fitted.current)^2)
		nu2 <- 1/rtrunc(n=1, spec="gamma", a=(1/prior.max.nu2), b=Inf,  shape=nu2.posterior.shape, scale=(1/nu2.posterior.scale))
		


		####################
		## Sample from phi
		####################
		Q.temp <- Q / tau2
		data.mean.phi <- Y - as.numeric(X.standardised %*% beta) - offset    
		data.precision.phi <- diag(rep((1/nu2), blocksize.phi))
		b <- rnorm(n)
		
		     for(r in 1:n.phi.block)
		     {
		     ## Create the prior  mean and variance
		     n.current <- length(phi.beg[r]:phi.fin[r])     
		     Q.current <- Q.temp[phi.beg[r]:phi.fin[r], phi.beg[r]:phi.fin[r]]
		     U <- chol(Q.current)
		     Uinv <- backsolve(U, diag(rep(1,n.current)))
		     block.mean <- - Uinv %*% (t(Uinv) %*% Q.temp[phi.beg[r]:phi.fin[r], -(phi.beg[r]:phi.fin[r])] %*% phi[-(phi.beg[r]:phi.fin[r])])
		     
		     ## Create the full conditional
		     U2 <- chol((Q.current + data.precision.phi[1:n.current, 1:n.current])) 
		     U2inv <- backsolve(U2, diag(rep(1,n.current)))
		     fc.mean.phi <- U2inv %*% t(U2inv) %*% (Q.current %*% block.mean + data.precision.phi[1:n.current, 1:n.current] %*% data.mean.phi[phi.beg[r]:phi.fin[r]])
		     
		     ## Update phi
		     phi[phi.beg[r]:phi.fin[r]] <- fc.mean.phi + U2inv %*% b[phi.beg[r]:phi.fin[r]]
		     }   
		phi <- phi - mean(phi)
		
    

		##################
		## Sample from tau2
		##################
		tau2.posterior.scale <- 0.5 * sum(phi * (Q %*% phi))
		tau2 <- 1/rtrunc(n=1, spec="gamma", a=(1/prior.max.tau2), b=Inf,  shape=tau2.posterior.shape, scale=(1/tau2.posterior.scale))
		
    
 
		##################
		## Sample from rho
		##################
		#### Propose a value
		proposal.rho <- rtrunc(n=1, spec="norm", a=0, b=1,  mean=rho, sd=proposal.sd.rho)    
		
		#### Calculate the acceptance probability
		proposal.Q <- proposal.rho * W.star + (1-proposal.rho) * I.n
		spam.Q.proposal@entries <- as.numeric(proposal.Q)[indices.Q]   
		proposal.det.Q <- sum(log(diag(chol.spam(spam.Q.proposal))))
		logprob.current <- det.Q - tau2.posterior.scale / tau2
		logprob.proposal <- proposal.det.Q - 0.5 * sum(phi * (spam.Q.proposal %*% phi)) / tau2
		prob <- exp(logprob.proposal - logprob.current)
		
		#### Accept or reject the proposal
		     if(prob > runif(1))
		     {
		     rho <- proposal.rho
		     Q <- proposal.Q
		     spam.Q <- spam.Q.proposal
		     det.Q <- proposal.det.Q
		     accept[3] <- accept[3] + 1   
		     }else
		     {
		     }
		accept[4] <- accept[4] + 1   
		

   
		######################
		#### Sample from alpha
		######################
		proposal.alpha <- alpha
		     for(r in 1:q)
		     {
		     proposal.alpha[r] <- rtrunc(n=1, spec="norm", a=0, b=prior.max.alpha[r],  mean=alpha[r], sd=proposal.sd.alpha[r])
		     }
		
		#### Calculate Q.proposal
		Z.combined.proposal <- array(0, c(n,n))
		     for(r in 1:q)
		     {
		     Z.combined.proposal <- Z.combined.proposal + proposal.alpha[r] * Z[[r]]
		     }
		
		W.temp.proposal <- array(as.numeric(exp(-Z.combined.proposal)>=0.5), c(n,n)) * W
		W.star.proposal <- -W.temp.proposal
		diag(W.star.proposal) <- apply(W.temp.proposal, 1, sum)
		proposal.Q <- rho * W.star.proposal + (1 - rho) * I.n
		spam.Q.proposal@entries <- as.numeric(proposal.Q)[indices.Q]   
		proposal.det.Q <- sum(log(diag(chol.spam(spam.Q.proposal))))
		
		#### Calculate the acceptance probability
		logprob.current <- det.Q - 0.5 * sum(phi * (spam.Q %*% phi)) / tau2
		logprob.proposal <- proposal.det.Q - 0.5 * sum(phi * (spam.Q.proposal %*% phi)) / tau2
		prob <- exp(logprob.proposal - logprob.current)
		
		#### Accept or reject the proposed value
		     if(prob > runif(1))
		     {
		     alpha <- proposal.alpha
		     Q <- proposal.Q
		     det.Q <- proposal.det.Q 
		     spam.Q <- spam.Q.proposal
		     W.star <- W.star.proposal 
		     accept[1] <- accept[1] + 1   
		     }else
		     {
		     }
		accept[2] <- accept[2] + 1   
		
		
             

    	     #########################
    	     ## Calculate the deviance
    	     #########################
    	     fitted <- as.numeric(X.standardised %*% beta) + phi + offset
    	     deviance <- -2 * sum(dnorm(Y, mean = fitted, sd = rep(sqrt(nu2),n), log = TRUE))



    	     ###################
    	     ## Save the results
    	     ###################
               if(j > burnin & (j-burnin)%%thin==0)
               {
               ele <- (j - burnin) / thin
               samples.beta[ele, ] <- beta
               samples.phi[ele, ] <- phi
               samples.nu2[ele, ] <- nu2
               samples.tau2[ele, ] <- tau2
               samples.rho[ele, ] <- rho
               samples.alpha[ele, ] <- alpha
               samples.deviance[ele, ] <- deviance
               }else
               {
               }

    
               
     #######################################
     #### Update the acceptance rate for rho
     #######################################
      k <- j/100
          if(ceiling(k)==floor(k))
          {
          #### Determine the acceptance probabilities
          accept.rho <- 100 * accept[3] / accept[4]
          accept.all <- accept.all + accept
          accept <- c(0,0,0,0)
                    
               #### rho tuning parameter
               if(accept.rho > 70)
               {
               proposal.sd.rho <- min(2 * proposal.sd.rho, 10)
               }else if(accept.rho < 50)              
               {
               proposal.sd.rho <- 0.5 * proposal.sd.rho
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
		
	



###################################
#### Summarise and save the results 
###################################
## Acceptance rates
     if(fix.rho)
     {
     accept.alpha <- 100 * accept[1] / accept[2]    
     }else
     {
     accept.alpha <- 100 * accept.all[1] / accept.all[2]
     }

## Deviance information criterion (DIC)
median.beta <- apply(samples.beta, 2, median)
median.phi <- apply(samples.phi, 2, median)
fitted.median <- X.standardised %*% median.beta + median.phi + offset
nu2.median <- median(samples.nu2)
deviance.fitted <- -2 * sum(dnorm(Y, mean = fitted.median, sd = rep(sqrt(nu2.median),n), log = TRUE))
p.d <- mean(samples.deviance) - deviance.fitted
DIC <- 2 * mean(samples.deviance) - deviance.fitted




#### transform the parameters back to the origianl covariate scale.
samples.beta.orig <- samples.beta
number.cts <- sum(X.indicator==1)     
if(number.cts>0)
{
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
}else
{
}



#### Create a summary object
samples.beta.orig <- mcmc(samples.beta.orig)
summary.beta <- t(apply(samples.beta.orig, 2, quantile, c(0.5, 0.025, 0.975))) 
summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(100,p))
rownames(summary.beta) <- colnames(X)
colnames(summary.beta) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept")

samples.alpha <- mcmc(samples.alpha)
summary.alpha <- t(apply(samples.alpha, 2, quantile, c(0.5, 0.025, 0.975))) 
summary.alpha <- cbind(summary.alpha, rep(n.keep, q), rep(accept.alpha,q))
colnames(summary.alpha) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept")

	if(!is.null(names(Z)))
 	{
 	rownames(summary.alpha) <- names(Z)
	}else
	{
	names.Z <- rep(NA,q)
		for(j in 1:q)
		{
		names.Z[j] <- paste("Z[[",j, "]]", sep="")
		}	
	rownames(summary.alpha) <- names.Z	
	}


	if(fix.rho)
	{
	summary.hyper <- array(NA, c(2 ,5))
	summary.hyper[1, 1:3] <- quantile(samples.nu2, c(0.5, 0.025, 0.975))
	summary.hyper[1, 4:5] <- c(n.keep, 100)
	summary.hyper[2, 1:3] <- quantile(samples.tau2, c(0.5, 0.025, 0.975))
	summary.hyper[2, 4:5] <- c(n.keep, 100)

	summary.results <- rbind(summary.beta, summary.hyper, summary.alpha)
	alpha.min <- c(rep(NA, (p+2)), alpha.threshold)
	summary.results <- cbind(summary.results, alpha.min)
	rownames(summary.results)[(p+1):(p+2)] <- c("nu2", "tau2")
	summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
	summary.results[ , 4:5] <- round(summary.results[ , 4:5], 1)
	summary.results[ , 6] <- round(summary.results[ , 6], 4)
	}else
	{
  accept.rho <- 100 * accept.all[3] / accept.all[4]  
	summary.hyper <- array(NA, c(3 ,5))
	summary.hyper[1, 1:3] <- quantile(samples.nu2, c(0.5, 0.025, 0.975))
	summary.hyper[1, 4:5] <- c(n.keep, 100)
	summary.hyper[2, 1:3] <- quantile(samples.tau2, c(0.5, 0.025, 0.975))
	summary.hyper[2, 4:5] <- c(n.keep, 100)
	summary.hyper[3, 1:3] <- quantile(samples.rho, c(0.5, 0.025, 0.975))
	summary.hyper[3, 4:5] <- c(n.keep, accept.rho)

	summary.results <- rbind(summary.beta, summary.hyper, summary.alpha)
	alpha.min <- c(rep(NA, (p+3)), alpha.threshold)
	summary.results <- cbind(summary.results, alpha.min)
	rownames(summary.results)[(p+1):(p+3)] <- c("nu2", "tau2", "rho")
	summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
	summary.results[ , 4:5] <- round(summary.results[ , 4:5], 1)
	summary.results[ , 6] <- round(summary.results[ , 6], 4)
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
residuals <- array(NA, c(n, 5))
colnames(fitted.values) <- c("Mean", "Sd", "Median", "2.5%", "97.5%")
colnames(residuals) <- c("Mean", "Sd", "Median", "2.5%", "97.5%")
fitted.temp <- array(NA, c(nrow(samples.beta), n))
residuals.temp <- array(NA, c(nrow(samples.beta), n)) 
    for(i in 1:nrow(samples.alpha))
    {
    temp <- X.standardised %*% samples.beta[i, ] + samples.phi[i, ] + offset
    fitted.temp[i, ] <- temp
    residuals.temp[i, ] <- Y - temp    
    }
fitted.values[ ,1] <- apply(fitted.temp, 2, mean)
fitted.values[ ,2] <- apply(fitted.temp, 2, sd)
fitted.values[ ,3:5] <- t(apply(fitted.temp, 2, quantile, c(0.5, 0.025, 0.975)))
fitted.values <- round(fitted.values, 4)
residuals[ ,1] <- apply(residuals.temp, 2, mean)
residuals[ ,2] <- apply(residuals.temp, 2, sd)
residuals[ ,3:5] <- t(apply(residuals.temp, 2, quantile, c(0.5, 0.025, 0.975)))
residuals <- round(residuals, 4)


     
#### Create the posterior medians for the neighbourhood matrix W
W.posterior <- array(NA, c(n,n))
W.border.prob <- array(NA, c(n,n))
	for(i in 1:n)
	{
		for(j in 1:n)
		{
			if(W[i,j]==1)
			{
			z.temp <- NA
				for(k in 1:q)
				{
				z.temp <- c(z.temp, Z[[k]][i,j])
				}	
			z.temp <- z.temp[-1]
			w.temp <- exp(-samples.alpha %*% z.temp)
			w.posterior <- as.numeric(w.temp>=0.5)
			W.posterior[i,j] <- ceiling(median(w.posterior))
			W.border.prob[i,j] <- (1 - sum(w.posterior) / length(w.posterior))
			}else
			{
			}	
		}	
	}



#### Print a summary of the results to the screen
	if(fix.rho)
	{
	cat("\n#################\n")
	cat("#### Model fitted\n")
	cat("#################\n\n")
	cat("Likelihood model - Gaussian (identity link function) \n")
	cat("Random effects model - Localised CAR binary weights\n")
	cat("Regression equation - ")
	print(formula)
	cat("Dissimilarity metrics - ")
	cat(rownames(summary.alpha), sep=", ")

	cat("\n\n############\n")
	cat("#### Results\n")
	cat("############\n\n")

	cat("Posterior quantiles and acceptance rates\n\n")
	print(summary.results)
	cat("\n\n")
	cat("The global spatial correlation parameter rho is fixed at ", rho,"\n\n", sep="")
	cat("Acceptance rate for the random effects is ", 100, "%","\n\n", sep="")
	cat("DIC = ", DIC, "     ", "p.d = ", p.d, "\n")
	cat("\n")
	}else
	{
	cat("\n#################\n")
	cat("#### Model fitted\n")
	cat("#################\n\n")
	cat("Likelihood model - Gaussian (identity link function) \n")
	cat("Random effects model - Localised CAR binary weights\n")
	cat("Regression equation - ")
	print(formula)
	cat("Dissimilarity metrics - ")
	cat(rownames(summary.alpha), sep=", ")

	cat("\n\n############\n")
	cat("#### Results\n")
	cat("############\n\n")

	cat("Posterior quantiles and acceptance rates\n\n")
	print(summary.results)
	cat("\n\n")
	cat("Acceptance rate for the random effects is ", 100, "%","\n\n", sep="")
	cat("DIC = ", DIC, "     ", "p.d = ", p.d, "\n")
	cat("\n")
     }

## Compile and return the results
	if(fix.rho)
	{
	results <- list(formula=formula, samples.beta=samples.beta.orig, samples.phi=mcmc(samples.phi), samples.nu2=mcmc(samples.nu2), samples.tau2=mcmc(samples.tau2), samples.alpha=mcmc(samples.alpha), fitted.values=fitted.values, random.effects=random.effects, W.posterior=W.posterior, W.border.prob=W.border.prob, residuals=residuals, DIC=DIC, p.d=p.d, summary.results=summary.results)
	}else
	{
	results <- list(formula=formula, samples.beta=samples.beta.orig, samples.phi=mcmc(samples.phi), samples.nu2=mcmc(samples.nu2), samples.tau2=mcmc(samples.tau2), samples.rho=mcmc(samples.rho), samples.alpha=mcmc(samples.alpha), fitted.values=fitted.values, random.effects=random.effects, W.posterior=W.posterior, W.border.prob=W.border.prob, residuals=residuals, DIC=DIC, p.d=p.d, summary.results=summary.results)	}

return(results)
}