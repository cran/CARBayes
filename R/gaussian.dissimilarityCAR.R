gaussian.dissimilarityCAR <-
function(formula, data=NULL, beta=NULL, phi=NULL, nu2=NULL, tau2=NULL, rho=0.99, alpha=NULL, W, Z, burnin=0, n.sample=1000, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.nu2=NULL, prior.tau2=NULL)
{
cat("Setting up the model\n")
a<-proc.time()
##############################################
#### Format the arguments and check for errors
##############################################
#### Overall formula object
frame <- try(suppressWarnings(model.frame(formula, data=data, na.action=na.pass)), silent=TRUE)
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
    if(length(rho)!= 1) stop("rho is the wrong length.", call.=FALSE)
    if(sum(is.na(rho))>0) stop("rho has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(rho)) stop("rho has non-numeric values.", call.=FALSE)
    if(rho < 0.5 | rho >=1) stop("rho is outside the interval [0.5,1).", call.=FALSE)

## Covariance parameters alpha
     if(is.null(alpha)) alpha <- runif(n=q, min=rep(0,q), max=(alpha.max/(2+q)))
    if(length(alpha)!= q) stop("alpha is the wrong length.", call.=FALSE)
    if(sum(is.na(alpha))>0) stop("alpha has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(alpha)) stop("alpha has non-numeric values.", call.=FALSE)


#### Priors
## Put in default priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(1000, p)
    if(is.null(prior.tau2)) prior.tau2 <- c(0.001, 0.001)
    if(is.null(prior.nu2)) prior.nu2 <- c(0.001, 0.001)

    
## Checks    
    if(length(prior.mean.beta)!=p) stop("the vector of prior means for beta is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.mean.beta)) stop("the vector of prior means for beta is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.mean.beta))!=0) stop("the vector of prior means for beta has missing values.", call.=FALSE)    
 
    if(length(prior.var.beta)!=p) stop("the vector of prior variances for beta is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.var.beta)) stop("the vector of prior variances for beta is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.var.beta))!=0) stop("the vector of prior variances for beta has missing values.", call.=FALSE)    
    if(min(prior.var.beta) <=0) stop("the vector of prior variances has elements less than zero", call.=FALSE)
 
    if(length(prior.tau2)!=2) stop("the prior value for tau2 is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.tau2)) stop("the prior value for tau2 is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.tau2))!=0) stop("the prior value for tau2 has missing values.", call.=FALSE)    
   
    if(length(prior.nu2)!=2) stop("the prior value for nu2 is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.nu2)) stop("the prior value for nu2 is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.nu2))!=0) stop("the prior value for nu2 has missing values.", call.=FALSE)    


#### MCMC quantities
## Checks
    if(!is.numeric(burnin)) stop("burn-in is not a number", call.=FALSE)
    if(!is.numeric(n.sample)) stop("n.sample is not a number", call.=FALSE)    
    if(n.sample <= 0) stop("n.sample is less than or equal to zero.", call.=FALSE)
    if(burnin < 0) stop("burn-in is less than zero.", call.=FALSE)
    if(n.sample <= burnin)  stop("Burn-in is greater than n.sample.", call.=FALSE)
    if(!is.numeric(thin)) stop("thin is not a number", call.=FALSE)
    if(thin <= 0) stop("thin is less than or equal to zero.", call.=FALSE)


## Matrices to store samples
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, p))
samples.phi <- array(NA, c(n.keep, n))
samples.nu2 <- array(NA, c(n.keep, 1))
samples.tau2 <- array(NA, c(n.keep, 1))
samples.alpha <- array(NA, c(n.keep, q))
samples.deviance <- array(NA, c(n.keep, 1))


     
## Metropolis quantities
accept <- c(0,0)
proposal.sd.alpha <- 0.02 * alpha.max
tau2.posterior.shape <- prior.tau2[1] + 0.5*(n-1)
nu2.posterior.shape <- prior.nu2[1] + 0.5*n


#### Checks for the original W matrix
    if(!is.matrix(W)) stop("W is not a matrix.", call.=FALSE)
    if(nrow(W)!= n) stop("W has the wrong number of rows.", call.=FALSE)
    if(ncol(W)!= n) stop("W has the wrong number of columns.", call.=FALSE)
    if(sum(is.na(W))>0) stop("W has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(W)) stop("W has non-numeric values.", call.=FALSE)
    if(!sum(names(table(W))==c(0,1))==2) stop("W has non-binary (zero and one) values.", call.=FALSE)


## Ensure the W matrix is symmetric
Wnew <- array(0, c(n,n))
     for(i in 1:n)
     {
          for(j in 1:n)
          {
               if(i>j)
               {
               temp <- W[i,j]
               Wnew[i,j] <- temp
               Wnew[j,i] <- temp
               }else{}
          }
     }
W <- Wnew  
n.neighbours <- apply(W, 2, sum)
spam.W <- as.spam(W)
     
     
## Create the triplet object
W.triplet <- c(NA, NA, NA)
     for(i in 1:n)
     {
          for(j in 1:n)
          {
               if(W[i,j]==1)
               {
               W.triplet <- rbind(W.triplet, c(i,j, NA))     
               }else{}
          }
     }
W.triplet <- W.triplet[-1, ]     
n.triplet <- nrow(W.triplet) 

     
## Create the start and finish points for W updating
W.begfin <- array(NA, c(n, 2))     
temp <- 1
     for(i in 1:n)
     {
     W.begfin[i, ] <- c(temp, (temp + n.neighbours[i]-1))
     temp <- temp + n.neighbours[i]
     }
     
     
## Create the Z triplet form
Z.triplet <- array(NA, c(n.triplet, q))
     for(i in 1:n.triplet)
     {
     row <- W.triplet[i,1]
     col <- W.triplet[i,2]
          for(j in 1:q)
          {
          Z.triplet[i,j] <- Z[[j]][row, col]     
          }     
     }
W.triplet[ ,3] <- as.numeric(exp(-Z.triplet %*% alpha)>=0.5)
spam.W@entries <- W.triplet[ ,3]      
spam.Wprop <- spam.W     
W.tripletprop <- W.triplet
     
     
#### Create the matrix form of Q     
Q <- -rho * spam.W 
diag(Q) <- rho * rowSums(spam.W) + 1-rho
det.Q <- sum(log(diag(chol.spam(Q))))     


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
## Start timer
cat("Collecting", n.sample, "samples\n", sep = " ")
progressBar <- txtProgressBar(style = 3)
percentage.points<-round((1:100/100)*n.sample)

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
          nu2.posterior.scale <- prior.nu2[2] + 0.5 * sum((Y - fitted.current)^2)
          nu2 <- 1 / rgamma(1, nu2.posterior.shape, scale=(1/nu2.posterior.scale))    
		


		####################
		## Sample from phi
		####################
          offset.phi <- (Y - as.numeric(X.standardised %*% beta) - offset) / nu2    
          phi <- gaussiandissimilaritycarupdate(Wtriplet=W.triplet, Wbegfin=W.begfin, nsites=n, phi=phi, tau2=tau2, rho=rho, nu2=nu2, offset=offset.phi)
          phi <- phi - mean(phi)
		
    
    
		##################
		## Sample from tau2
		##################
          temp2 <- quadformW(Wtriplet=W.triplet, Wbegfin=W.begfin, n_triplet=n.triplet, nsites=n, phi=phi, rho=rho)
          tau2.posterior.scale <- temp2 + prior.tau2[2] 
          tau2 <- 1 / rgamma(1, tau2.posterior.shape, scale=(1/tau2.posterior.scale))		
    
 
        	######################
		#### Sample from alpha
		######################
    	     ## Propose a value
          proposal.alpha <- alpha
    	          for(r in 1:q)
    	          {
    	          proposal.alpha[r] <- rtrunc(n=1, spec="norm", a=0, b=alpha.max[r],  mean=alpha[r], sd=proposal.sd.alpha[r])
    	          }
               
          ## Create the proposal values for W and Q
          W.tripletprop[ ,3] <- as.numeric(exp(-Z.triplet %*% proposal.alpha)>=0.5)
          spam.Wprop@entries <- W.tripletprop[ ,3]     
          Qprop <- -rho * spam.Wprop 
          diag(Qprop) <- rho * rowSums(spam.Wprop) + 1-rho
          det.Qprop <- sum(log(diag(chol.spam(Qprop))))     
          temp3 <- quadformW(Wtriplet=W.tripletprop, Wbegfin=W.begfin, n_triplet=n.triplet, nsites=n, phi=phi, rho=rho)
               
          #### Calculate the acceptance probability
    	     logprob.current <- det.Q - temp2 / tau2
    	     logprob.proposal <- det.Qprop - temp3 / tau2
          prob <- exp(logprob.proposal - logprob.current)
    	     
    	     #### Accept or reject the proposed value
    	          if(prob > runif(1))
    	          {
    	          alpha <- proposal.alpha
    	          det.Q <- det.Qprop 
               W.triplet[ ,3] <- W.tripletprop[ ,3] 
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

    
    
          ################################       
          ## print progress to the console
          ################################
               if(j %in% percentage.points)
               {
               setTxtProgressBar(progressBar, j/n.sample)
               }
   		}


# end timer
cat("\nSummarising results")
close(progressBar)
###################################
#### Summarise and save the results 
###################################
## Acceptance rates
accept.alpha <- 100 * accept[1] / accept[2]
accept.final <- c(100, 100, 100, 100, accept.alpha)
names(accept.final) <- c("beta", "phi", "nu2", "tau2", "alpha")
 
        
     
## Deviance information criterion (DIC)
median.beta <- apply(samples.beta, 2, median)
median.phi <- apply(samples.phi, 2, median)
fitted.median <- X.standardised %*% median.beta + median.phi + offset
nu2.median <- median(samples.nu2)
deviance.fitted <- -2 * sum(dnorm(Y, mean = fitted.median, sd = rep(sqrt(nu2.median),n), log = TRUE))
p.d <- mean(samples.deviance) - deviance.fitted
DIC <- 2 * mean(samples.deviance) - deviance.fitted


#### Compute the Conditional Predictive Ordinate
CPO.temp <- array(NA, c(nrow(samples.phi), n))
    for(i in 1:nrow(samples.phi))
    {
    temp.fitted <- samples.phi[i, ] + X.standardised %*% samples.beta[i, ] + offset
    CPO.temp[i, ] <- 1 / dnorm(x=Y, mean=temp.fitted, sd=sqrt(samples.nu2[i,1]))
    }
CPO <- 1/apply(CPO.temp, 2, mean)
MPL <- sum(log(CPO)) 
     

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





## Compile and return the results
model.string <- c("Likelihood model - Gaussian (identity link function)", "\nRandom effects model - Localised CAR", "\nDissimilarity metrics - ", rownames(summary.alpha), "\n")     
samples <- list(beta=samples.beta.orig, phi=mcmc(samples.phi), tau2=mcmc(samples.tau2), nu2=mcmc(samples.nu2), alpha=mcmc(samples.alpha))
W.summary <- list(W.posterior=W.posterior, W.border.prob=W.border.prob)
results <- list(formula=formula, samples=samples, fitted.values=fitted.values, random.effects=random.effects, residuals=residuals, W.summary=W.summary, DIC=DIC, p.d=p.d,  MPL=MPL, summary.results=summary.results, model=model.string, accept=accept.final)
class(results) <- "carbayes"
b<-proc.time()
cat(" finished in ", round(b[3]-a[3], 1), "seconds")   
return(results)
}
