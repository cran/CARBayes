gaussian.independent <-
function(formula, data=NULL, beta=NULL, theta=NULL, nu2=NULL, sigma2=NULL, burnin=0, n.sample=1000, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.max.nu2=NULL, prior.max.sigma2=NULL)
{
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
    
## Random effects theta
    if(is.null(theta)) theta <- rnorm(n=n, mean=rep(0,n), sd=rep(0.1, n))    
    if(length(theta)!= n) stop("theta is the wrong length.", call.=FALSE)
    if(sum(is.na(theta))>0) stop("theta has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(theta)) stop("theta has non-numeric values.", call.=FALSE)

## Random effects variance sigma2
    if(is.null(sigma2)) sigma2 <- runif(1)
    if(length(sigma2)!= 1) stop("sigma2 is the wrong length.", call.=FALSE)
    if(sum(is.na(sigma2))>0) stop("sigma2 has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(sigma2)) stop("sigma2 has non-numeric values.", call.=FALSE)
    if(sigma2 <= 0) stop("sigma2 is negative or zero.", call.=FALSE)


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
samples.theta <- array(NA, c(n.keep, n))
samples.nu2 <- array(NA, c(n.keep, 1))
samples.sigma2 <- array(NA, c(n.keep, 1))
samples.deviance <- array(NA, c(n.keep, 1))

## Metropolis quantities
nu2.posterior.shape <- 0.5*n-1
sigma2.posterior.shape <- 0.5*n-1

#### Priors
## Put in default priors
## N(0, 100) for beta 
## U(0, 10) for sigma2
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(1000, p)
    if(is.null(prior.max.sigma2)) prior.max.sigma2 <- 1000
    if(is.null(prior.max.nu2)) prior.max.nu2 <- 1000
    
    
## Checks    
    if(length(prior.mean.beta)!=p) stop("the vector of prior means for beta is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.mean.beta)) stop("the vector of prior means for beta is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.mean.beta))!=0) stop("the vector of prior means for beta has missing values.", call.=FALSE)    
 
    if(length(prior.var.beta)!=p) stop("the vector of prior variances for beta is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.var.beta)) stop("the vector of prior variances for beta is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.var.beta))!=0) stop("the vector of prior variances for beta has missing values.", call.=FALSE)    
    if(min(prior.var.beta) <=0) stop("the vector of prior variances has elements less than zero", call.=FALSE)

    if(length(prior.max.sigma2)!=1) stop("the maximum prior value for sigma2 is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.max.sigma2)) stop("the maximum prior value for sigma2 is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.max.sigma2))!=0) stop("the maximum prior value for sigma2 has missing values.", call.=FALSE)    
    if(min(prior.max.sigma2) <=0) stop("the maximum prior value for sigma2 is less than zero", call.=FALSE)

    if(length(prior.max.nu2)!=1) stop("the maximum prior value for nu2 is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.max.nu2)) stop("the maximum prior value for nu2 is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.max.nu2))!=0) stop("the maximum prior value for nu2 has missing values.", call.=FALSE)    
    if(min(prior.max.nu2) <=0) stop("the maximum prior value for nu2 is less than zero", call.=FALSE)


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
    for(j in 1:n.sample)
    {
    ####################
    ## Sample from beta
    ####################
    data.mean.beta <- data.temp.beta %*% (Y - theta - offset)
    U <- chol((prior.precision.beta + data.precision.beta / nu2))
    Uinv <- backsolve(U, diag(rep(1,p)))
    fc.mean.beta <- Uinv %*% (t(Uinv) %*% (prior.precision.beta %*% prior.mean.beta + (data.precision.beta / nu2) %*% data.mean.beta))
    beta <- fc.mean.beta + Uinv %*% rnorm(p)    
    

    
    ##################
    ## Sample from nu2
    ##################
    fitted.current <-  as.numeric(X.standardised %*% beta) + theta + offset
    nu2.posterior.scale <- 0.5 * sum((Y - fitted.current)^2)
    nu2 <- 1/rtrunc(n=1, spec="gamma", a=(1/prior.max.nu2), b=Inf,  shape=nu2.posterior.shape, scale=(1/nu2.posterior.scale))
    


    ####################
    ## Sample from theta
    ####################
    data.mean.theta <- Y - as.numeric(X.standardised %*% beta) - offset    
    fc.variance.scalar <- 1/((1/nu2) + (1/sigma2))
    fc.mean <- fc.variance.scalar * (data.mean.theta / nu2)
    theta <- rnorm(n=n, mean=fc.mean, sd=sqrt(fc.variance.scalar))
    theta <- theta - mean(theta)
    
    
    
    ####################
    ## Sample from sigma2
    ####################
    sigma2.posterior.scale <- 0.5*sum(theta^2)
    sigma2 <- 1/rtrunc(n=1, spec="gamma", a=(1/prior.max.sigma2), b=Inf,  shape=sigma2.posterior.shape, scale=(1/sigma2.posterior.scale))
    
    
            
    #########################
    ## Calculate the deviance
    #########################
    fitted <- as.numeric(X.standardised %*% beta) + theta + offset
    deviance <- -2 * sum(dnorm(Y, mean = fitted, sd = rep(sqrt(nu2),n), log = TRUE))



    ###################
    ## Save the results
    ###################
        if(j > burnin & (j-burnin)%%thin==0)
        {
        ele <- (j - burnin) / thin
        samples.beta[ele, ] <- beta
        samples.theta[ele, ] <- theta
        samples.nu2[ele, ] <- nu2
        samples.sigma2[ele, ] <- sigma2
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


accept.final <- rep(100, 4)
names(accept.final) <- c("beta", "theta", "nu2", "sigma")
     
     
###################################
#### Summarise and save the results 
###################################
## Deviance information criterion (DIC)
median.beta <- apply(samples.beta, 2, median)
median.theta <- apply(samples.theta, 2, median)
fitted.median <- X.standardised %*% median.beta + median.theta + offset
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
summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(100, p))
rownames(summary.beta) <- colnames(X)
colnames(summary.beta) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept")

summary.hyper <- array(NA, c(2 ,5))
summary.hyper[1, 1:3] <- quantile(samples.nu2, c(0.5, 0.025, 0.975))
summary.hyper[1, 4:5] <- c(n.keep,100)
summary.hyper[2, 1:3] <- quantile(samples.sigma2, c(0.5, 0.025, 0.975))
summary.hyper[2, 4:5] <- c(n.keep, 100)

summary.results <- rbind(summary.beta, summary.hyper)
rownames(summary.results)[(nrow(summary.results)-1):nrow(summary.results)] <- c("nu2", "sigma2")
summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
summary.results[ , 4:5] <- round(summary.results[ , 4:5], 1)



#### Create the random effects summary
random.effects <- array(NA, c(n, 5))
colnames(random.effects) <- c("Mean", "Sd", "Median", "2.5%", "97.5%")
random.effects[ ,1] <- apply(samples.theta, 2, mean)
random.effects[ ,2] <- apply(samples.theta, 2, sd)
random.effects[ ,3:5] <- t(apply(samples.theta, 2, quantile, c(0.5, 0.025, 0.975)))
random.effects <- round(random.effects, 4)


#### Create the Fitted values
fitted.values <- array(NA, c(n, 5))
residuals <- array(NA, c(n, 5))
colnames(fitted.values) <- c("Mean", "Sd", "Median", "2.5%", "97.5%")
colnames(residuals) <- c("Mean", "Sd", "Median", "2.5%", "97.5%")
fitted.temp <- array(NA, c(nrow(samples.beta), n))
residuals.temp <- array(NA, c(nrow(samples.beta), n))    
    for(i in 1:nrow(samples.beta))
    {
    temp <- X.standardised %*% samples.beta[i, ] + samples.theta[i, ] + offset
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

     




## Compile and return the results
model.string <- c("Likelihood model - Gaussian (identity link function)", "\nRandom effects model - Independent\n")     
samples <- list(beta=samples.beta.orig, theta=mcmc(samples.theta), sigma2=mcmc(samples.sigma2), nu2=mcmc(samples.nu2))
results <- list(formula=formula, samples=samples, fitted.values=fitted.values, random.effects=random.effects, residuals=residuals, W.summary=NULL, DIC=DIC, p.d=p.d, summary.results=summary.results, model=model.string, accept=accept.final)
class(results) <- "carbayes"
return(results)
}

