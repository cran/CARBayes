gaussian.iarCAR <-
function(formula, beta=NULL, phi=NULL, nu2=NULL, tau2=NULL, W, burnin=0, n.sample=1000, blocksize.phi=10, prior.mean.beta=NULL, prior.var.beta=NULL, prior.max.nu2=NULL, prior.max.tau2=NULL)
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


#### MCMC quantities
## Checks
    if(!is.numeric(burnin)) stop("burn-in is not a number", call.=FALSE)
    if(!is.numeric(n.sample)) stop("n.sample is not a number", call.=FALSE)    
    if(n.sample <= 0) stop("n.sample is less than or equal to zero.", call.=FALSE)
    if(burnin < 0) stop("burn-in is less than zero.", call.=FALSE)
    if(n.sample <= burnin)  stop("Burn-in is greater than n.sample.", call.=FALSE)

    if(!is.numeric(blocksize.phi)) stop("blocksize.phi is not a number", call.=FALSE)
    if(blocksize.phi <= 0) stop("blocksize.phi is less than or equal to zero", call.=FALSE)
    if(!(floor(blocksize.phi)==ceiling(blocksize.phi))) stop("blocksize.phi has non-integer values.", call.=FALSE)

## Matrices to store samples
samples.beta <- array(NA, c((n.sample-burnin), p))
samples.phi <- array(NA, c((n.sample-burnin), n))
samples.nu2 <- array(NA, c((n.sample-burnin), 1))
samples.tau2 <- array(NA, c((n.sample-burnin), 1))
samples.deviance <- array(NA, c((n.sample-burnin), 1))


#### Priors
## Put in default priors
## N(0, 100) for beta 
## U(0, 10) for tau2
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(1000, p)
    if(is.null(prior.max.tau2)) prior.max.tau2 <- 1000
    if(is.null(prior.max.nu2)) prior.max.nu2 <- 1000
    
    
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


#### CAR quantities
    if(!is.matrix(W)) stop("W is not a matrix.", call.=FALSE)
    if(nrow(W)!= n) stop("W has the wrong number of rows.", call.=FALSE)
    if(ncol(W)!= n) stop("W has the wrong number of columns.", call.=FALSE)
    if(sum(is.na(W))>0) stop("W has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(W)) stop("W has non-numeric values.", call.=FALSE)
    if(!sum(names(table(W))==c(0,1))==2) stop("W has non-binary (zero and one) values.", call.=FALSE)

n.neighbours <- as.numeric(apply(W, 1, sum))
Q <- diag(n.neighbours)  - W


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
    #### Calculate the full conditional mean and variance
    data.mean.beta <- data.temp.beta %*% (Y - phi - offset)
    fc.variance.beta <- solve((prior.precision.beta + data.precision.beta / nu2))
    fc.mean.beta <- fc.variance.beta %*% (prior.precision.beta %*% prior.mean.beta + (data.precision.beta / nu2) %*% data.mean.beta)
    
    #### Update beta by Gibbs sampling  
    beta <- mvrnorm(n=1, mu=fc.mean.beta, Sigma=fc.variance.beta)
    
    
    
    ##################
    ## Sample from nu2
    ##################
    fitted.current <-  as.numeric(X.standardised %*% beta) + phi + offset
	 nu2.posterior.scale <- 0.5 * sum((Y - fitted.current)^2)
    nu2 <- rinvgamma(n=1, shape=(0.5*n-1), scale=nu2.posterior.scale)
            while(nu2 > prior.max.nu2)
            {
            nu2 <- rinvgamma(n=1, shape=(0.5*n-1), scale=nu2.posterior.scale)
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
    data.mean.phi <- Y - as.numeric(X.standardised %*% beta) - offset    
        
        for(r in 1:n.block)
        {
        ## Create the prior and data means and variances
        prior.precision.phi <- as.matrix(Q.temp[beg[r]:fin[r], beg[r]:fin[r]])
        block.length <- nrow(prior.precision.phi)
        prior.var.phi <- chol2inv(chol(prior.precision.phi))
        prior.mean.phi <- - prior.var.phi %*% Q.temp[beg[r]:fin[r], -(beg[r]:fin[r])] %*% phi[-(beg[r]:fin[r])]
        data.precision.phi <- diag(rep((1/nu2),(block.length+1)))
		  data.precision.phi <- data.precision.phi[1:block.length, 1:block.length]

        ## Create the full conditional
        fc.variance.phi <- solve((prior.precision.phi + data.precision.phi))
        fc.mean.phi <- fc.variance.phi %*% (prior.precision.phi %*% prior.mean.phi + data.precision.phi %*% data.mean.phi[beg[r]:fin[r]])
        
        ## Update phi
        phi[beg[r]:fin[r]] <- mvrnorm(n=1, mu=fc.mean.phi, Sigma=fc.variance.phi)
        }
        
    phi <- phi - mean(phi)
    
    

    ##################
    ## Sample from tau2
    ##################
    tau2.posterior.scale <- 0.5 * t(phi) %*% Q %*% phi
    tau2 <- rinvgamma(n=1, shape=(0.5*(n-3)), scale=tau2.posterior.scale)
            while(tau2 > prior.max.tau2)
            {
            tau2 <- rinvgamma(n=1, shape=(0.5*(n-3)), scale=tau2.posterior.scale)
            }
    
            
    
    #########################
    ## Calculate the deviance
    #########################
    fitted <- as.numeric(X.standardised %*% beta) + phi + offset
    deviance <- -2 * sum(dnorm(Y, mean = fitted, sd = rep(sqrt(nu2),n), log = TRUE))



    ###################
    ## Save the results
    ###################
        if(j > burnin)
        {
        ele <- j - burnin
        samples.beta[ele, ] <- beta
        samples.phi[ele, ] <- phi
        samples.nu2[ele, ] <- nu2
        samples.tau2[ele, ] <- tau2
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



###################################
#### Summarise and save the results 
###################################
## Deviance information criterion (DIC)
median.beta <- apply(samples.beta, 2, median)
median.phi <- apply(samples.phi, 2, median)
fitted.median <- X.standardised %*% median.beta + median.phi + offset
nu2.median <- median(samples.nu2)
deviance.fitted <- -2 * sum(dnorm(Y, mean = fitted.median, sd = rep(sqrt(nu2.median),n), log = TRUE))
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

summary.hyper <- array(NA, c(2 ,5))
summary.hyper[1, 1:3] <- quantile(samples.nu2, c(0.5, 0.025, 0.975))
summary.hyper[1, 4:5] <- c((n.sample-burnin), as.numeric(100 * (1-rejectionRate(mcmc(samples.nu2)))))
summary.hyper[2, 1:3] <- quantile(samples.tau2, c(0.5, 0.025, 0.975))
summary.hyper[2, 4:5] <- c((n.sample-burnin), as.numeric(100 * (1-rejectionRate(mcmc(samples.tau2)))))

summary.results <- rbind(summary.beta, summary.hyper)
rownames(summary.results)[(nrow(summary.results)-1):nrow(summary.results)] <- c("nu2", "tau2")
summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
summary.results[ , 4:5] <- round(summary.results[ , 4:5], 1)



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
    fitted.temp[i, ] <- X.standardised %*% samples.beta[i, ] + samples.phi[i, ] + offset
    }
fitted.values[ ,1] <- apply(fitted.temp, 2, mean)
fitted.values[ ,2] <- apply(fitted.temp, 2, sd)
fitted.values[ ,3:5] <- t(apply(fitted.temp, 2, quantile, c(0.5, 0.025, 0.975)))
fitted.values <- round(fitted.values, 4)



#### Print a summary of the results to the screen
cat("\n#################\n")
cat("#### Model fitted\n")
cat("#################\n\n")
cat("Likelihood model - Gaussian (identity link function) \n")
cat("Random effects model - Intrinsic CAR\n")
cat("Regression equation - ")
print(formula)

cat("\n\n############\n")
cat("#### Results\n")
cat("############\n\n")

cat("Posterior quantiles and acceptance rates\n\n")
print(summary.results)
cat("\n\n")
cat("Acceptance rate for the random effects is ", 100, "%","\n\n", sep="")
cat("DIC = ", DIC, "     ", "p.d = ", p.d, "\n")


## Compile and return the results
results <- list(formula=formula, samples.beta=samples.beta.orig, samples.phi=mcmc(samples.phi), samples.nu2=mcmc(samples.nu2), samples.tau2=mcmc(samples.tau2), fitted.values=fitted.values, random.effects=random.effects, residuals=residuals, DIC=DIC, p.d=p.d, summary.results=summary.results)
return(results)
}
