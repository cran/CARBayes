poisson.independent <-
function(formula, beta=NULL, theta=NULL, sigma2=NULL, burnin=0, n.sample=1000, blocksize.beta=5, blocksize.theta=10, prior.mean.beta=NULL, prior.var.beta=NULL, prior.max.sigma2=NULL)
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

    if(!is.numeric(blocksize.beta)) stop("blocksize.beta is not a number", call.=FALSE)
    if(blocksize.beta <= 0) stop("blocksize.beta is less than or equal to zero", call.=FALSE)
    if(!(floor(blocksize.beta)==ceiling(blocksize.beta))) stop("blocksize.beta has non-integer values.", call.=FALSE)
     if(!is.numeric(blocksize.theta)) stop("blocksize.theta is not a number", call.=FALSE)
    if(blocksize.theta <= 0) stop("blocksize.theta is less than or equal to zero", call.=FALSE)
    if(!(floor(blocksize.theta)==ceiling(blocksize.theta))) stop("blocksize.theta has non-integer values.", call.=FALSE)

## Matrices to store samples
samples.beta <- array(NA, c((n.sample-burnin), p))
samples.theta <- array(NA, c((n.sample-burnin), n))
samples.sigma2 <- array(NA, c((n.sample-burnin), 1))
samples.deviance <- array(NA, c((n.sample-burnin), 1))

## Metropolis quantities
accept.all <- rep(0,4)
accept <- accept.all
proposal.sd.beta <- 0.01
proposal.sd.theta <- 0.01
proposal.corr.beta <- solve(t(X.standardised) %*% X.standardised)


#### Priors
## Put in default priors
## N(0, 100) for beta 
## U(0, 10) for sigma2
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(1000, p)
    if(is.null(prior.max.sigma2)) prior.max.sigma2 <- 1000

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


###########################
#### Run the Bayesian model
###########################
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
        fitted.proposal <- exp(as.numeric(X.standardised %*% proposal.beta) + theta + offset)
        fitted.current <- exp(as.numeric(X.standardised %*% beta) + theta + offset)

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
    ## Sample from theta
    ####################
    #### Create the blocking structure
    if(blocksize.theta >= n)
    {
    n.block <- 1
    beg <- 1
    fin <- n
    }else
    {
    init <- sample(1:blocksize.theta,  1)
    n.standard <- floor((n-init) / blocksize.theta)
    remainder <- n - (init + n.standard * blocksize.theta)
        
        if(n.standard==0)
        {
        beg <- c(1,(init+1))
        fin <- c(init,n)
        }else if(remainder==0)
        {
        beg <- c(1,seq((init+1), n, blocksize.theta))
        fin <- c(init, seq((init+blocksize.theta), n, blocksize.theta))
        }else
        {
        beg <- c(1, seq((init+1), n, blocksize.theta))
        fin <- c(init, seq((init+blocksize.theta), n, blocksize.theta), n)
        }
    n.block <- length(beg)
    }
    

    #### Update the parameters in blocks
    proposal.theta <- theta    
    beta.offset <- X.standardised %*% beta + offset
        
        for(r in 1:n.block)
          {
        ## Propose a value
        n.current <- length(beg[r]:fin[r])
        proposal.theta[beg[r]:fin[r]] <- rnorm(n=n.current, mean=theta[beg[r]:fin[r]], sd=rep(proposal.sd.theta, n.current))
        fitted.proposal <- exp(beta.offset[beg[r]:fin[r]] + proposal.theta[beg[r]:fin[r]])
        fitted.current <- exp(beta.offset[beg[r]:fin[r]] + theta[beg[r]:fin[r]])

        ## Calculate the acceptance probability
        prob1 <- sum(Y[beg[r]:fin[r]] * (log(fitted.proposal) - log(fitted.current)) + fitted.current - fitted.proposal)
        prob2 <- sum(theta[beg[r]:fin[r]]^2  - proposal.theta[beg[r]:fin[r]]^2) / (2 * sigma2)
        prob <- exp(prob1 + prob2)

        ## Accept or reject the value
            if(prob > runif(1))
            {
            theta[beg[r]:fin[r]] <- proposal.theta[beg[r]:fin[r]]
            accept[3] <- accept[3] + 1  
            accept[4] <- accept[4] + 1 
            }else
            {
            proposal.theta[beg[r]:fin[r]] <- theta[beg[r]:fin[r]]
            accept[4] <- accept[4] + 1 
            }
        }
    theta <- theta - mean(theta)
    
    
    
    ####################
    ## Sample from sigma2
    ####################
    sigma2 <- rinvgamma(n=1, shape=(0.5*n-1), scale=(0.5*sum(theta^2)))
            while(sigma2 > prior.max.sigma2)
            {
            sigma2 <- rinvgamma(n=1, shape=(0.5*n-1), scale=(0.5*sum(theta^2)))
            }
    
            
    
    #########################
    ## Calculate the deviance
    #########################
    fitted <- exp(as.numeric(X.standardised %*% beta) + theta + offset)
    deviance <- -2 * sum(Y * log(fitted) -  fitted - lfactorial(Y))



    ###################
    ## Save the results
    ###################
        if(j > burnin)
        {
        ele <- j - burnin
        samples.beta[ele, ] <- beta
        samples.theta[ele, ] <- theta
        samples.sigma2[ele, ] <- sigma2
        samples.deviance[ele, ] <- deviance
        }else
        {
        }


    ########################################
    ## Self tune the acceptance probabilties
    ########################################
    k <- j/100
        if(ceiling(k)==floor(k))
        {
        #### Determine the acceptance probabilities
        accept.beta <- 100 * accept[1] / accept[2]
        accept.theta <- 100 * accept[3] / accept[4]
        accept.all <- accept.all + accept
        accept <- c(0,0,0,0)
            
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
            
        #### theta tuning parameter
            if(accept.theta > 40)
            {
            proposal.sd.theta <- 2 * proposal.sd.theta
            }else if(accept.theta < 30)              
            {
            proposal.sd.theta <- 0.5 * proposal.sd.theta
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



###################################
#### Summarise and save the results 
###################################
## Deviance information criterion (DIC)
median.beta <- apply(samples.beta, 2, median)
median.theta <- apply(samples.theta, 2, median)
fitted.median <- exp(X.standardised %*% median.beta + median.theta + offset)
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

summary.hyper <- quantile(samples.sigma2, c(0.5, 0.025, 0.975))
summary.hyper <- c(summary.hyper, (n.sample-burnin), as.numeric(100 * (1-rejectionRate(mcmc(samples.sigma2)))))

summary.results <- rbind(summary.beta, summary.hyper)
rownames(summary.results)[nrow(summary.results)] <- "sigma2"
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
colnames(fitted.values) <- c("Mean", "Sd", "Median", "2.5%", "97.5%")
fitted.temp <- array(NA, c(nrow(samples.beta), n))
    for(i in 1:nrow(samples.beta))
    {
    fitted.temp[i, ] <- exp(X.standardised %*% samples.beta[i, ] + samples.theta[i, ] + offset)
    }
fitted.values[ ,1] <- apply(fitted.temp, 2, mean)
fitted.values[ ,2] <- apply(fitted.temp, 2, sd)
fitted.values[ ,3:5] <- t(apply(fitted.temp, 2, quantile, c(0.5, 0.025, 0.975)))
fitted.values <- round(fitted.values, 4)


#### Print a summary of the results to the screen
cat("\n#################\n")
cat("#### Model fitted\n")
cat("#################\n\n")
cat("Likelihood model - Poisson (log link function) \n")
cat("Random effects model - Independent\n")
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


## Compile and return the results
results <- list(formula=formula, samples.beta=samples.beta.orig, samples.theta=mcmc(samples.theta), samples.sigma2=mcmc(samples.sigma2), fitted.values=fitted.values, random.effects=random.effects, residuals=residuals, DIC=DIC, p.d=p.d, summary.results=summary.results)
return(results)
}
