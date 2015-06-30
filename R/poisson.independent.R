poisson.independent <- function(formula, data=NULL, burnin, n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.sigma2=NULL, verbose=TRUE)
{
#### Check on the verbose option
     if(is.null(verbose)) verbose=TRUE     
     if(!is.logical(verbose)) stop("the verbose option is not logical.", call.=FALSE)

     if(verbose)
     {
     cat("Setting up the model\n")
     a<-proc.time()
     }else{}

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
beta <- glm(Y~X.standardised-1, offset=offset, family=poisson)$coefficients
theta <- rnorm(n=n, mean=rep(0,n), sd=rep(0.1, n))    
sigma2 <- runif(1)


#### Priors
## Put in default priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(1000, p)
    if(is.null(prior.sigma2)) prior.sigma2 <- c(0.001, 0.001)
     
## Checks    
    if(length(prior.mean.beta)!=p) stop("the vector of prior means for beta is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.mean.beta)) stop("the vector of prior means for beta is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.mean.beta))!=0) stop("the vector of prior means for beta has missing values.", call.=FALSE)    
 
    if(length(prior.var.beta)!=p) stop("the vector of prior variances for beta is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.var.beta)) stop("the vector of prior variances for beta is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.var.beta))!=0) stop("the vector of prior variances for beta has missing values.", call.=FALSE)    
    if(min(prior.var.beta) <=0) stop("the vector of prior variances has elements less than zero", call.=FALSE)

    if(length(prior.sigma2)!=2) stop("the prior value for sigma2 is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.sigma2)) stop("the prior value for sigma2 is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.sigma2))!=0) stop("the prior value for sigma2 has missing values.", call.=FALSE)    

     
#### MCMC quantities
## Checks
if(is.null(burnin)) stop("the burnin argument is missing", call.=FALSE)
if(is.null(n.sample)) stop("the n.sample argument is missing", call.=FALSE)
if(!is.numeric(burnin)) stop("burn-in is not a number", call.=FALSE)
if(!is.numeric(n.sample)) stop("n.sample is not a number", call.=FALSE) 
if(!is.numeric(thin)) stop("thin is not a number", call.=FALSE)
if(n.sample <= 0) stop("n.sample is less than or equal to zero.", call.=FALSE)
if(burnin < 0) stop("burn-in is less than zero.", call.=FALSE)
if(thin <= 0) stop("thin is less than or equal to zero.", call.=FALSE)
if(n.sample <= burnin)  stop("Burn-in is greater than n.sample.", call.=FALSE)
if(n.sample <= thin)  stop("thin is greater than n.sample.", call.=FALSE)
if(burnin!=round(burnin)) stop("burnin is not an integer.", call.=FALSE) 
if(n.sample!=round(n.sample)) stop("n.sample is not an integer.", call.=FALSE) 
if(thin!=round(thin)) stop("thin is not an integer.", call.=FALSE) 


## Compute the blocking structure for beta     
blocksize.beta <- 5 
     if(blocksize.beta >= p)
     {
     n.beta.block <- 1
     beta.beg <- 1
     beta.fin <- p
     }else
     {
     n.standard <- 1 + floor((p-blocksize.beta) / blocksize.beta)
     remainder <- p - n.standard * blocksize.beta
     
          if(remainder==0)
          {
          beta.beg <- c(1,seq((blocksize.beta+1), p, blocksize.beta))
          beta.fin <- seq(blocksize.beta, p, blocksize.beta)
          n.beta.block <- length(beta.beg)
          }else
          {
          beta.beg <- c(1, seq((blocksize.beta+1), p, blocksize.beta))
          beta.fin <- c(seq((blocksize.beta), p, blocksize.beta), p)
          n.beta.block <- length(beta.beg)
          }
     }         




## Matrices to store samples
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, p))
samples.theta <- array(NA, c(n.keep, n))
samples.sigma2 <- array(NA, c(n.keep, 1))
samples.deviance <- array(NA, c(n.keep, 1))
samples.fitted <- array(NA, c(n.keep, n))

## Metropolis quantities
accept.all <- rep(0,4)
accept <- accept.all
proposal.sd.beta <- 0.01
proposal.sd.theta <- 0.01
proposal.corr.beta <- solve(t(X.standardised) %*% X.standardised)
chol.proposal.corr.beta <- chol(proposal.corr.beta) 
sigma2.posterior.shape <- prior.sigma2[1] + 0.5 * n




###########################
#### Run the Bayesian model
###########################
## Start timer
     if(verbose)
     {
     cat("Collecting", n.sample, "samples\n", sep = " ")
     progressBar <- txtProgressBar(style = 3)
     percentage.points<-round((1:100/100)*n.sample)
     }else
     {
     percentage.points<-round((1:100/100)*n.sample)     
     }
     

     
     for(j in 1:n.sample)
    {
    ####################
    ## Sample from beta
    ####################
    proposal <- beta + (sqrt(proposal.sd.beta)* t(chol.proposal.corr.beta)) %*% rnorm(p)
    proposal.beta <- beta
    offset.temp <- theta + offset

       for(r in 1:n.beta.block)
       {
       proposal.beta[beta.beg[r]:beta.fin[r]] <- proposal[beta.beg[r]:beta.fin[r]]
       prob <- poissonbetaupdate(X.standardised, n, p, beta, proposal.beta, offset.temp, Y, prior.mean.beta, prior.var.beta)
            if(prob > runif(1))
            {
            beta[beta.beg[r]:beta.fin[r]] <- proposal.beta[beta.beg[r]:beta.fin[r]]
            accept[1] <- accept[1] + 1  
            }else
            {
            proposal.beta[beta.beg[r]:beta.fin[r]] <- beta[beta.beg[r]:beta.fin[r]]
            }
        }

    accept[2] <- accept[2] + n.beta.block    

    

    
    ####################
    ## Sample from theta
    ####################
    beta.offset <- as.numeric(X.standardised %*% beta) + offset        
    temp1 <- poissonindepupdate(nsites=n, theta=theta, sigma2=sigma2, y=Y, theta_tune=proposal.sd.theta, offset=beta.offset) 
    theta <- temp1[[1]]
    theta <- theta - mean(theta)    
    accept[3] <- accept[3] + temp1[[2]]
    accept[4] <- accept[4] + n
    
    
    
    ####################
    ## Sample from sigma2
    ####################
    sigma2.posterior.scale <- prior.sigma2[2] + 0.5*sum(theta^2)
    sigma2 <- 1 / rgamma(1, sigma2.posterior.shape, scale=(1/sigma2.posterior.scale))    
    
            
    
    #########################
    ## Calculate the deviance
    #########################
    fitted <- exp(as.numeric(X.standardised %*% beta) + theta + offset)
    deviance.all <- dpois(x=as.numeric(Y), lambda=fitted, log=TRUE)
    deviance <- -2 * sum(deviance.all)  



    ###################
    ## Save the results
    ###################
        if(j > burnin & (j-burnin)%%thin==0)
        {
        ele <- (j - burnin) / thin
        samples.beta[ele, ] <- beta
        samples.theta[ele, ] <- theta
        samples.sigma2[ele, ] <- sigma2
        samples.deviance[ele, ] <- deviance
        samples.fitted[ele, ] <- fitted
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
            }else if(accept.beta < 20)              
            {
            proposal.sd.beta <- 0.5 * proposal.sd.beta
            }else
            {
            }
            
        #### theta tuning parameter
            if(accept.theta > 50)
            {
            proposal.sd.theta <- 2 * proposal.sd.theta
            }else if(accept.theta < 40)              
            {
            proposal.sd.theta <- 0.5 * proposal.sd.theta
            }else
            {
            }
        }else
        {   
        }

    
    
    ################################       
    ## print progress to the console
    ################################
          if(j %in% percentage.points & verbose)
          {
          setTxtProgressBar(progressBar, j/n.sample)
          }
     }


# end timer
     if(verbose)
     {
     cat("\nSummarising results")
     close(progressBar)
     }else
     {}
     
 ###################################
#### Summarise and save the results 
###################################
## Compute the acceptance rates
accept.beta <- 100 * accept.all[1] / accept.all[2]
accept.theta <- 100 * accept.all[3] / accept.all[4]
accept.sigma2 <- 100
accept.final <- c(accept.beta, accept.theta, accept.sigma2)
names(accept.final) <- c("beta", "theta", "sigma2")

## Deviance information criterion (DIC)
median.beta <- apply(samples.beta, 2, median)
median.theta <- apply(samples.theta, 2, median)
fitted.median <- exp(X.standardised %*% median.beta + median.theta + offset)
deviance.fitted <- -2 * sum(dpois(x=Y, lambda=fitted.median, log=TRUE))
p.d <- median(samples.deviance) - deviance.fitted
DIC <- 2 * median(samples.deviance) - deviance.fitted    
     
     
#### Compute the Conditional Predictive Ordinate
CPO <- rep(NA, n)
     for(j in 1:n)
     {
     CPO[j] <- 1/median((1 / dpois(x=Y[j], lambda=samples.fitted[ ,j])))    
     }
LMPL <- sum(log(CPO))  


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
summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta,p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
rownames(summary.beta) <- colnames(X)
colnames(summary.beta) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")


summary.hyper <- quantile(samples.sigma2, c(0.5, 0.025, 0.975))
summary.hyper <- c(summary.hyper, n.keep, accept.sigma2, effectiveSize(samples.sigma2), geweke.diag(samples.sigma2)$z)

summary.results <- rbind(summary.beta, summary.hyper)
rownames(summary.results)[nrow(summary.results)] <- "sigma2"
summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)


#### Create the Fitted values and residuals
fitted.values <- apply(samples.fitted, 2, median)
residuals <- as.numeric(Y) - fitted.values


## Compile and return the results
modelfit <- c(DIC, p.d, LMPL)
names(modelfit) <- c("DIC", "p.d", "LMPL")
model.string <- c("Likelihood model - Poisson (log link function)", "\nRandom effects model - Independent\n")
samples <- list(beta=samples.beta.orig, theta=mcmc(samples.theta), sigma2=mcmc(samples.sigma2), fitted=mcmc(samples.fitted))
results <- list(summary.results=summary.results, samples=samples, fitted.values=fitted.values, residuals=residuals, modelfit=modelfit, accept=accept.final, localised.structure=NULL,  formula=formula, model=model.string, X=X)
class(results) <- "carbayes"

     if(verbose)
     {
     b<-proc.time()
     cat(" finished in ", round(b[3]-a[3], 1), "seconds")
     }else
     {}
return(results)
}

