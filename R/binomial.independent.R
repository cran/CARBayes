binomial.independent <- function(formula, data=NULL, beta=NULL, theta=NULL, sigma2=NULL, trials, burnin=0, n.sample=1000, thin=1, blocksize.beta=5, prior.mean.beta=NULL, prior.var.beta=NULL, prior.sigma2=NULL)
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



#### Response variable and trials
## Create the response
Y <- model.response(frame)
failures <- trials - Y
     
## Check for errors
    if(sum(is.na(trials))>0) stop("the numbers of trials has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(trials)) stop("the numbers of trials has non-numeric values.", call.=FALSE)
int.check <- n-sum(ceiling(trials)==floor(trials))
    if(int.check > 0) stop("the numbers of trials has non-integer values.", call.=FALSE)
    if(min(trials)<=0) stop("the numbers of trials has zero or negative values.", call.=FALSE)

    if(sum(is.na(Y))>0) stop("the response has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(Y)) stop("the response variable has non-numeric values.", call.=FALSE)
int.check <- n-sum(ceiling(Y)==floor(Y))
    if(int.check > 0) stop("the respons variable has non-integer values.", call.=FALSE)
    if(min(Y)<0) stop("the response variable has negative values.", call.=FALSE)
    if(sum(Y>trials)>0) stop("the response variable has larger values that the numbers of trials.", call.=FALSE)



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
dat <- cbind(Y, failures)

    if(is.null(beta)) beta <- glm(dat~X.standardised-1, offset=offset, family=binomial)$coefficients
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
    if(!is.numeric(burnin)) stop("burn-in is not a number", call.=FALSE)
    if(!is.numeric(n.sample)) stop("n.sample is not a number", call.=FALSE)    
    if(n.sample <= 0) stop("n.sample is less than or equal to zero.", call.=FALSE)
    if(burnin < 0) stop("burn-in is less than zero.", call.=FALSE)
    if(n.sample <= burnin)  stop("Burn-in is greater than n.sample.", call.=FALSE)
    if(!is.numeric(thin)) stop("thin is not a number", call.=FALSE)
    if(thin <= 0) stop("thin is less than or equal to zero.", call.=FALSE)

    if(!is.numeric(blocksize.beta)) stop("blocksize.beta is not a number", call.=FALSE)
    if(blocksize.beta <= 0) stop("blocksize.beta is less than or equal to zero", call.=FALSE)
    if(!(floor(blocksize.beta)==ceiling(blocksize.beta))) stop("blocksize.beta has non-integer values.", call.=FALSE)


## Compute the blocking structure for beta     
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
cat("Collecting", n.sample, "samples\n", sep = " ")
progressBar <- txtProgressBar(style = 3)
percentage.points<-round((1:100/100)*n.sample)

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
       prob <- binomialbetaupdate(X.standardised, n, p, beta, proposal.beta, offset.temp, Y, failures, prior.mean.beta, prior.var.beta)
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
    temp1 <- binomialindepupdate(nsites=n, theta=theta, sigma2=sigma2, y=Y, failures=failures, theta_tune=proposal.sd.theta, offset=beta.offset) 
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
    logit <- as.numeric(X.standardised %*% beta) + theta + offset    
    prob <- exp(logit)  / (1 + exp(logit))
    deviance <- -2 * sum(dbinom(x=Y, size=trials, prob=prob, log=TRUE))



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
            if(accept.beta > 50)
            {
            proposal.sd.beta <- 2 * proposal.sd.beta
            }else if(accept.beta < 70)              
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
## Compute the acceptance rates
accept.beta <- 100 * accept.all[1] / accept.all[2]
accept.theta <- 100 * accept.all[3] / accept.all[4]
accept.sigma2 <- 100
accept.final <- c(accept.beta, accept.theta, accept.sigma2)
names(accept.final) <- c("beta", "theta", "sigma2")
     
     
## Deviance information criterion (DIC)
median.beta <- apply(samples.beta, 2, median)
median.theta <- apply(samples.theta, 2, median)
median.logit <- as.numeric(X.standardised %*% median.beta) + median.theta + offset    
median.prob <- exp(median.logit)  / (1 + exp(median.logit))
fitted.median <- trials * median.prob
deviance.fitted <- -2 * sum(dbinom(x=Y, size=trials, prob=median.prob, log=TRUE))
p.d <- mean(samples.deviance) - deviance.fitted
DIC <- 2 * mean(samples.deviance) - deviance.fitted

#### Compute the Conditional Predictive Ordinate
CPO.temp <- array(NA, c(nrow(samples.theta), n))
    for(i in 1:nrow(samples.theta))
    {
    temp.logit <- samples.theta[i, ]  + X.standardised %*% samples.beta[i, ] + offset
    temp <- exp(temp.logit)  / (1 + exp(temp.logit))
    CPO.temp[i, ] <- 1 / dbinom(x=Y, size=trials, prob=temp)
    }
CPO <- 1/apply(CPO.temp, 2, mean)
MPL <- sum(log(CPO))     

     

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
summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta, p))
rownames(summary.beta) <- colnames(X)
colnames(summary.beta) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept")

summary.hyper <- quantile(samples.sigma2, c(0.5, 0.025, 0.975))
summary.hyper <- c(summary.hyper, n.keep, accept.sigma2)

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
residuals <- array(NA, c(n, 5))
colnames(fitted.values) <- c("Mean", "Sd", "Median", "2.5%", "97.5%")
colnames(residuals) <- c("Mean", "Sd", "Median", "2.5%", "97.5%")
fitted.temp <- array(NA, c(nrow(samples.beta), n))
residuals.temp <- array(NA, c(nrow(samples.beta), n))
    for(i in 1:nrow(samples.beta))
    {
    temp.logit <- X.standardised %*% samples.beta[i, ] + samples.theta[i, ] + offset    
    temp <- trials * exp(temp.logit)  / (1 + exp(temp.logit))
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
model.string <- c("Likelihood model - Binomial (logit link function)", "\nRandom effects model - Independent\n")
samples <- list(beta=samples.beta.orig, theta=mcmc(samples.theta), sigma2=mcmc(samples.sigma2))
results <- list(formula=formula, samples=samples, fitted.values=fitted.values, random.effects=random.effects, residuals=residuals, W.summary=NULL, DIC=DIC, p.d=p.d, MPL=MPL, summary.results=summary.results, model=model.string, accept=accept.final)
class(results) <- "carbayes"
b<-proc.time()
cat(" finished in ", round(b[3]-a[3], 1), "seconds")
return(results)
}

