binomial.bymCAR <-
function(formula, data=NULL, beta=NULL, phi=NULL, theta=NULL, tau2=NULL, sigma2=NULL, trials, W, burnin=0, n.sample=1000, thin=1, blocksize.beta=5, blocksize.phi=10, blocksize.theta=10, prior.mean.beta=NULL, prior.var.beta=NULL, prior.max.tau2=NULL, prior.max.sigma2=NULL)
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



#### Response variable and trials
## Create the response
Y <- model.response(frame)
    
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
dat <- cbind(Y, trials-Y)

    if(is.null(beta)) beta <- glm(dat~X.standardised-1, offset=offset, family=binomial)$coefficients
    if(length(beta)!= p) stop("beta is the wrong length.", call.=FALSE)
    if(sum(is.na(beta))>0) stop("beta has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(beta)) stop("beta has non-numeric values.", call.=FALSE)

## Random effects phi
    if(is.null(phi)) phi <- rnorm(n=n, mean=rep(0,n), sd=rep(0.1, n))
    if(length(phi)!= n) stop("phi is the wrong length.", call.=FALSE)
    if(sum(is.na(phi))>0) stop("phi has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(phi)) stop("phi has non-numeric values.", call.=FALSE)

## Random effects theta
    if(is.null(theta)) theta <- rnorm(n=n, mean=rep(0,n), sd=rep(0.1, n))
    if(length(theta)!= n) stop("theta is the wrong length.", call.=FALSE)
    if(sum(is.na(theta))>0) stop("theta has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(theta)) stop("theta has non-numeric values.", call.=FALSE)
    
## Random effects variance tau2
    if(is.null(tau2)) tau2 <- runif(1)
    if(length(tau2)!= 1) stop("tau2 is the wrong length.", call.=FALSE)
    if(sum(is.na(tau2))>0) stop("tau2 has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(tau2)) stop("tau2 has non-numeric values.", call.=FALSE)
    if(tau2 <= 0) stop("tau2 is negative or zero.", call.=FALSE)

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

    if(!is.numeric(blocksize.beta)) stop("blocksize.beta is not a number", call.=FALSE)
    if(blocksize.beta <= 0) stop("blocksize.beta is less than or equal to zero", call.=FALSE)
    if(!(floor(blocksize.beta)==ceiling(blocksize.beta))) stop("blocksize.beta has non-integer values.", call.=FALSE)
    if(!is.numeric(blocksize.phi)) stop("blocksize.phi is not a number", call.=FALSE)
    if(blocksize.phi <= 0) stop("blocksize.phi is less than or equal to zero", call.=FALSE)
    if(!(floor(blocksize.phi)==ceiling(blocksize.phi))) stop("blocksize.phi has non-integer values.", call.=FALSE)
    if(!is.numeric(blocksize.theta)) stop("blocksize.theta is not a number", call.=FALSE)
    if(blocksize.theta <= 0) stop("blocksize.theta is less than or equal to zero", call.=FALSE)
    if(!(floor(blocksize.theta)==ceiling(blocksize.theta))) stop("blocksize.theta has non-integer values.", call.=FALSE)


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
          beta.fin <- c(blocksize.beta, seq((blocksize.beta+blocksize.beta), p, blocksize.beta))
          n.beta.block <- length(beta.beg)
          }else
          {
          beta.beg <- c(1, seq((blocksize.beta+1), p, blocksize.beta))
          beta.fin <- c(blocksize.beta, seq((blocksize.beta+blocksize.beta), p, blocksize.beta), p)
          n.beta.block <- length(beta.beg)
          }
     }         


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
samples.tau2 <- array(NA, c(n.keep, 1))
samples.theta <- array(NA, c(n.keep, n))
samples.sigma2 <- array(NA, c(n.keep, 1))
samples.deviance <- array(NA, c(n.keep, 1))

## Metropolis quantities
accept.all <- rep(0,6)
accept <- accept.all
proposal.sd.beta <- 0.01
proposal.sd.phi <- 0.1
proposal.sd.theta <- 0.1
proposal.corr.beta <- solve(t(X.standardised) %*% X.standardised)
chol.proposal.corr.beta <- chol(proposal.corr.beta) 
tau2.posterior.shape <- 0.5 * n - 1
sigma2.posterior.shape <- 0.5 * n - 1


#### Priors
## Put in default priors
## N(0, 100) for beta 
## U(0, 10) for variance parameters
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(1000, p)
    if(is.null(prior.max.tau2)) prior.max.tau2 <- 1000
    if(is.null(prior.max.sigma2)) prior.max.sigma2 <- 1000
    
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
    
    if(length(prior.max.sigma2)!=1) stop("the maximum prior value for sigma2 is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.max.sigma2)) stop("the maximum prior value for sigma2 is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.max.sigma2))!=0) stop("the maximum prior value for sigma2 has missing values.", call.=FALSE)    
    if(min(prior.max.sigma2) <=0) stop("the maximum prior value for sigma2 is less than zero", call.=FALSE)


#### CAR quantities
    if(!is.matrix(W)) stop("W is not a matrix.", call.=FALSE)
    if(nrow(W)!= n) stop("W has the wrong number of rows.", call.=FALSE)
    if(ncol(W)!= n) stop("W has the wrong number of columns.", call.=FALSE)
    if(sum(is.na(W))>0) stop("W has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(W)) stop("W has non-numeric values.", call.=FALSE)
    if(min(W)<0) stop("W has negative elements.", call.=FALSE)
    if(sum(W!=t(W))>0) stop("W is not symmetric.", call.=FALSE)

n.neighbours <- as.numeric(apply(W, 1, sum))
Q <- diag(n.neighbours)  - W


## quantities required in updating phi              
block.mean.part <- as.list(rep(0,n.phi.block))
block.var.chol <- as.list(rep(0,n.phi.block))

     for(r in 1:n.phi.block)
     {
     Q.current <- Q[phi.beg[r]:phi.fin[r], phi.beg[r]:phi.fin[r]]
     block.var <- chol2inv(chol(Q.current))
     block.mean.part[[r]] <- - block.var %*% Q[phi.beg[r]:phi.fin[r], -(phi.beg[r]:phi.fin[r])]
     block.var.chol[[r]] <- chol(block.var)
     }


#### Other quantities needed for the MCMC algorithm
failures <- trials - Y



###########################
#### Run the Bayesian model
###########################
    for(j in 1:n.sample)
    {
    ####################
    ## Sample from beta
    ####################
    proposal <- beta + (sqrt(proposal.sd.beta)* t(chol.proposal.corr.beta)) %*% rnorm(p)
    proposal.beta <- beta    
    phi.offset <- phi + theta + offset
         
         for(r in 1:n.beta.block)
         {
         ## Propose a value
         proposal.beta[beta.beg[r]:beta.fin[r]] <- proposal[beta.beg[r]:beta.fin[r]]
         logit.proposal <- as.numeric(X.standardised %*% proposal.beta) + phi.offset
         logit.current <- as.numeric(X.standardised %*% beta) + phi.offset    
         prob.proposal <- exp(logit.proposal)  / (1 + exp(logit.proposal))
         prob.current <- exp(logit.current)  / (1 + exp(logit.current))
              
         ## Calculate the acceptance probability
         prob1 <- sum(Y * (log(prob.proposal) - log(prob.current)) + failures * (log(1-prob.proposal) - log(1-prob.current)))          
         prob2 <- sum(((beta[beta.beg[r]:beta.fin[r]] - prior.mean.beta[beta.beg[r]:beta.fin[r]])^2 - (proposal.beta[beta.beg[r]:beta.fin[r]] - prior.mean.beta[beta.beg[r]:beta.fin[r]])^2) / prior.var.beta[beta.beg[r]:beta.fin[r]])
         prob <- exp(prob1 + prob2)
              
              ## Accept or reject the value
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
    ## Sample from phi
    ####################
    Q.temp <- Q / tau2
    beta.offset <- as.numeric(X.standardised %*% beta) + theta + offset        
    b <- rnorm(n)
    
         for(r in 1:n.phi.block)
         {
         ## Propose a value
         Q.current <- Q.temp[phi.beg[r]:phi.fin[r], phi.beg[r]:phi.fin[r]]
         block.mean <- block.mean.part[[r]] %*% phi[-(phi.beg[r]:phi.fin[r])]
         proposal.phi <- phi[phi.beg[r]:phi.fin[r]] + (sqrt(proposal.sd.phi) * sqrt(tau2) * t(block.var.chol[[r]])) %*% b[phi.beg[r]:phi.fin[r]]
         logit.proposal <- beta.offset[phi.beg[r]:phi.fin[r]] + proposal.phi
         logit.current <- beta.offset[phi.beg[r]:phi.fin[r]] + phi[phi.beg[r]:phi.fin[r]]    
         prob.proposal <- exp(logit.proposal)  / (1 + exp(logit.proposal))
         prob.current <- exp(logit.current)  / (1 + exp(logit.current))
         
         ## Calculate the acceptance probability
         prob1 <- sum(Y[phi.beg[r]:phi.fin[r]] * (log(prob.proposal) - log(prob.current)) + failures[phi.beg[r]:phi.fin[r]] * (log(1-prob.proposal) - log(1-prob.current)))          
         prob2 <- t(phi[phi.beg[r]:phi.fin[r]] - block.mean) %*% Q.current %*% (phi[phi.beg[r]:phi.fin[r]] - block.mean) - t(proposal.phi - block.mean) %*% Q.current %*% (proposal.phi - block.mean)
         prob <- exp(prob1 + 0.5 * prob2)
         
         ## Accept or reject the value
              if(prob > runif(1))
              {
              phi[phi.beg[r]:phi.fin[r]] <- proposal.phi
              accept[3] <- accept[3] + 1  
              }else
              {
              }
         }                      
    accept[4] <- accept[4] + n.phi.block
    phi <- phi - mean(phi)
    
    

    ####################
    ## Sample from theta
    ####################
    proposal.theta <- rnorm(n, mean=theta, sd=proposal.sd.theta)
    beta.offset <- as.numeric(X.standardised %*% beta) + phi + offset        
    logit.proposal <- beta.offset + proposal.theta
    logit.current <- beta.offset + theta    
    prob.proposal <- exp(logit.proposal)  / (1 + exp(logit.proposal))
    prob.current <- exp(logit.current)  / (1 + exp(logit.current))
    
         for(r in 1:n.phi.block)
         {
         ## Calculate the acceptance probability
         prob1 <- sum(Y[phi.beg[r]:phi.fin[r]] * (log(prob.proposal[phi.beg[r]:phi.fin[r]]) - log(prob.current[phi.beg[r]:phi.fin[r]])) + failures[phi.beg[r]:phi.fin[r]] * (log(1-prob.proposal[phi.beg[r]:phi.fin[r]]) - log(1-prob.current[phi.beg[r]:phi.fin[r]])))          
         prob2 <- sum(theta[phi.beg[r]:phi.fin[r]]^2  - proposal.theta[phi.beg[r]:phi.fin[r]]^2) / sigma2
         prob <- exp(prob1 + 0.5 * prob2)
         
         ## Accept or reject the value
              if(prob > runif(1))
              {
              theta[phi.beg[r]:phi.fin[r]] <- proposal.theta[phi.beg[r]:phi.fin[r]]
              accept[5] <- accept[5] + 1  
              }else
              {
              }
         }
    accept[6] <- accept[6] + n.phi.block
    theta <- theta - mean(theta)
    


    ###################
    ## Sample from tau2
    ###################
    tau2.posterior.scale <- 0.5 * sum(phi * (Q %*% phi))
    tau2 <- 1/rtrunc(n=1, spec="gamma", a=(1/prior.max.tau2), b=Inf,  shape=tau2.posterior.shape, scale=(1/tau2.posterior.scale))
    
    
    
    #####################
    ## Sample from sigma2
    #####################
    sigma2.posterior.scale <- 0.5*sum(theta^2)
    sigma2 <- 1/rtrunc(n=1, spec="gamma", a=(1/prior.max.sigma2), b=Inf,  shape=sigma2.posterior.shape, scale=(1/sigma2.posterior.scale))
    
    
    
    #########################
    ## Calculate the deviance
    #########################
    logit <- as.numeric(X.standardised %*% beta) + phi + theta + offset    
    prob <- exp(logit)  / (1 + exp(logit))
    deviance <- -2 * sum(dbinom(x=Y, size=trials, prob=prob, log=TRUE))



    ###################
    ## Save the results
    ###################
        if(j > burnin & (j-burnin)%%thin==0)
        {
        ele <- (j - burnin) / thin
        samples.beta[ele, ] <- beta
        samples.phi[ele, ] <- phi
        samples.tau2[ele, ] <- tau2
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
        accept.phi <- 100 * accept[3] / accept[4]
        accept.theta <- 100 * accept[5] / accept[6]
        accept.all <- accept.all + accept
        accept <- c(0,0,0,0,0,0)
            
        #### beta tuning parameter
            if(accept.beta > 70)
            {
            proposal.sd.beta <- 2 * proposal.sd.beta
            }else if(accept.beta < 50)              
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
## Acceptance rates
accept.beta <- 100 * accept.all[1] / accept.all[2]
accept.phi <- 100 * accept.all[3] / accept.all[4]
accept.theta <- 100 * accept.all[5] / accept.all[6]
accept.tau2 <- 100
accept.sigma2 <- 100
accept.final <- c(accept.beta, accept.phi, accept.theta, accept.tau2, accept.sigma2)
names(accept.final) <- c("beta", "phi", "theta", "tau2", "sigma2")

## Deviance information criterion (DIC)
median.beta <- apply(samples.beta, 2, median)
median.phi <- apply(samples.phi, 2, median)
median.theta <- apply(samples.theta, 2, median)
median.logit <- as.numeric(X.standardised %*% median.beta) + median.phi + median.theta + offset    
median.prob <- exp(median.logit)  / (1 + exp(median.logit))
fitted.median <- trials * median.prob
deviance.fitted <- -2 * sum(dbinom(x=Y, size=trials, prob=median.prob, log=TRUE))
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
summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta,p))
rownames(summary.beta) <- colnames(X)
colnames(summary.beta) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept")

summary.hyper <- array(NA, c(2,5))
summary.hyper[1, 1:3] <- quantile(samples.tau2, c(0.5, 0.025, 0.975))
summary.hyper[1, 4:5] <- c(n.keep, accept.tau2)
summary.hyper[2, 1:3] <- quantile(samples.sigma2, c(0.5, 0.025, 0.975))
summary.hyper[2, 4:5] <- c(n.keep, accept.sigma2)

summary.results <- rbind(summary.beta, summary.hyper)
rownames(summary.results)[(nrow(summary.results)-1):(nrow(summary.results))] <- c("tau2", "sigma2")
summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
summary.results[ , 4:5] <- round(summary.results[ , 4:5], 1)



#### Create the random effects summary
random.effects <- array(NA, c(n, 5))
colnames(random.effects) <- c("Mean", "Sd", "Median", "2.5%", "97.5%")
samples.re <- samples.phi + samples.theta
random.effects[ ,1] <- apply(samples.re, 2, mean)
random.effects[ ,2] <- apply(samples.re, 2, sd)
random.effects[ ,3:5] <- t(apply(samples.re, 2, quantile, c(0.5, 0.025, 0.975)))
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
    temp.logit <- X.standardised %*% samples.beta[i, ] + samples.re[i, ] + offset    
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
model.string <- c("Likelihood model - Binomial (logit link function)", "\nRandom effects model - BYM CAR\n")
samples <- list(beta=samples.beta.orig, phi=mcmc(samples.phi), theta=mcmc(samples.theta), tau2=mcmc(samples.tau2), sigma2=mcmc(samples.sigma2))
results <- list(formula=formula, samples=samples, fitted.values=fitted.values, random.effects=random.effects, residuals=residuals, W.summary=W, DIC=DIC, p.d=p.d, summary.results=summary.results, model=model.string, accept=accept.final)
class(results) <- "carbayes"
return(results)
}

