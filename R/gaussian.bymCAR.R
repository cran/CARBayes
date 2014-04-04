gaussian.bymCAR <-
function(formula, data=NULL,  W, burnin=0, n.sample=1000, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.nu2=NULL, prior.tau2=NULL, prior.sigma2=NULL, verbose=TRUE)
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



#### Offset variable
## Create the offset
offset <- try(model.offset(frame), silent=TRUE)

## Check for errors
    if(class(offset)=="try-error")   stop("the offset is not numeric.", call.=FALSE)
    if(is.null(offset))  offset <- rep(0,n)
    if(sum(is.na(offset))>0) stop("the offset has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(offset)) stop("the offset variable has non-numeric values.", call.=FALSE)
    


#### Initial parameter values
beta <- lm(Y~X.standardised-1, offset=offset)$coefficients
nu2 <- runif(1)
phi <- rnorm(n=n, mean=rep(0,n), sd=rep(0.1, n))
theta <- rnorm(n=n, mean=rep(0,n), sd=rep(0.1, n))
tau2 <- runif(1)
sigma2 <- runif(1)


#### Priors
## Put in default priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(1000, p)
    if(is.null(prior.tau2)) prior.tau2 <- c(0.001, 0.001)
    if(is.null(prior.nu2)) prior.nu2 <- c(0.001, 0.001)
    if(is.null(prior.sigma2)) prior.sigma2 <- c(0.001, 0.001)        
        
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



## Matrices to store samples
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, p))
samples.phi <- array(NA, c(n.keep, n))
samples.nu2 <- array(NA, c(n.keep, 1))
samples.tau2 <- array(NA, c(n.keep, 1))
samples.theta <- array(NA, c(n.keep, n))
samples.sigma2 <- array(NA, c(n.keep, 1))
samples.deviance <- array(NA, c(n.keep, n))
samples.fitted <- array(NA, c(n.keep, n))

## Metropolis quantities
sigma2.posterior.shape <- prior.sigma2[1] + 0.5*n
tau2.posterior.shape <- prior.tau2[1] + 0.5*(n-1)
nu2.posterior.shape <- prior.nu2[1] + 0.5*n
     
     
#### CAR quantities
    if(!is.matrix(W)) stop("W is not a matrix.", call.=FALSE)
    if(nrow(W)!= n) stop("W has the wrong number of rows.", call.=FALSE)
    if(ncol(W)!= n) stop("W has the wrong number of columns.", call.=FALSE)
    if(sum(is.na(W))>0) stop("W has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(W)) stop("W has non-numeric values.", call.=FALSE)
    if(min(W)<0) stop("W has negative elements.", call.=FALSE)
    if(sum(W!=t(W))>0) stop("W is not symmetric.", call.=FALSE)

### Create the duplet form
n.neighbours <- as.numeric(apply(W, 1, sum))
W.duplet <- c(NA, NA)
     for(i in 1:n)
     {
          for(j in 1:n)
          {
               if(W[i,j]==1)
               {
               W.duplet <- rbind(W.duplet, c(i,j))     
               }else{}
          }
     }
W.duplet <- W.duplet[-1, ]     
n.duplet <- nrow(W.duplet) 


## Create the list object
Wlist <- as.list(rep(NA,n))     
     for(i in 1:n)
     {
     Wlist[[i]] <- which(W[i, ]==1)     
     }


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
    data.mean.beta <- data.temp.beta %*% (Y - phi - theta - offset)
    U <- chol((prior.precision.beta + data.precision.beta / nu2))
    Uinv <- backsolve(U, diag(rep(1,p)))
    fc.mean.beta <- Uinv %*% (t(Uinv) %*% (prior.precision.beta %*% prior.mean.beta + (data.precision.beta / nu2) %*% data.mean.beta))
    beta <- fc.mean.beta + Uinv %*% rnorm(p)
    
    
    
    ##################
    ## Sample from nu2
    ##################
    fitted.current <-  as.numeric(X.standardised %*% beta) + phi + theta + offset
    nu2.posterior.scale <- prior.nu2[2] + 0.5 * sum((Y - fitted.current)^2)
    nu2 <- 1 / rgamma(1, nu2.posterior.shape, scale=(1/nu2.posterior.scale))    


    ####################
    ## Sample from phi
    ####################
    offset.phi <- (Y - as.numeric(X.standardised %*% beta) - theta - offset) / nu2    
    phi <- gaussiancarupdate(W_list=Wlist, nsites=n, phi=phi, nneighbours=n.neighbours, tau2=tau2, rho_num=1, rho_den=1, nu2=nu2, offset=offset.phi)
    phi <- phi - mean(phi)
    

    
    ####################
    ## Sample from theta
    ####################
    data.mean.theta <- Y - as.numeric(X.standardised %*% beta) - phi - offset    
    fc.variance.scalar <- 1/((1/nu2) + (1/sigma2))
    fc.mean <- fc.variance.scalar * (data.mean.theta / nu2)
    theta <- rnorm(n=n, mean=fc.mean, sd=sqrt(fc.variance.scalar))
    theta <- theta - mean(theta)



    ###################
    ## Sample from tau2
    ###################
    temp2 <- quadform(W_duplet1=W.duplet[ ,1], W_duplet2=W.duplet[ ,2], n_duplet=n.duplet,  nsites=n, phi=phi, nneighbours=n.neighbours, diagonal=1, offdiagonal=1)      
    tau2.posterior.scale <- temp2 + prior.tau2[2] 
    tau2 <- 1 / rgamma(1, tau2.posterior.shape, scale=(1/tau2.posterior.scale))
    
    
    
    #####################
    ## Sample from sigma2
    #####################
    sigma2.posterior.scale <- prior.sigma2[2] + 0.5*sum(theta^2)
    sigma2 <- 1 / rgamma(1, sigma2.posterior.shape, scale=(1/sigma2.posterior.scale))
    
    
    
    #########################
    ## Calculate the deviance
    #########################
    fitted <- as.numeric(X.standardised %*% beta) + phi + theta + offset
    deviance <- dnorm(Y, mean = fitted, sd = rep(sqrt(nu2),n))



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
        samples.theta[ele, ] <- theta
        samples.sigma2[ele, ] <- sigma2
        samples.deviance[ele, ] <- deviance
        samples.fitted[ele, ] <- fitted
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
accept.final <- rep(100, 6)
names(accept.final) <- c("beta", "phi", "nu2", "tau2", "theta", "sigma")
     
     
## Deviance information criterion (DIC)
median.beta <- apply(samples.beta, 2, median)
median.phi <- apply(samples.phi, 2, median)
median.theta <- apply(samples.theta, 2, median)
fitted.median <- X.standardised %*% median.beta + median.phi + median.theta + offset
nu2.median <- median(samples.nu2)
deviance.fitted <- -2 * sum(dnorm(Y, mean = fitted.median, sd = rep(sqrt(nu2.median),n), log = TRUE))
deviance.sum <- apply(-2 * log(samples.deviance), 1, sum)
p.d <- mean(deviance.sum) - deviance.fitted
DIC <- 2 * mean(deviance.sum) - deviance.fitted
like.fitted <- apply(samples.deviance, 2, mean)
DIC3 <- 2 * mean(deviance.sum)   + 2 * sum(log(like.fitted))     

     
#### Compute the Conditional Predictive Ordinate
CPO.temp <- 1 / samples.deviance
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
summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(100, p))
rownames(summary.beta) <- colnames(X)
colnames(summary.beta) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept")

summary.hyper <- array(NA, c(3,5))
summary.hyper[1, 1:3] <- quantile(samples.nu2, c(0.5, 0.025, 0.975))
summary.hyper[1, 4:5] <- c(n.keep, 100)
summary.hyper[2, 1:3] <- quantile(samples.tau2, c(0.5, 0.025, 0.975))
summary.hyper[2, 4:5] <- c(n.keep, 100)
summary.hyper[3, 1:3] <- quantile(samples.sigma2, c(0.5, 0.025, 0.975))
summary.hyper[3, 4:5] <- c(n.keep, 100)

summary.results <- rbind(summary.beta, summary.hyper)
rownames(summary.results)[(nrow(summary.results)-2):(nrow(summary.results))] <- c("nu2", "tau2", "sigma2")
summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
summary.results[ , 4:5] <- round(summary.results[ , 4:5], 1)


#### Create the Fitted values and residuals
fitted.values <- array(NA, c(n, 5))
colnames(fitted.values) <- c("Mean", "Sd", "Median", "2.5%", "97.5%")
fitted.values[ ,1] <- apply(samples.fitted, 2, mean)
fitted.values[ ,2] <- apply(samples.fitted, 2, sd)
fitted.values[ ,3:5] <- t(apply(samples.fitted, 2, quantile, c(0.5, 0.025, 0.975)))
fitted.values <- round(fitted.values, 4)

residuals <- array(NA, c(n, 5))
colnames(residuals) <- c("Mean", "Sd", "Median", "2.5%", "97.5%")
residuals.temp <- array(NA, c(nrow(samples.beta), n))
     for(i in 1:nrow(samples.beta))
     {
     residuals.temp[i, ] <- as.numeric(Y) - samples.fitted[i, ]
     }
residuals[ ,1] <- apply(residuals.temp, 2, mean)
residuals[ ,2] <- apply(residuals.temp, 2, sd)
residuals[ ,3:5] <- t(apply(residuals.temp, 2, quantile, c(0.5, 0.025, 0.975)))
residuals <- round(residuals, 4)


## Compile and return the results
modelfit <- c(DIC, p.d, DIC3, MPL)
names(modelfit) <- c("DIC", "p.d", "DIC3", "MPL")
model.string <- c("Likelihood model - Gaussian (identity link function)", "\nRandom effects model - BYM CAR\n")     
samples <- list(beta=samples.beta.orig, phi=mcmc(samples.phi), theta=mcmc(samples.theta), tau2=mcmc(samples.tau2), nu2=mcmc(samples.nu2), sigma2=mcmc(samples.sigma2))
results <- list(formula=formula, samples=samples, fitted.values=fitted.values, residuals=residuals, W.summary=W, modelfit=modelfit, summary.results=summary.results, model=model.string, accept=accept.final)
class(results) <- "carbayes"

     if(verbose)
     {
     b<-proc.time()
     cat(" finished in ", round(b[3]-a[3], 1), "seconds")
     }else
     {}
return(results)
}

