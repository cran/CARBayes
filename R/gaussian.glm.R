gaussian.glm <- function(formula, data=NULL, burnin, n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.nu2=NULL, verbose=TRUE)
{
##############################################
#### Format the arguments and check for errors
##############################################
#### Verbose
a <- common.verbose(verbose)
    
    
#### Frame object
frame.results <- common.frame(formula, data, "gaussian")
n <- frame.results$n
p <- frame.results$p
X <- frame.results$X
X.standardised <- frame.results$X.standardised
X.sd <- frame.results$X.sd
X.mean <- frame.results$X.mean
X.indicator <- frame.results$X.indicator 
offset <- frame.results$offset
Y <- frame.results$Y
Y.miss <- frame.results$Y.miss
which.miss <- frame.results$which.miss
n.miss <- frame.results$n.miss  
Y.short <- Y[which.miss==1]
X.short <- X.standardised[which.miss==1, ]
offset.short <- offset[which.miss==1]
    
    
#### Priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(100000, p)
    if(is.null(prior.nu2)) prior.nu2 <- c(1, 0.01)
common.prior.beta.check(prior.mean.beta, prior.var.beta, p)
common.prior.var.check(prior.nu2)  

    
#### MCMC quantities - burnin, n.sample, thin
common.burnin.nsample.thin.check(burnin, n.sample, thin)      
    
    
    
#############################
#### Initial parameter values
#############################
#### Initial parameter values
mod.glm <- lm(Y~X.standardised-1, offset=offset)
beta.mean <- mod.glm$coefficients
beta.sd <- sqrt(diag(summary(mod.glm)$cov.unscaled)) * summary(mod.glm)$sigma
beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)
res.temp <- Y - X.standardised %*% beta.mean - offset
nu2 <- var(as.numeric(res.temp), na.rm=TRUE)

    
    
###############################    
#### Set up the MCMC quantities    
###############################
#### Matrices to store samples
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, p))
samples.nu2 <- array(NA, c(n.keep, 1))
samples.deviance <- array(NA, c(n.keep, 1))
samples.like <- array(NA, c(n.keep, n))
samples.fitted <- array(NA, c(n.keep, n))
    if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))
    
    
#### Metropolis quantities
nu2.posterior.shape <- prior.nu2[1] + 0.5*(n-n.miss)

    
#### Beta update quantities
data.precision.beta <- t(X.short) %*% X.short
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
#### Start timer
    if(verbose)
    {
    cat("Generating", n.keep, "post burnin and thinned (if requested) samples\n", sep = " ")
    progressBar <- txtProgressBar(style = 3)
    percentage.points<-round((1:100/100)*n.sample)
    }else
    {
    percentage.points<-round((1:100/100)*n.sample)     
    }
    
    
#### Create the MCMC samples    
    for(j in 1:n.sample)
    {
    ####################
    ## Sample from beta
    ####################
    fc.precision <- prior.precision.beta + data.precision.beta / nu2
    fc.var <- solve(fc.precision)
    beta.offset <- as.numeric(Y.short - offset.short)
    beta.offset2 <- t(X.short) %*% beta.offset / nu2 + prior.precision.beta %*% prior.mean.beta
    fc.mean <- fc.var %*% beta.offset2
    chol.var <- t(chol(fc.var))
    beta <- fc.mean + chol.var %*% rnorm(p)        
        
        
        
    ##################
    ## Sample from nu2
    ##################
    fitted.current <-  as.numeric(X.short %*% beta) + offset.short
    nu2.posterior.scale <- prior.nu2[2] + 0.5 * sum((Y.short - fitted.current)^2)
    nu2 <- 1 / rgamma(1, nu2.posterior.shape, scale=(1/nu2.posterior.scale))    

        
    #########################
    ## Calculate the deviance
    #########################
    fitted <- as.numeric(X.standardised %*% beta) + offset
    deviance.all <- dnorm(Y, mean = fitted, sd = rep(sqrt(nu2),n), log=TRUE)
    like <- exp(deviance.all)
    deviance <- -2 * sum(deviance.all, na.rm=TRUE)  
        
        
        
    ###################
    ## Save the results
    ###################
        if(j > burnin & (j-burnin)%%thin==0)
        {
        ele <- (j - burnin) / thin
        samples.beta[ele, ] <- beta
        samples.nu2[ele, ] <- nu2
        samples.deviance[ele, ] <- deviance
        samples.like[ele, ] <- like
        samples.fitted[ele, ] <- fitted
            if(n.miss>0) samples.Y[ele, ] <- rnorm(n=n.miss, mean=fitted[which.miss==0], sd=sqrt(nu2))
        }else
        {}

        
        
    ################################       
    ## print progress to the console
    ################################
        if(j %in% percentage.points & verbose)
        {
        setTxtProgressBar(progressBar, j/n.sample)
        }
    }
    
    
##### end timer
    if(verbose)
    {
    cat("\nSummarising results")
    close(progressBar)
    }else
    {}
    

    
###################################
#### Summarise and save the results 
###################################
#### Compute the acceptance rates
accept.final <- rep(100, 2)
names(accept.final) <- c("beta", "nu2")
    
    
#### Deviance information criterion (DIC)
median.beta <- apply(samples.beta, 2, median)
fitted.median <- X.standardised %*% median.beta + offset
nu2.median <- median(samples.nu2)
deviance.fitted <- -2 * sum(dnorm(Y, mean = fitted.median, sd = rep(sqrt(nu2.median),n), log = TRUE), na.rm=TRUE)
p.d <- median(samples.deviance) - deviance.fitted
DIC <- 2 * median(samples.deviance) - deviance.fitted    
    
    
#### Watanabe-Akaike Information Criterion (WAIC)
LPPD <- sum(log(apply(samples.like,2,mean)), na.rm=TRUE)
p.w <- sum(apply(log(samples.like),2,var), na.rm=TRUE)
WAIC <- -2 * (LPPD - p.w)
    
    
#### Compute the Conditional Predictive Ordinate  
CPO <- rep(NA, n)
    for(j in 1:n)
    {
    CPO[j] <- 1/median((1 / dnorm(Y[j], mean=samples.fitted[ ,j], sd=sqrt(samples.nu2))))    
    }
LMPL <- sum(log(CPO), na.rm=TRUE)       
    
    
#### Compute the % deviance explained
deviance.null <- sum((Y - mean(Y, na.rm=T))^2, na.rm=T)
deviance.fit <- sum((Y - fitted.median)^2, na.rm=T)
percent_dev_explained <- 100 * (deviance.null - deviance.fit) / deviance.null


#### transform the parameters back to the origianl covariate scale.
samples.beta.orig <- common.betatransform(samples.beta, X.indicator, X.mean, X.sd, p, FALSE)
    
    
#### Create a summary object
samples.beta.orig <- mcmc(samples.beta.orig)
summary.beta <- t(apply(samples.beta.orig, 2, quantile, c(0.5, 0.025, 0.975))) 
summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(100,p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
rownames(summary.beta) <- colnames(X)
colnames(summary.beta) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
    
summary.hyper <- array(NA, c(1 ,7))
summary.hyper[1, 1:3] <- quantile(samples.nu2, c(0.5, 0.025, 0.975))
summary.hyper[1, 4:7] <- c(n.keep, 100, effectiveSize(samples.nu2), geweke.diag(samples.nu2)$z)
summary.results <- rbind(summary.beta, summary.hyper)
rownames(summary.results)[nrow(summary.results)] <- c("nu2")
summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)

    
#### Create the Fitted values and residuals
fitted.values <- apply(samples.fitted, 2, median)
response.residuals <- as.numeric(Y) - fitted.values
pearson.residuals <- response.residuals /sqrt(nu2.median)
deviance.residuals <- sign(response.residuals) * sqrt((Y-fitted.values)^2/nu2.median)
residuals <- data.frame(response=response.residuals, pearson=pearson.residuals, deviance=deviance.residuals)
    
    
#### Compile and return the results
loglike <- (-0.5 * deviance.fitted)
modelfit <- c(DIC, p.d, WAIC, p.w, LMPL, loglike, percent_dev_explained)
names(modelfit) <- c("DIC", "p.d", "WAIC", "p.w", "LMPL", "loglikelihood", "Percentage deviance explained")
model.string <- c("Likelihood model - Gaussian (identity link function)", "\nRandom effects model - None\n")
    if(n.miss==0) samples.Y = NA
 
samples <- list(beta=samples.beta.orig, nu2=mcmc(samples.nu2), fitted=mcmc(samples.fitted), Y=mcmc(samples.Y))
    results <- list(summary.results=summary.results, samples=samples, fitted.values=fitted.values, residuals=residuals, modelfit=modelfit, accept=accept.final, localised.structure=NULL,  formula=formula, model=model.string, X=X)
    class(results) <- "CARBayes"
    
    #### Finish by stating the time taken    
    if(verbose)
    {
        b<-proc.time()
        cat(" finished in ", round(b[3]-a[3], 1), "seconds")
    }else
    {}
    return(results)
}

