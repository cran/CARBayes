poisson.glm <- function(formula, data=NULL,  burnin, n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, verbose=TRUE)
{
##############################################
#### Format the arguments and check for errors
##############################################
#### Verbose
a <- common.verbose(verbose)
    
    
#### Frame object
frame.results <- common.frame(formula, data, "poisson")
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
    

#### Priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(100000, p)
common.prior.beta.check(prior.mean.beta, prior.var.beta, p)

    
## Compute the blocking structure for beta     
block.temp <- common.betablock(p)
beta.beg  <- block.temp[[1]]
beta.fin <- block.temp[[2]]
n.beta.block <- block.temp[[3]]
list.block <- as.list(rep(NA, n.beta.block*2))
    for(r in 1:n.beta.block)
    {
    list.block[[r]] <- beta.beg[r]:beta.fin[r]-1
    list.block[[r+n.beta.block]] <- length(list.block[[r]])
    }
    
    
#### MCMC quantities - burnin, n.sample, thin
common.burnin.nsample.thin.check(burnin, n.sample, thin)      
    
    
    
#############################
#### Initial parameter values
#############################
#### Initial parameter values
mod.glm <- glm(Y~X.standardised-1, offset=offset, family="quasipoisson")
beta.mean <- mod.glm$coefficients
beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)
 
    
    
###############################    
#### Set up the MCMC quantities    
###############################
#### Matrices to store samples   
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, p))
samples.deviance <- array(NA, c(n.keep, 1))
samples.like <- array(NA, c(n.keep, n))
samples.fitted <- array(NA, c(n.keep, n))
    if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))
    
    
#### Metropolis quantities
accept.all <- rep(0,2)
accept <- accept.all
proposal.sd.beta <- 0.01

    

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
    offset.temp <- offset
        if(p>2)
        {
        temp <- poissonbetaupdateMALA(X.standardised, n, p, beta, offset.temp, Y.miss, prior.mean.beta, prior.var.beta, which.miss, n.beta.block, proposal.sd.beta, list.block)
        }else
        {
        temp <- poissonbetaupdateRW(X.standardised, n, p, beta, offset.temp, Y.miss, prior.mean.beta, prior.var.beta, which.miss, proposal.sd.beta)
        }
    beta <- temp[[1]]
    accept[1] <- accept[1] + temp[[2]]
    accept[2] <- accept[2] + n.beta.block  
        
   
         
    #########################
    ## Calculate the deviance
    #########################
    lp <- as.numeric(X.standardised %*% beta) + offset
    fitted <- exp(lp)
    deviance.all <- dpois(x=as.numeric(Y), lambda=fitted, log=TRUE)
    like <- exp(deviance.all)
    deviance <- -2 * sum(deviance.all, na.rm=TRUE)  
        
        
        
    ###################
    ## Save the results
    ###################
        if(j > burnin & (j-burnin)%%thin==0)
        {
        ele <- (j - burnin) / thin
        samples.beta[ele, ] <- beta
        samples.deviance[ele, ] <- deviance
        samples.like[ele, ] <- like
        samples.fitted[ele, ] <- fitted
            if(n.miss>0) samples.Y[ele, ] <- rpois(n=n.miss, lambda=fitted[which.miss==0])
        }else
        {}
        
        
    
    ########################################
    ## Self tune the acceptance probabilties
    ########################################
    k <- j/100
        if(ceiling(k)==floor(k))
        {
        #### Update the proposal sds
            if(p>2)
            {
            proposal.sd.beta <- common.accceptrates1(accept[1:2], proposal.sd.beta, 40, 50)
            }else
            {
            proposal.sd.beta <- common.accceptrates1(accept[1:2], proposal.sd.beta, 30, 40)    
            }
            accept.all <- accept.all + accept
            accept <- rep(0,2)
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
    
    
#### end timer
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
accept.beta <- 100 * accept.all[1] / accept.all[2]
accept.final <- c(accept.beta)
names(accept.final) <- c("beta")
    
    
#### Deviance information criterion (DIC)
median.beta <- apply(samples.beta, 2, median)
fitted.median <- exp(X.standardised %*% median.beta  + offset)
deviance.fitted <- -2 * sum(dpois(x=Y, lambda=fitted.median, log=TRUE), na.rm=TRUE)
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
    CPO[j] <- 1/median((1 / dpois(x=Y[j], lambda=samples.fitted[ ,j])))    
    }
LMPL <- sum(log(CPO), na.rm=TRUE)  
    
    
#### Compute the % deviance explained
fit.null <- glm(Y~1, offset=offset, family="poisson")$fitted.values
deviance.null <- -2 * sum(dpois(x=Y[which(!is.na(Y))], lambda=fit.null, log=TRUE), na.rm=TRUE)
percent_dev_explained <- 100 * (deviance.null - deviance.fitted) / deviance.null


#### transform the parameters back to the origianl covariate scale.
samples.beta.orig <- common.betatransform(samples.beta, X.indicator, X.mean, X.sd, p, FALSE)
    
    
#### Create a summary object
samples.beta.orig <- mcmc(samples.beta.orig)
summary.beta <- t(apply(samples.beta.orig, 2, quantile, c(0.5, 0.025, 0.975))) 
summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta,p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
rownames(summary.beta) <- colnames(X)
colnames(summary.beta) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
summary.results <- summary.beta
summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)

    
#### Create the Fitted values and residuals
fitted.values <- apply(samples.fitted, 2, median)
response.residuals <- as.numeric(Y) - fitted.values
pearson.residuals <- response.residuals /sqrt(fitted.values)
deviance.residuals <- sign(response.residuals) * sqrt(2 * (Y * log(Y/fitted.values) + fitted.values - Y))
residuals <- data.frame(response=response.residuals, pearson=pearson.residuals, deviance=deviance.residuals)
    
    
#### Compile and return the results
loglike <- (-0.5 * deviance.fitted)
modelfit <- c(DIC, p.d, WAIC, p.w, LMPL, loglike, percent_dev_explained)
names(modelfit) <- c("DIC", "p.d", "WAIC", "p.w", "LMPL", "loglikelihood", "Percentage deviance explained")
model.string <- c("Likelihood model - Poisson (log link function)", "\nRandom effects model - None\n")
    if(n.miss==0) samples.Y = NA

samples <- list(beta=samples.beta.orig, fitted=mcmc(samples.fitted), Y=mcmc(samples.Y))
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
