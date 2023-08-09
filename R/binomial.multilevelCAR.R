binomial.multilevelCAR <- function(formula, data=NULL, trials, W, ind.area, burnin, n.sample, thin=1, n.chains=1, n.cores=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.tau2=NULL, rho=NULL,  MALA=TRUE,  verbose=TRUE)
{
##############################################
#### Format the arguments and check for errors
##############################################
#### Verbose
a <- common.verbose(verbose)
    
    
#### Frame object
frame.results <- common.frame(formula, data, "binomial")
n <- frame.results$n
p <- frame.results$p
X <- frame.results$X
X.standardised <- frame.results$X.standardised
X.sd <- frame.results$X.sd
X.mean <- frame.results$X.mean
X.indicator <- frame.results$X.indicator 
offset <- frame.results$offset
Y <- frame.results$Y
which.miss <- frame.results$which.miss
n.miss <- frame.results$n.miss  
K <- length(unique(ind.area))
   
 
#### Check on MALA argument
if(length(MALA)!=1) stop("MALA is not length 1.", call.=FALSE)
if(!is.logical(MALA)) stop("MALA is not logical.", call.=FALSE)  


#### Check and format the trials argument
    if(sum(is.na(trials))>0) stop("the numbers of trials has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(trials)) stop("the numbers of trials has non-numeric values.", call.=FALSE)
int.check <- n-sum(ceiling(trials)==floor(trials))
    if(int.check > 0) stop("the numbers of trials has non-integer values.", call.=FALSE)
    if(min(trials)<=0) stop("the numbers of trials has zero or negative values.", call.=FALSE)
failures <- trials - Y
    if(sum(Y>trials, na.rm=TRUE)>0) stop("the response variable has larger values that the numbers of trials.", call.=FALSE)


#### rho
    if(is.null(rho))
    {
    rho <- runif(1)
    fix.rho <- FALSE   
    }else
    {
    fix.rho <- TRUE    
    }
    if(!is.numeric(rho) ) stop("rho is fixed but is not numeric.", call.=FALSE)  
    if(rho<0 ) stop("rho is outside the range [0, 1].", call.=FALSE)  
    if(rho>1 ) stop("rho is outside the range [0, 1].", call.=FALSE)  


#### Checks and formatting for ind.area
    if(!is.vector(ind.area)) stop("ind.area is not a vector.", call.=FALSE)
    if(sum(ceiling(ind.area)==floor(ind.area))!=n) stop("ind.area does not have all integer values.", call.=FALSE)    
    if(min(ind.area)!=1) stop("the minimum value in ind.area is not 1.", call.=FALSE)    
    if(max(ind.area)!=K) stop("the maximum value in ind.area is not equal to the number of spatial areal units.", call.=FALSE)    
    if(length(table(ind.area))!=K) stop("the number of unique areas in ind.area does not equal K.", call.=FALSE)    


#### Priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(100000, p)
    if(is.null(prior.tau2)) prior.tau2 <- c(1, 0.01)
common.prior.beta.check(prior.mean.beta, prior.var.beta, p)
common.prior.var.check(prior.tau2)    


#### Compute the blocking structure for beta     
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



########################
#### Run the MCMC chains
########################
if(n.chains==1)
{
    #### Only 1 chain
    results <- binomial.multilevelCARMCMC(Y=Y, failures=failures, trials=trials, offset=offset, X.standardised=X.standardised, W=W, ind.area=ind.area, rho=rho, fix.rho=fix.rho, n=n, K=K, p=p, which.miss=which.miss, n.miss=n.miss, burnin=burnin, n.sample=n.sample, thin=thin, MALA=MALA, n.beta.block=n.beta.block, list.block=list.block, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.tau2=prior.tau2, verbose=verbose, chain=1)
}else if(n.chains > 1 & ceiling(n.chains)==floor(n.chains) & n.cores==1)
{
    #### Multiple chains in  series
    results <- as.list(rep(NA, n.chains))
    for(i in 1:n.chains)
    {
        results[[i]] <- binomial.multilevelCARMCMC(Y=Y, failures=failures, trials=trials, offset=offset, X.standardised=X.standardised, W=W, ind.area=ind.area, rho=rho, fix.rho=fix.rho, n=n, K=K, p=p, which.miss=which.miss, n.miss=n.miss, burnin=burnin, n.sample=n.sample, thin=thin, MALA=MALA, n.beta.block=n.beta.block, list.block=list.block, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.tau2=prior.tau2, verbose=verbose, chain=i)
    }
}else if(n.chains > 1 & ceiling(n.chains)==floor(n.chains) & n.cores>1 & ceiling(n.cores)==floor(n.cores))   
{   
    #### Multiple chains in parallel
    results <- as.list(rep(NA, n.chains))
    if(verbose)
    {
        compclust <- makeCluster(n.cores, outfile="CARBayesprogress.txt")
        cat("The current progress of the model fitting algorithm has been output to CARBayesprogress.txt in the working directory")
    }else
    {
        compclust <- makeCluster(n.cores)
    }
    results <- clusterCall(compclust, fun=binomial.multilevelCARMCMC, Y=Y, failures=failures, trials=trials, offset=offset, X.standardised=X.standardised, W=W, ind.area=ind.area, rho=rho, fix.rho=fix.rho, n=n, K=K, p=p, which.miss=which.miss, n.miss=n.miss, burnin=burnin, n.sample=n.sample, thin=thin, MALA=MALA, n.beta.block=n.beta.block, list.block=list.block, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.tau2=prior.tau2, verbose=verbose, chain="all")
    stopCluster(compclust)
}else
{
    stop("n.chains or n.cores are not positive integers.", call.=FALSE)  
}


#### end timer
if(verbose)
{
    cat("\nSummarising results.\n")
}else
{}



###################################
#### Summarise and save the results 
###################################
if(n.chains==1)
{
    ## Compute the acceptance rates
    accept.beta <- 100 * results$accept[1] / results$accept[2]
    accept.phi <- 100 * results$accept[3] / results$accept[4]
    accept.tau2 <- 100
    if(!fix.rho)
    {
        accept.rho <- 100 * results$accept[5] / results$accept[6]
    }else
    {
        accept.rho <- NA    
    }
    accept.final <- c(accept.beta, accept.phi, accept.rho, accept.tau2)
    names(accept.final) <- c("beta", "phi", "rho", "tau2")
    
    ## Compute the model fit criterion
    mean.beta <- apply(results$samples.beta, 2, mean)
    mean.phi <- apply(results$samples.phi, 2, mean)
    mean.phi.extend <-  mean.phi[ind.area]
    mean.logit <- as.numeric(X.standardised %*% mean.beta) + mean.phi.extend + offset    
    mean.prob <- exp(mean.logit)  / (1 + exp(mean.logit))
    fitted.mean <- trials * mean.prob
    deviance.fitted <- -2 * sum(dbinom(x=Y, size=trials, prob=mean.prob, log=TRUE), na.rm=TRUE)
    modelfit <- common.modelfit(results$samples.loglike, deviance.fitted)
    
    ## Create the Fitted values and residuals
    fitted.values <- apply(results$samples.fitted, 2, mean)
    response.residuals <- as.numeric(Y) - fitted.values
    pearson.residuals <- response.residuals /sqrt(fitted.values * (1 - mean.prob))
    residuals <- data.frame(response=response.residuals, pearson=pearson.residuals)
    
    ## Create MCMC objects and back transform the regression parameters
    samples.beta.orig <- common.betatransform(results$samples.beta, X.indicator, X.mean, X.sd, p, FALSE)
    samples <- list(beta=mcmc(samples.beta.orig), phi=mcmc(results$samples.phi), rho=mcmc(results$samples.rho), tau2=mcmc(results$samples.tau2), fitted=mcmc(results$samples.fitted), Y=mcmc(results$samples.Y))
    
    #### Create a summary object
    n.keep <- floor((n.sample - burnin)/thin)
    summary.beta <- t(rbind(apply(samples$beta, 2, mean), apply(samples$beta, 2, quantile, c(0.025, 0.975)))) 
    summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta,p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
    rownames(summary.beta) <- colnames(X)
    colnames(summary.beta) <- c("Mean", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
    
    summary.hyper <- array(NA, c(2 ,7))
    summary.hyper[1, 1:3] <- c(mean(samples$tau2), quantile(samples$tau2, c(0.025, 0.975)))
    summary.hyper[1, 4:7] <- c(n.keep, accept.tau2, effectiveSize(samples$tau2), geweke.diag(samples$tau2)$z)
    if(!fix.rho)
    {
        summary.hyper[2, 1:3] <- c(mean(samples$rho), quantile(samples$rho, c(0.025, 0.975)))
        summary.hyper[2, 4:7] <- c(n.keep, accept.rho, effectiveSize(samples$rho), geweke.diag(samples$rho)$z)
    }else
    {
        summary.hyper[2, 1:3] <- c(rho, rho, rho)
        summary.hyper[2, 4:7] <- rep(NA, 4)
    }
    summary.results <- rbind(summary.beta, summary.hyper)
    rownames(summary.results)[(nrow(summary.results)-1):nrow(summary.results)] <- c("tau2", "rho")
    summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
    summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
}else
{
    ## Compute the acceptance rates
    accept.temp <- lapply(results, function(l) l[["accept"]])
    accept.temp2 <- do.call(what=rbind, args=accept.temp)
    accept.beta <- 100 * sum(accept.temp2[ ,1]) / sum(accept.temp2[ ,2])
    accept.phi <- 100 * sum(accept.temp2[ ,3]) / sum(accept.temp2[ ,4])
    accept.tau2 <- 100
    if(!fix.rho)
    {
        accept.rho <- 100 * sum(accept.temp2[ ,5]) / sum(accept.temp2[ ,6])
    }else
    {
        accept.rho <- NA    
    }
    accept.final <- c(accept.beta, accept.phi, accept.rho, accept.tau2)
    names(accept.final) <- c("beta", "phi", "rho", "tau2")
    
    ## Extract the samples into separate lists
    samples.beta.list <- lapply(results, function(l) l[["samples.beta"]])
    samples.phi.list <- lapply(results, function(l) l[["samples.phi"]])
    samples.rho.list <- lapply(results, function(l) l[["samples.rho"]])
    samples.tau2.list <- lapply(results, function(l) l[["samples.tau2"]])
    samples.loglike.list <- lapply(results, function(l) l[["samples.loglike"]])
    samples.fitted.list <- lapply(results, function(l) l[["samples.fitted"]])
    samples.Y.list <- lapply(results, function(l) l[["samples.Y"]])
    
    ## Convert the samples into separate matrix objects
    samples.beta.matrix <- do.call(what=rbind, args=samples.beta.list)
    samples.phi.matrix <- do.call(what=rbind, args=samples.phi.list)
    samples.rho.matrix <- do.call(what=rbind, args=samples.rho.list)
    samples.tau2.matrix <- do.call(what=rbind, args=samples.tau2.list)
    samples.loglike.matrix <- do.call(what=rbind, args=samples.loglike.list)
    samples.fitted.matrix <- do.call(what=rbind, args=samples.fitted.list)
    
    ## Compute the model fit criteria
    mean.beta <- apply(samples.beta.matrix, 2, mean)
    mean.phi <- apply(samples.phi.matrix, 2, mean)
    mean.phi.extend <-  mean.phi[ind.area]
    mean.logit <- as.numeric(X.standardised %*% mean.beta) + mean.phi.extend + offset    
    mean.prob <- exp(mean.logit)  / (1 + exp(mean.logit))
    fitted.mean <- trials * mean.prob
    deviance.fitted <- -2 * sum(dbinom(x=Y, size=trials, prob=mean.prob, log=TRUE), na.rm=TRUE)
    modelfit <- common.modelfit(samples.loglike.matrix, deviance.fitted)
    
    ## Create the Fitted values and residuals
    fitted.values <- apply(samples.fitted.matrix, 2, mean)
    response.residuals <- as.numeric(Y) - fitted.values
    pearson.residuals <- response.residuals /sqrt(fitted.values * (1 - mean.prob))
    residuals <- data.frame(response=response.residuals, pearson=pearson.residuals)
    
    ## Backtransform the regression parameters
    samples.beta.list <- samples.beta.list
    for(j in 1:n.chains)
    {
        samples.beta.list[[j]] <- common.betatransform(samples.beta.list[[j]], X.indicator, X.mean, X.sd, p, FALSE)  
    }
    samples.beta.matrix <- do.call(what=rbind, args=samples.beta.list)
    
    ## Create MCMC objects
    beta.temp <- samples.beta.list
    phi.temp <- samples.phi.list
    rho.temp <- samples.rho.list
    tau2.temp <- samples.tau2.list
    loglike.temp <- samples.loglike.list
    fitted.temp <- samples.fitted.list
    Y.temp <- samples.Y.list
    for(j in 1:n.chains)
    {
        beta.temp[[j]] <- mcmc(samples.beta.list[[j]])   
        phi.temp[[j]] <- mcmc(samples.phi.list[[j]])   
        rho.temp[[j]] <- mcmc(samples.rho.list[[j]])   
        tau2.temp[[j]] <- mcmc(samples.tau2.list[[j]])   
        loglike.temp[[j]] <- mcmc(samples.loglike.list[[j]])   
        fitted.temp[[j]] <- mcmc(samples.fitted.list[[j]])   
        Y.temp[[j]] <- mcmc(samples.Y.list[[j]])   
    }
    beta.mcmc <- as.mcmc.list(beta.temp)
    phi.mcmc <- as.mcmc.list(phi.temp)
    rho.mcmc <- as.mcmc.list(rho.temp)
    tau2.mcmc <- as.mcmc.list(tau2.temp)
    fitted.mcmc <- as.mcmc.list(fitted.temp)
    Y.mcmc <- as.mcmc.list(Y.temp)
    samples <- list(beta=beta.mcmc, phi=phi.mcmc, rho=rho.mcmc, tau2=tau2.mcmc, fitted=fitted.mcmc, Y=Y.mcmc)
    
    ## Create a summary object
    n.keep <- floor((n.sample - burnin)/thin)
    summary.beta <- t(rbind(apply(samples.beta.matrix, 2, mean), apply(samples.beta.matrix, 2, quantile, c(0.025, 0.975)))) 
    summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta,p), effectiveSize(beta.mcmc), gelman.diag(beta.mcmc)$psrf[ ,2])
    rownames(summary.beta) <- colnames(X)
    colnames(summary.beta) <- c("Mean", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "PSRF (upper 95% CI)")
    summary.hyper <- array(NA, c(2 ,7))
    summary.hyper[1, 1:3] <- c(mean(samples.tau2.matrix), quantile(samples.tau2.matrix, c(0.025, 0.975)))
    summary.hyper[1, 4:7] <- c(n.keep, accept.tau2, effectiveSize(tau2.mcmc),  gelman.diag(tau2.mcmc)$psrf[ ,2])
    if(!fix.rho)
    {
        summary.hyper[2, 1:3] <- c(mean(samples.rho.matrix), quantile(samples.rho.matrix, c(0.025, 0.975)))
        summary.hyper[2, 4:7] <- c(n.keep, accept.rho, effectiveSize(rho.mcmc),  gelman.diag(rho.mcmc)$psrf[ ,2])
    }else
    {
        summary.hyper[2, 1:3] <- c(rho, rho, rho)
        summary.hyper[2, 4:7] <- rep(NA, 4)
    }
    summary.results <- rbind(summary.beta, summary.hyper)
    rownames(summary.results)[(nrow(summary.results)-1):nrow(summary.results)] <- c("tau2", "rho")
    summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
    summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
}



###################################
#### Compile and return the results
###################################
model.string <- c("Likelihood model - Binomial (logit link function)", "\nRandom effects model - Multilevel Leroux CAR\n")
n.total <- floor((n.sample - burnin) / thin) * n.chains
mcmc.info <- c(n.total, n.sample, burnin, thin, n.chains)
names(mcmc.info) <- c("Total samples", "n.sample", "burnin", "thin", "n.chains")
results <- list(summary.results=summary.results, samples=samples, fitted.values=fitted.values, residuals=residuals, modelfit=modelfit, accept=accept.final, localised.structure=NULL,  formula=formula, model=model.string, mcmc.info=mcmc.info, X=X)
class(results) <- "CARBayes"
if(verbose)
{
    b<-proc.time()
    cat("Finished in ", round(b[3]-a[3], 1), "seconds.\n")
}else
{}
return(results)
}
