binomial.MVlerouxCAR <- function(formula, data=NULL, trials, W, burnin, n.sample, thin=1, n.chains=1, n.cores=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.Sigma.df=NULL, prior.Sigma.scale=NULL, rho=NULL, MALA=TRUE, verbose=TRUE)
{
##############################################
#### Format the arguments and check for errors
##############################################
#### Verbose
a <- common.verbose(verbose)
    
    
#### Frame object
frame.results <- common.frame(formula, data, "binomial")
K <- frame.results$n
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
J <- ncol(Y)
N.all <- K * J


#### Create a missing list
    if(n.miss>0)
    {
    miss.locator <- array(NA, c(n.miss, 2))
    colnames(miss.locator) <- c("row", "column")
    locations <- which(t(which.miss)==0)
    miss.locator[ ,1] <- ceiling(locations/J)
    miss.locator[ ,2] <- locations - (miss.locator[ ,1]-1) * J
    }else
    {}


#### Check on MALA argument
if(length(MALA)!=1) stop("MALA is not length 1.", call.=FALSE)
if(!is.logical(MALA)) stop("MALA is not logical.", call.=FALSE)  


#### Check and format the trials argument
    if(ncol(trials)!=J) stop("trials has the wrong number of columns.", call.=FALSE)
    if(nrow(trials)!=K) stop("trials has the wrong number of rows.", call.=FALSE)
    if(sum(is.na(trials))>0) stop("the numbers of trials has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(trials)) stop("the numbers of trials has non-numeric values.", call.=FALSE)
int.check <- N.all-sum(ceiling(trials)==floor(trials))
    if(int.check > 0) stop("the numbers of trials has non-integer values.", call.=FALSE)
    if(min(trials)<=0) stop("the numbers of trials has zero or negative values.", call.=FALSE)
failures <- trials - Y
    if(sum(Y>trials, na.rm=TRUE)>0) stop("the response variable has larger values that the numbers of trials.", call.=FALSE)


#### W matrix
    if(!is.matrix(W)) stop("W is not a matrix.", call.=FALSE)
    if(ceiling(N.all/K)!= floor(N.all/K)) stop("The number of data points divided by the number of rows in W is not a whole number.", call.=FALSE)


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



#### Priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(100000, p)
    if(is.null(prior.Sigma.df)) prior.Sigma.df <- 2
    if(is.null(prior.Sigma.scale)) prior.Sigma.scale <- rep(100000, J)
common.prior.beta.check(prior.mean.beta, prior.var.beta, p)
    if(!is.numeric(prior.Sigma.scale)) stop("prior.Sigma.scale has non-numeric values.", call.=FALSE)    
    if(sum(is.na(prior.Sigma.scale))!=0) stop("prior.Sigma.scale has missing values.", call.=FALSE)   


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
    results <- binomial.MVlerouxCARMCMC(Y=Y, failures=failures, trials=trials, offset=offset, X.standardised=X.standardised, W=W, rho=rho, fix.rho=fix.rho, K=K, p=p, J=J, N.all=N.all, which.miss=which.miss, n.miss=n.miss, miss.locator=miss.locator, burnin=burnin, n.sample=n.sample, thin=thin, n.beta.block=n.beta.block, list.block=list.block, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.Sigma.df=prior.Sigma.df, prior.Sigma.scale=prior.Sigma.scale, MALA=MALA, verbose=verbose, chain=1)
}else if(n.chains > 1 & ceiling(n.chains)==floor(n.chains) & n.cores==1)
{
    #### Multiple chains in  series
    results <- as.list(rep(NA, n.chains))
    for(i in 1:n.chains)
    {
        results[[i]] <- binomial.MVlerouxCARMCMC(Y=Y, failures=failures, trials=trials, offset=offset, X.standardised=X.standardised, W=W, rho=rho, fix.rho=fix.rho, K=K, p=p, J=J, N.all=N.all, which.miss=which.miss, n.miss=n.miss, miss.locator=miss.locator, burnin=burnin, n.sample=n.sample, thin=thin, n.beta.block=n.beta.block, list.block=list.block, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.Sigma.df=prior.Sigma.df, prior.Sigma.scale=prior.Sigma.scale, MALA=MALA, verbose=verbose, chain=i)
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
    results <- clusterCall(compclust, fun=binomial.MVlerouxCARMCMC, Y=Y, failures=failures, trials=trials, offset=offset, X.standardised=X.standardised, W=W, rho=rho, fix.rho=fix.rho, K=K, p=p, J=J, N.all=N.all, which.miss=which.miss, n.miss=n.miss, miss.locator=miss.locator, burnin=burnin, n.sample=n.sample, thin=thin, n.beta.block=n.beta.block, list.block=list.block, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.Sigma.df=prior.Sigma.df, prior.Sigma.scale=prior.Sigma.scale, MALA=MALA, verbose=verbose, chain="all")
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
accept.beta <- 100 * sum(results$accept.beta[1:J]) / sum(results$accept.beta[(J+1):(2*J)])
accept.phi <- 100 * results$accept[1] / results$accept[2]
accept.Sigma <- 100
    if(!fix.rho)
    {
    accept.rho <- 100 * results$accept[3] / results$accept[4]
    }else
    {
    accept.rho <- NA    
    }
accept.final <- c(accept.beta, accept.phi, accept.rho, accept.Sigma)
names(accept.final) <- c("beta", "phi", "rho", "Sigma")

    ## Compute the model fit criterion
    mean.beta <- matrix(apply(results$samples.beta, 2, mean), nrow=p, ncol=J, byrow=F)
    mean.phi <- matrix(apply(results$samples.phi, 2, mean), nrow=K, ncol=J, byrow=T)
    mean.logit <- X.standardised %*% mean.beta + mean.phi + offset
    mean.prob <- exp(mean.logit)  / (1 + exp(mean.logit))
    fitted.mean <- trials * mean.prob
    deviance.fitted <- -2 * sum(dbinom(x=as.numeric(t(Y)), size=as.numeric(t(trials)), prob=as.numeric(t(mean.prob)), log=TRUE), na.rm=TRUE)
    modelfit <- common.modelfit(results$samples.loglike, deviance.fitted)

    ## Create the Fitted values and residuals
    fitted.values <- matrix(apply(results$samples.fitted, 2, mean), nrow=K, ncol=J, byrow=T)
    response.residuals <- Y - fitted.values
    pearson.residuals <- response.residuals /sqrt(fitted.values * (1 - mean.prob))
    residuals <- list(response=response.residuals, pearson=pearson.residuals)

    ## Create MCMC objects and back transform the regression parameters
    samples.beta.orig <- results$samples.beta
        for(r in 1:J)
        {
        samples.beta.orig[ ,((r-1)*p+1):(r*p)] <- common.betatransform(results$samples.beta[ ,((r-1)*p+1):(r*p) ], X.indicator, X.mean, X.sd, p, FALSE)
        }
    samples <- list(beta=mcmc(samples.beta.orig), phi=mcmc(results$samples.phi), Sigma=results$samples.Sigma, rho=mcmc(results$samples.rho),  fitted=mcmc(results$samples.fitted), Y=mcmc(results$samples.Y))

    #### Create a summary object
    n.keep <- floor((n.sample - burnin)/thin)
    summary.beta <- t(rbind(apply(samples$beta, 2, mean), apply(samples$beta, 2, quantile, c(0.025, 0.975)))) 
    summary.beta <- cbind(summary.beta, rep(n.keep, J*p), rep(accept.beta,J*p), effectiveSize(samples$beta), geweke.diag(samples$beta)$z)
    col.name <- rep(NA, p*J)
        if(is.null(colnames(Y)))
        {
            for(r in 1:J)
            {
            col.name[((r-1)*p+1):(r*p)] <- paste("Variable ", r,  " - ", colnames(X), sep="")   
            }
        }else
        {
            for(r in 1:J)
            {
            col.name[((r-1)*p+1):(r*p)] <- paste(colnames(Y)[r],  " - ", colnames(X), sep="")   
            }
        }
    rownames(summary.beta) <- col.name
    colnames(summary.beta) <- c("Mean", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")

    summary.hyper <- array(NA, c((J+1) ,7))
    summary.hyper[1:J, 1] <- diag(apply(samples$Sigma, c(2,3), mean))
    summary.hyper[1:J, 2] <- diag(apply(samples$Sigma, c(2,3), quantile, c(0.025)))
    summary.hyper[1:J, 3] <- diag(apply(samples$Sigma, c(2,3), quantile, c(0.975)))
    summary.hyper[1:J, 4] <- n.keep
    summary.hyper[1:J, 5] <- accept.Sigma
    summary.hyper[1:J, 6] <- diag(apply(samples$Sigma, c(2,3), effectiveSize))
        for(r in 1:J)
        {
        summary.hyper[r, 7] <- geweke.diag(samples$Sigma[ ,r,r])$z    
        }
    
        if(!fix.rho)
        {
        summary.hyper[(J+1), 1:3] <- c(mean(samples$rho), quantile(samples$rho, c(0.025, 0.975)))
        summary.hyper[(J+1), 4:7] <- c(n.keep, accept.rho, effectiveSize(samples$rho), geweke.diag(samples$rho)$z)
        }else
        {
        summary.hyper[(J+1), 1:3] <- c(rho, rho, rho)
        summary.hyper[(J+1), 4:7] <- rep(NA, 4)
        }
    
    summary.results <- rbind(summary.beta, summary.hyper)
    rownames(summary.results)[((J*p)+1): nrow(summary.results)] <- c(paste(rep("Sigma",J), 1:J, 1:J, sep=""), "rho")
    summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
    summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
    }else
    {
    ## Compute the acceptance rates
    accept.temp <- lapply(results, function(l) l[["accept.beta"]])
    accept.temp2 <- do.call(what=rbind, args=accept.temp)
    accept.beta <- 100 * sum(accept.temp2[ ,1:J]) / sum(accept.temp2[ ,(J+1):(2*J)])
    accept.temp <- lapply(results, function(l) l[["accept"]])
    accept.temp2 <- do.call(what=rbind, args=accept.temp)
    accept.phi <- 100 * sum(accept.temp2[ ,1]) / sum(accept.temp2[ ,2])
    accept.Sigma <- 100
        if(!fix.rho)
        {
        accept.rho <- 100 * sum(accept.temp2[ ,3]) / sum(accept.temp2[ ,4])
        }else
        {
        accept.rho <- NA    
        }
    accept.final <- c(accept.beta, accept.phi, accept.rho, accept.Sigma)
    names(accept.final) <- c("beta", "phi", "rho", "Sigma")

    ## Extract the samples into separate lists
    samples.beta.list <- lapply(results, function(l) l[["samples.beta"]])
    samples.phi.list <- lapply(results, function(l) l[["samples.phi"]])
    samples.Sigma.list <- lapply(results, function(l) l[["samples.Sigma"]])
    samples.rho.list <- lapply(results, function(l) l[["samples.rho"]])    
    samples.loglike.list <- lapply(results, function(l) l[["samples.loglike"]])
    samples.fitted.list <- lapply(results, function(l) l[["samples.fitted"]])
    samples.Y.list <- lapply(results, function(l) l[["samples.Y"]])

    ## Convert the samples into separate matrix objects
    samples.beta.matrix <- do.call(what=rbind, args=samples.beta.list)
    samples.phi.matrix <- do.call(what=rbind, args=samples.phi.list)
    samples.rho.matrix <- do.call(what=rbind, args=samples.rho.list)
    samples.loglike.matrix <- do.call(what=rbind, args=samples.loglike.list)
    samples.fitted.matrix <- do.call(what=rbind, args=samples.fitted.list)
    
    ## Compute the model fit criteria
    mean.beta <- matrix(apply(samples.beta.matrix, 2, mean), nrow=p, ncol=J, byrow=F)
    mean.phi <- matrix(apply(samples.phi.matrix, 2, mean), nrow=K, ncol=J, byrow=T)
    mean.logit <- X.standardised %*% mean.beta + mean.phi + offset
    mean.prob <- exp(mean.logit)  / (1 + exp(mean.logit))
    fitted.mean <- trials * mean.prob
    deviance.fitted <- -2 * sum(dbinom(x=as.numeric(t(Y)), size=as.numeric(t(trials)), prob=as.numeric(t(mean.prob)), log=TRUE), na.rm=TRUE)
    modelfit <- common.modelfit(samples.loglike.matrix, deviance.fitted)

    ## Create the Fitted values and residuals
    fitted.values <- matrix(apply(samples.fitted.matrix, 2, mean), nrow=K, ncol=J, byrow=T)
    response.residuals <- Y - fitted.values
    pearson.residuals <- response.residuals /sqrt(fitted.values * (1 - mean.prob))
    residuals <- list(response=response.residuals, pearson=pearson.residuals)

    ## Backtransform the regression parameters
    samples.beta.list <- samples.beta.list
    for(j in 1:n.chains)
    {
        for(r in 1:J)
        {
            samples.beta.list[[j]][ ,((r-1)*p+1):(r*p)] <- common.betatransform(samples.beta.list[[j]][ ,((r-1)*p+1):(r*p)], X.indicator, X.mean, X.sd, p, FALSE)  
        }
    }
    samples.beta.matrix <- do.call(what=rbind, args=samples.beta.list)

    ## Create MCMC objects
    beta.temp <- samples.beta.list
    phi.temp <- samples.phi.list
    rho.temp <- samples.rho.list
    loglike.temp <- samples.loglike.list
    fitted.temp <- samples.fitted.list
    Y.temp <- samples.Y.list
        for(j in 1:n.chains)
        {
        beta.temp[[j]] <- mcmc(samples.beta.list[[j]])   
        phi.temp[[j]] <- mcmc(samples.phi.list[[j]])   
        rho.temp[[j]] <- mcmc(samples.rho.list[[j]])   
        loglike.temp[[j]] <- mcmc(samples.loglike.list[[j]])   
        fitted.temp[[j]] <- mcmc(samples.fitted.list[[j]])   
        Y.temp[[j]] <- mcmc(samples.Y.list[[j]])   
        }
    beta.mcmc <- as.mcmc.list(beta.temp)
    phi.mcmc <- as.mcmc.list(phi.temp)
    rho.mcmc <- as.mcmc.list(rho.temp)
    fitted.mcmc <- as.mcmc.list(fitted.temp)
    Y.mcmc <- as.mcmc.list(Y.temp)
    samples <- list(beta=beta.mcmc, phi=phi.mcmc, rho=rho.mcmc, Sigma=samples.Sigma.list, fitted=fitted.mcmc, Y=Y.mcmc)
    
    ## Create a summary object
    n.keep <- floor((n.sample - burnin)/thin)
    summary.beta <- t(rbind(apply(samples.beta.matrix, 2, mean), apply(samples.beta.matrix, 2, quantile, c(0.025, 0.975)))) 
    summary.beta <- cbind(summary.beta, rep(n.keep, J*p), rep(accept.beta,J*p), effectiveSize(beta.mcmc), gelman.diag(beta.mcmc)$psrf[ ,2])
    col.name <- rep(NA, p*J)
    
    if(is.null(colnames(Y)))
    {
        for(r in 1:J)
        {
            col.name[((r-1)*p+1):(r*p)] <- paste("Category ", r,  " - ", colnames(X), sep="")   
        }
    }else
    {
        for(r in 1:J)
        {
            col.name[((r-1)*p+1):(r*p)] <- paste(colnames(Y)[r],  " - ", colnames(X), sep="")   
        }
    }
    rownames(summary.beta) <- col.name
    colnames(summary.beta) <- c("Mean", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "PSRF (upper 95% CI)")
    
    summary.hyper <- array(NA, c((J+1) ,7))
    summary.hyper[1:J, 4] <- rep(n.keep, J)
    summary.hyper[1:J, 5] <- rep(accept.Sigma, J)
        for(r in 1:J)
        {
        test.vec <- samples.Sigma.list[[1]][ , r, r]
        test.list <- as.list(rep(NA, n.chains))
        test.list[[1]] <- mcmc(samples.Sigma.list[[1]][ , r, r])
            for(i in 2:n.chains)
            {
            test.vec <- c(test.vec, samples.Sigma.list[[i]][ , r, r])
            test.list[[i]] <- mcmc(samples.Sigma.list[[i]][ , r, r])
            }
        test.mcmc <- as.mcmc.list(test.list)
        summary.hyper[r,1]  <- mean(test.vec)   
        summary.hyper[r,2:3]  <- quantile(test.vec, c(0.025, 0.975))
        summary.hyper[r,6]  <- effectiveSize(test.mcmc)    
        summary.hyper[r,7]  <- gelman.diag(test.mcmc)$psrf[ ,2]    
        }
    if(!fix.rho)
    {
        summary.hyper[(J+1), 1:3] <- c(mean(samples.rho.matrix), quantile(samples.rho.matrix, c(0.025, 0.975)))
        summary.hyper[(J+1), 4:7] <- c(n.keep, accept.rho, effectiveSize(rho.mcmc), gelman.diag(rho.mcmc)$psrf[ ,2])
    }else
    {
        summary.hyper[(J+1), 1:3] <- c(rho, rho, rho)
        summary.hyper[(J+1), 4:7] <- rep(NA, 4)
    }
    summary.results <- rbind(summary.beta, summary.hyper)
    rownames(summary.results)[((J*p)+1): nrow(summary.results)] <- c(paste(rep("Sigma",J), 1:J, 1:J, sep=""), "rho")
    summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
    summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
}



###################################
#### Compile and return the results
###################################
model.string <- c("Likelihood model - Binomial (logit link function)", "\nRandom effects model - Leroux MCAR\n")
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