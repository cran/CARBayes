zip.bymCAR <- function(formula, formula.omega, data=NULL, W, burnin, n.sample, thin=1, n.chains=1, n.cores=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.tau2=NULL, prior.sigma2=NULL, prior.mean.delta=NULL, prior.var.delta=NULL, MALA=TRUE, verbose=TRUE)
{
##############################################
#### Format the arguments and check for errors
##############################################
#### Verbose
a <- common.verbose(verbose)
    

#### Frame object
frame.results <- common.frame(formula, data, "poisson")
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



#### Frame object for the omega model
## Create the matrix
frame.omega <- try(suppressWarnings(model.frame(formula.omega, data=data, na.action=na.pass)), silent=TRUE)
    if(class(frame.omega)[1]=="try-error") stop("the formula.omega inputted contains an error.", call.=FALSE)
V <- try(suppressWarnings(model.matrix(object=attr(frame.omega, "terms"), data=frame.omega)), silent=TRUE)
    if(class(V)[1]=="try-error") stop("the covariate matrix for the zero probabilities contains inappropriate values.", call.=FALSE)
    if(length(V)==0)
    {
    V <- matrix(rep(1, K), nrow=K, ncol=1, byrow=FALSE)    
    }else
    {}
    if(sum(is.na(V))>0) stop("the covariate matrix for the zero probabilities contains missing 'NA' values.", call.=FALSE)
    if(nrow(V)!=nrow(X)) stop("the two matrices of covariates don't have the same length.", call.=FALSE)
q <- ncol(V)

## Check for linearly related columns
cor.V <- suppressWarnings(cor(V))
diag(cor.V) <- 0
    if(max(cor.V, na.rm=TRUE)==1) stop("the covariate matrix for the zero probabilities has two exactly linearly related columns.", call.=FALSE)
    if(min(cor.V, na.rm=TRUE)==-1) stop("the covariate matrix for the zero probabilities has two exactly linearly related columns.", call.=FALSE)
    if(q>1)
    {
        if(sort(apply(V, 2, sd))[2]==0) stop("the covariate matrix for the zero probabilities has two intercept terms.", call.=FALSE)
    }else
    {}

## Standardise the matrix
V.standardised <- V
V.sd <- apply(V, 2, sd)
V.mean <- apply(V, 2, mean)
V.indicator <- rep(NA, q)       # To determine which parameter estimates to transform back

    for(j in 1:q)
    {
        if(length(table(V[ ,j]))>2)
        {
        V.indicator[j] <- 1
        V.standardised[ ,j] <- (V[ ,j] - mean(V[ ,j])) / sd(V[ ,j])
        }else if(length(table(V[ ,j]))==1)
        {
        V.indicator[j] <- 2
        }else
        {
        V.indicator[j] <- 0
        }
    }

## Check for an offset term
offset.omega <- try(model.offset(frame.omega), silent=TRUE)
    if(class(offset.omega)[1]=="try-error")   stop("the offset for the probability of being a zero is not numeric.", call.=FALSE)
    if(is.null(offset.omega))  offset.omega <- rep(0,K)
    if(sum(is.na(offset.omega))>0) stop("the offset for the probability of being a zero has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(offset.omega)) stop("the offset for the probability of being a zero variable has non-numeric values.", call.=FALSE)


#### Set up which elements are zero
which.zero <- which(Y==0)
n.zero <- length(which.zero)


#### Check on MALA argument
    if(length(MALA)!=1) stop("MALA is not length 1.", call.=FALSE)
    if(!is.logical(MALA)) stop("MALA is not logical.", call.=FALSE)  


#### Priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(100000, p)
    if(is.null(prior.tau2)) prior.tau2 <- c(1, 0.01)
    if(is.null(prior.sigma2)) prior.sigma2 <- c(1, 0.01)
    if(is.null(prior.mean.delta)) prior.mean.delta <- rep(0, q)
    if(is.null(prior.var.delta)) prior.var.delta <- rep(100000, q)
common.prior.beta.check(prior.mean.beta, prior.var.beta, p)
common.prior.var.check(prior.tau2)
common.prior.var.check(prior.sigma2)
common.prior.beta.check(prior.mean.delta, prior.var.delta, q)


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


## Compute the blocking structure for delta     
block.temp <- common.betablock(q)
delta.beg  <- block.temp[[1]]
delta.fin <- block.temp[[2]]
n.delta.block <- block.temp[[3]]
list.block.delta <- as.list(rep(NA, n.delta.block*2))
    for(r in 1:n.delta.block)
    {
    list.block.delta[[r]] <- delta.beg[r]:delta.fin[r]-1
    list.block.delta[[r+n.delta.block]] <- length(list.block.delta[[r]])
    }    



#### MCMC quantities - burnin, n.sample, thin
common.burnin.nsample.thin.check(burnin, n.sample, thin)



########################
#### Run the MCMC chains
########################
   if(n.chains==1)
   {
   #### Only 1 chain
   results <- zip.bymCARMCMC(Y=Y, offset=offset, offset.omega=offset.omega, X.standardised=X.standardised, V.standardised=V.standardised, W=W, K=K, p=p, q=q, which.miss=which.miss, n.miss=n.miss, which.zero=which.zero, n.zero=n.zero, burnin=burnin, n.sample=n.sample, thin=thin, MALA=MALA, n.beta.block=n.beta.block, list.block=list.block, n.delta.block=n.delta.block, list.block.delta=list.block.delta, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.mean.delta=prior.mean.delta, prior.var.delta=prior.var.delta, prior.tau2=prior.tau2, prior.sigma2=prior.sigma2, verbose=verbose, chain=1)
   }else if(n.chains > 1 & ceiling(n.chains)==floor(n.chains) & n.cores==1)
   {
   #### Multiple chains in  series
   results <- as.list(rep(NA, n.chains))
         for(i in 1:n.chains)
         {
         results[[i]] <- zip.bymCARMCMC(Y=Y, offset=offset, offset.omega=offset.omega, X.standardised=X.standardised, V.standardised=V.standardised, W=W, K=K, p=p, q=q, which.miss=which.miss, n.miss=n.miss, which.zero=which.zero, n.zero=n.zero, burnin=burnin, n.sample=n.sample, thin=thin, MALA=MALA, n.beta.block=n.beta.block, list.block=list.block,  n.delta.block=n.delta.block, list.block.delta=list.block.delta, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.mean.delta=prior.mean.delta, prior.var.delta=prior.var.delta, prior.tau2=prior.tau2, prior.sigma2=prior.sigma2, verbose=verbose, chain=i)
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
   results <- clusterCall(compclust, fun=zip.bymCARMCMC, Y=Y, offset=offset, offset.omega=offset.omega, X.standardised=X.standardised, V.standardised=V.standardised, W=W, K=K, p=p, q=q, which.miss=which.miss, n.miss=n.miss, which.zero=which.zero, n.zero=n.zero, burnin=burnin, n.sample=n.sample, thin=thin, MALA=MALA, n.beta.block=n.beta.block, list.block=list.block,  n.delta.block=n.delta.block, list.block.delta=list.block.delta, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.mean.delta=prior.mean.delta, prior.var.delta=prior.var.delta, prior.tau2=prior.tau2, prior.sigma2=prior.sigma2, verbose=verbose, chain="all")
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
    accept.delta <- 100 * results$accept[7] / results$accept[8]
    accept.phi <- 100 * results$accept[3] / results$accept[4]
    accept.theta <- 100 * results$accept[5] / results$accept[6]
    accept.tau2 <- 100
    accept.sigma2 <- 100
    accept.final <- c(accept.beta, accept.phi, accept.theta, accept.tau2, accept.sigma2, accept.delta)
    names(accept.final) <- c("beta", "phi", "theta", "tau2", "sigma2", "delta")

    ## Compute the model fit criterion
    mean.beta <- apply(results$samples.beta, 2, mean)
    mean.psi <- apply(results$samples.psi, 2, mean)
    mean.lp <- X.standardised %*% mean.beta  + mean.psi + offset
    mean.fitted <- exp(mean.lp)
    mean.Z <- round(apply(results$samples.Z,2,mean))
    mean.delta <- apply(results$samples.delta, 2, mean)
    mean.omega <- exp(V.standardised %*% mean.delta + offset.omega) / (1+exp(V.standardised %*% mean.delta + offset.omega))
    temp <- rep(0,K)
    temp[mean.Z==1] <- log(mean.omega[mean.Z==1])
    mean.deviance.all <- temp + (1-mean.Z) * (log(1-mean.omega) + dpois(x=as.numeric(Y), lambda=mean.fitted, log=T))
    deviance.fitted <- -2 * sum(mean.deviance.all, na.rm=TRUE)  
    modelfit <- common.modelfit(results$samples.loglike, deviance.fitted)

    ## Create the Fitted values and residuals
    fitted.values <- apply(results$samples.fitted, 2, mean)
    response.residuals <- as.numeric(Y) - fitted.values
    var.y <- fitted.values + (1 - mean.omega) * mean.omega * mean.fitted^2
    pearson.residuals <- response.residuals /sqrt(var.y)
    residuals <- data.frame(response=response.residuals, pearson=pearson.residuals)
    
    ## Create MCMC objects and back transform the regression parameters
    samples.beta.orig <- common.betatransform(results$samples.beta, X.indicator, X.mean, X.sd, p, FALSE)
    samples.delta.orig <- common.betatransform(results$samples.delta, V.indicator, V.mean, V.sd, q, FALSE)
    samples <- list(beta=mcmc(samples.beta.orig), psi=mcmc(results$samples.psi), tau2=mcmc(results$samples.tau2), sigma2=mcmc(results$samples.sigma2), delta=mcmc(samples.delta.orig), Z=mcmc(results$samples.Z), fitted=mcmc(results$samples.fitted), Y=mcmc(results$samples.Y))
    
    #### Create a summary object
    n.keep <- floor((n.sample - burnin)/thin)
    summary.beta <- t(rbind(apply(samples$beta, 2, mean), apply(samples$beta, 2, quantile, c(0.025, 0.975)))) 
    summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta,p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
    rownames(summary.beta) <- colnames(X)
    colnames(summary.beta) <- c("Mean", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
    summary.delta <- t(rbind(apply(samples$delta, 2, mean), apply(samples$delta, 2, quantile, c(0.025, 0.975)))) 
    summary.delta <- cbind(summary.delta, rep(n.keep, q), rep(accept.delta,q), effectiveSize(samples$delta), geweke.diag(samples$delta)$z)
    for(i in 1:q)
    {
        rownames(summary.delta)[i] <- paste("omega - ", colnames(V)[i])    
    }
    colnames(summary.delta) <- c("Mean", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
    
    summary.hyper <- array(NA, c(2 ,7))
    summary.hyper[1, 1:3] <- c(mean(samples$tau2), quantile(samples$tau2, c(0.025, 0.975)))
    summary.hyper[1, 4:7] <- c(n.keep, accept.tau2, effectiveSize(samples$tau2), geweke.diag(samples$tau2)$z)
    summary.hyper[2, 1:3] <- c(mean(samples$sigma2), quantile(samples$sigma2, c(0.025, 0.975)))
    summary.hyper[2, 4:7] <- c(n.keep, accept.sigma2, effectiveSize(samples$sigma2), geweke.diag(samples$sigma2)$z)
    summary.results <- rbind(summary.beta, summary.delta, summary.hyper)
    rownames(summary.results)[(nrow(summary.results)-1):nrow(summary.results)] <- c("tau2", "sigma2")
    summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
    summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
}else
{
    ## Compute the acceptance rates
    accept.temp <- lapply(results, function(l) l[["accept"]])
    accept.temp2 <- do.call(what=rbind, args=accept.temp)
    accept.beta <- 100 * sum(accept.temp2[ ,1]) / sum(accept.temp2[ ,2])
    accept.phi <- 100 * sum(accept.temp2[ ,3]) / sum(accept.temp2[ ,4])
    accept.delta <- 100 * sum(accept.temp2[ ,7]) / sum(accept.temp2[ ,8])
    accept.theta <- 100 * sum(accept.temp2[ ,5]) / sum(accept.temp2[ ,6])
    accept.tau2 <- 100
    accept.sigma2 <- 100
    accept.final <- c(accept.beta, accept.phi, accept.theta, accept.tau2, accept.sigma2, accept.delta)
    names(accept.final) <- c("beta", "phi", "theta", "tau2", "sigma2", "delta")
   
    ## Extract the samples into separate lists
    samples.beta.list <- lapply(results, function(l) l[["samples.beta"]])
    samples.psi.list <- lapply(results, function(l) l[["samples.psi"]])
    samples.sigma2.list <- lapply(results, function(l) l[["samples.sigma2"]])
    samples.tau2.list <- lapply(results, function(l) l[["samples.tau2"]])
    samples.delta.list <- lapply(results, function(l) l[["samples.delta"]])
    samples.Z.list <- lapply(results, function(l) l[["samples.Z"]])
    samples.loglike.list <- lapply(results, function(l) l[["samples.loglike"]])
    samples.fitted.list <- lapply(results, function(l) l[["samples.fitted"]])
    samples.Y.list <- lapply(results, function(l) l[["samples.Y"]])
    
    ## Convert the samples into separate matrix objects
    samples.beta.matrix <- do.call(what=rbind, args=samples.beta.list)
    samples.psi.matrix <- do.call(what=rbind, args=samples.psi.list)
    samples.sigma2.matrix <- do.call(what=rbind, args=samples.sigma2.list)
    samples.tau2.matrix <- do.call(what=rbind, args=samples.tau2.list)
    samples.delta.matrix <- do.call(what=rbind, args=samples.delta.list)
    samples.Z.matrix <- do.call(what=rbind, args=samples.Z.list)
    samples.loglike.matrix <- do.call(what=rbind, args=samples.loglike.list)
    samples.fitted.matrix <- do.call(what=rbind, args=samples.fitted.list)
    
    ## Compute the model fit criteria
    mean.beta <- apply(samples.beta.matrix, 2, mean)
    mean.psi <- apply(samples.psi.matrix, 2, mean)
    mean.lp <- X.standardised %*% mean.beta  + mean.psi + offset
    mean.fitted <- exp(mean.lp)
    mean.Z <- round(apply(samples.Z.matrix,2,mean))
    mean.delta <- apply(samples.delta.matrix, 2, mean)
    mean.omega <- exp(V.standardised %*% mean.delta + offset.omega) / (1+exp(V.standardised %*% mean.delta + offset.omega))
    temp <- rep(0,K)
    temp[mean.Z==1] <- log(mean.omega[mean.Z==1])
    mean.deviance.all <- temp + (1-mean.Z) * (log(1-mean.omega) + dpois(x=as.numeric(Y), lambda=mean.fitted, log=T))
    deviance.fitted <- -2 * sum(mean.deviance.all, na.rm=TRUE)  
    modelfit <- common.modelfit(samples.loglike.matrix, deviance.fitted)
    
    ## Create the Fitted values and residuals
    fitted.values <- apply(samples.fitted.matrix, 2, mean)
    response.residuals <- as.numeric(Y) - fitted.values
    var.y <- fitted.values + (1 - mean.omega) * mean.omega * mean.fitted^2
    pearson.residuals <- response.residuals /sqrt(var.y)
    residuals <- data.frame(response=response.residuals, pearson=pearson.residuals)
    
    ## Backtransform the regression parameters
    samples.beta.list <- samples.beta.list
    samples.delta.list <- samples.delta.list
    for(j in 1:n.chains)
    {
        samples.beta.list[[j]] <- common.betatransform(samples.beta.list[[j]], X.indicator, X.mean, X.sd, p, FALSE)  
        samples.delta.list[[j]] <- common.betatransform(samples.delta.list[[j]], V.indicator, V.mean, V.sd, q, FALSE)  
    }
    samples.beta.matrix <- do.call(what=rbind, args=samples.beta.list)
    samples.delta.matrix <- do.call(what=rbind, args=samples.delta.list)
    
    ## Create MCMC objects
    beta.temp <- samples.beta.list
    psi.temp <- samples.psi.list
    sigma2.temp <- samples.sigma2.list
    tau2.temp <- samples.tau2.list
    delta.temp <- samples.delta.list
    Z.temp <- samples.Z.list
    loglike.temp <- samples.loglike.list
    fitted.temp <- samples.fitted.list
    Y.temp <- samples.Y.list
    for(j in 1:n.chains)
    {
        beta.temp[[j]] <- mcmc(samples.beta.list[[j]])   
        psi.temp[[j]] <- mcmc(samples.psi.list[[j]])   
        sigma2.temp[[j]] <- mcmc(samples.sigma2.list[[j]])   
        tau2.temp[[j]] <- mcmc(samples.tau2.list[[j]]) 
        delta.temp[[j]] <- mcmc(samples.delta.list[[j]])   
        Z.temp[[j]] <- mcmc(samples.Z.list[[j]])  
        loglike.temp[[j]] <- mcmc(samples.loglike.list[[j]])   
        fitted.temp[[j]] <- mcmc(samples.fitted.list[[j]])   
        Y.temp[[j]] <- mcmc(samples.Y.list[[j]])   
    }
    beta.mcmc <- as.mcmc.list(beta.temp)
    psi.mcmc <- as.mcmc.list(psi.temp)
    sigma2.mcmc <- as.mcmc.list(sigma2.temp)
    tau2.mcmc <- as.mcmc.list(tau2.temp)
    delta.mcmc <- as.mcmc.list(delta.temp)
    Z.mcmc <- as.mcmc.list(Z.temp)
    fitted.mcmc <- as.mcmc.list(fitted.temp)
    Y.mcmc <- as.mcmc.list(Y.temp)
    samples <- list(beta=beta.mcmc, psi=psi.mcmc, sigma2=sigma2.mcmc, tau2=tau2.mcmc, delta=delta.mcmc, Z=Z.mcmc, fitted=fitted.mcmc, Y=Y.mcmc)
    
    ## Create a summary object
    n.keep <- floor((n.sample - burnin)/thin)
    summary.beta <- t(rbind(apply(samples.beta.matrix, 2, mean), apply(samples.beta.matrix, 2, quantile, c(0.025, 0.975)))) 
    summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta,p), effectiveSize(beta.mcmc), gelman.diag(beta.mcmc)$psrf[ ,2])
    rownames(summary.beta) <- colnames(X)
    colnames(summary.beta) <- c("Mean", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "PSRF (upper 95% CI)")
    summary.delta <- t(rbind(apply(samples.delta.matrix, 2, mean), apply(samples.delta.matrix, 2, quantile, c(0.025, 0.975)))) 
    summary.delta <- cbind(summary.delta, rep(n.keep, q), rep(accept.delta,q), effectiveSize(delta.mcmc), gelman.diag(delta.mcmc)$psrf[ ,2])
    for(i in 1:q)
    {
        rownames(summary.delta)[i] <- paste("omega - ", colnames(V)[i])    
    }
    colnames(summary.delta) <- c("Mean", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
    
    
    summary.hyper <- array(NA, c(2 ,7))
    summary.hyper[1, 1:3] <- c(mean(samples.tau2.matrix), quantile(samples.tau2.matrix, c(0.025, 0.975)))
    summary.hyper[1, 4:7] <- c(n.keep, accept.tau2, effectiveSize(tau2.mcmc),  gelman.diag(tau2.mcmc)$psrf[ ,2])
    summary.hyper[2, 1:3] <- c(mean(samples.sigma2.matrix), quantile(samples.sigma2.matrix, c(0.025, 0.975)))
    summary.hyper[2, 4:7] <- c(n.keep, accept.sigma2, effectiveSize(sigma2.mcmc),  gelman.diag(sigma2.mcmc)$psrf[ ,2])
    summary.results <- rbind(summary.beta, summary.delta, summary.hyper)
    rownames(summary.results)[(nrow(summary.results)-1):nrow(summary.results)] <- c("tau2", "sigma2")
    summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
    summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
}



###################################
#### Compile and return the results
###################################
model.string <- c("Likelihood model - Zero-Inflated Poisson (log link function)", "\nRandom effects model - BYM CAR\n")
n.total <- floor((n.sample - burnin) / thin) * n.chains
mcmc.info <- c(n.total, n.sample, burnin, thin, n.chains)
names(mcmc.info) <- c("Total samples", "n.sample", "burnin", "thin", "n.chains")
results <- list(summary.results=summary.results, samples=samples, fitted.values=fitted.values, residuals=residuals, modelfit=modelfit, accept=accept.final, localised.structure=NULL,  formula=c(formula, formula.omega), model=model.string, mcmc.info=mcmc.info, X=X)
class(results) <- "CARBayes"
if(verbose)
{
    b<-proc.time()
    cat("Finished in ", round(b[3]-a[3], 1), "seconds.\n")
}else
{}
return(results)
}