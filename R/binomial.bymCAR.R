binomial.bymCAR <- function(formula, data=NULL, trials, W, burnin, n.sample, thin=1, n.chains=1, n.cores=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.tau2=NULL, prior.sigma2=NULL, MALA=TRUE, verbose=TRUE)
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


    
#### Check on MALA argument
    if(length(MALA)!=1) stop("MALA is not length 1.", call.=FALSE)
    if(!is.logical(MALA)) stop("MALA is not logical.", call.=FALSE)  


#### Check and format the trials argument
    if(sum(is.na(trials))>0) stop("the numbers of trials has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(trials)) stop("the numbers of trials has non-numeric values.", call.=FALSE)
int.check <- K-sum(ceiling(trials)==floor(trials))
    if(int.check > 0) stop("the numbers of trials has non-integer values.", call.=FALSE)
    if(min(trials)<=0) stop("the numbers of trials has zero or negative values.", call.=FALSE)
failures <- trials - Y
    if(sum(Y>trials, na.rm=TRUE)>0) stop("the response variable has larger values that the numbers of trials.", call.=FALSE)


#### Priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(100000, p)
    if(is.null(prior.tau2)) prior.tau2 <- c(1, 0.01)
    if(is.null(prior.sigma2)) prior.sigma2 <- c(1, 0.01)
common.prior.beta.check(prior.mean.beta, prior.var.beta, p)
common.prior.var.check(prior.tau2)
common.prior.var.check(prior.sigma2)


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



########################
#### Run the MCMC chains
########################
   if(n.chains==1)
   {
   #### Only 1 chain
   results <- binomial.bymCARMCMC(Y=Y, failures=failures, trials=trials, offset=offset, X.standardised=X.standardised, W=W, K=K, p=p, which.miss=which.miss, n.miss=n.miss, burnin=burnin, n.sample=n.sample, thin=thin, MALA=MALA, n.beta.block=n.beta.block, list.block=list.block, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.tau2=prior.tau2, prior.sigma2=prior.sigma2, verbose=verbose, chain=1)
   }else if(n.chains > 1 & ceiling(n.chains)==floor(n.chains) & n.cores==1)
   {
   #### Multiple chains in  series
   results <- as.list(rep(NA, n.chains))
         for(i in 1:n.chains)
         {
         results[[i]] <- binomial.bymCARMCMC(Y=Y, failures=failures, trials=trials, offset=offset, X.standardised=X.standardised, W=W, K=K, p=p, which.miss=which.miss, n.miss=n.miss, burnin=burnin, n.sample=n.sample, thin=thin, MALA=MALA, n.beta.block=n.beta.block, list.block=list.block, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.tau2=prior.tau2, prior.sigma2=prior.sigma2, verbose=verbose, chain=i)
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
   results <- clusterCall(compclust, fun=binomial.bymCARMCMC, Y=Y, failures=failures, trials=trials, offset=offset, X.standardised=X.standardised, W=W, K=K, p=p, which.miss=which.miss, n.miss=n.miss, burnin=burnin, n.sample=n.sample, thin=thin, MALA=MALA, n.beta.block=n.beta.block, list.block=list.block, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.tau2=prior.tau2, prior.sigma2=prior.sigma2, verbose=verbose, chain="all")
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
   accept.theta <- 100 * results$accept[5] / results$accept[6]
   accept.tau2 <- 100
   accept.sigma2 <- 100
   accept.final <- c(accept.beta, accept.phi, accept.theta, accept.tau2, accept.sigma2)
   names(accept.final) <- c("beta", "phi", "theta", "tau2", "sigma2")

   ## Compute the model fit criterion
   mean.beta <- apply(results$samples.beta, 2, mean)
   mean.psi <- apply(results$samples.psi, 2, mean)
   mean.logit <- as.numeric(X.standardised %*% mean.beta) + mean.psi + offset    
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
   samples <- list(beta=mcmc(samples.beta.orig), psi=mcmc(results$samples.psi), tau2=mcmc(results$samples.tau2), sigma2=mcmc(results$samples.sigma2), fitted=mcmc(results$samples.fitted), Y=mcmc(results$samples.Y))
    
   #### Create a summary object
   n.keep <- floor((n.sample - burnin)/thin)
   summary.beta <- t(rbind(apply(samples$beta, 2, mean), apply(samples$beta, 2, quantile, c(0.025, 0.975)))) 
   summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta,p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
   rownames(summary.beta) <- colnames(X)
   colnames(summary.beta) <- c("Mean", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
   
   summary.hyper <- array(NA, c(2 ,7))
   summary.hyper[1, 1:3] <- c(mean(samples$tau2), quantile(samples$tau2, c(0.025, 0.975)))
   summary.hyper[1, 4:7] <- c(n.keep, accept.tau2, effectiveSize(samples$tau2), geweke.diag(samples$tau2)$z)
   summary.hyper[2, 1:3] <- c(mean(samples$sigma2), quantile(samples$sigma2, c(0.025, 0.975)))
   summary.hyper[2, 4:7] <- c(n.keep, accept.sigma2, effectiveSize(samples$sigma2), geweke.diag(samples$sigma2)$z)
   summary.results <- rbind(summary.beta, summary.hyper)
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
   accept.theta <- 100 * sum(accept.temp2[ ,5]) / sum(accept.temp2[ ,6])
   accept.tau2 <- 100
   accept.sigma2 <- 100   
   accept.final <- c(accept.beta, accept.phi, accept.theta, accept.tau2, accept.sigma2)
   names(accept.final) <- c("beta", "phi", "theta", "tau2", "sigma2")

   ## Extract the samples into separate lists
   samples.beta.list <- lapply(results, function(l) l[["samples.beta"]])
   samples.psi.list <- lapply(results, function(l) l[["samples.psi"]])
   samples.sigma2.list <- lapply(results, function(l) l[["samples.sigma2"]])
   samples.tau2.list <- lapply(results, function(l) l[["samples.tau2"]])
   samples.loglike.list <- lapply(results, function(l) l[["samples.loglike"]])
   samples.fitted.list <- lapply(results, function(l) l[["samples.fitted"]])
   samples.Y.list <- lapply(results, function(l) l[["samples.Y"]])

   ## Convert the samples into separate matrix objects
   samples.beta.matrix <- do.call(what=rbind, args=samples.beta.list)
   samples.psi.matrix <- do.call(what=rbind, args=samples.psi.list)
   samples.sigma2.matrix <- do.call(what=rbind, args=samples.sigma2.list)
   samples.tau2.matrix <- do.call(what=rbind, args=samples.tau2.list)
   samples.loglike.matrix <- do.call(what=rbind, args=samples.loglike.list)
   samples.fitted.matrix <- do.call(what=rbind, args=samples.fitted.list)

   ## Compute the model fit criteria
   mean.beta <- apply(samples.beta.matrix, 2, mean)
   mean.psi <- apply(samples.psi.matrix, 2, mean)
   mean.logit <- as.numeric(X.standardised %*% mean.beta) + mean.psi + offset    
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
   psi.temp <- samples.psi.list
   sigma2.temp <- samples.sigma2.list
   tau2.temp <- samples.tau2.list
   loglike.temp <- samples.loglike.list
   fitted.temp <- samples.fitted.list
   Y.temp <- samples.Y.list
      for(j in 1:n.chains)
      {
      beta.temp[[j]] <- mcmc(samples.beta.list[[j]])   
      psi.temp[[j]] <- mcmc(samples.psi.list[[j]])   
      sigma2.temp[[j]] <- mcmc(samples.sigma2.list[[j]])   
      tau2.temp[[j]] <- mcmc(samples.tau2.list[[j]])   
      loglike.temp[[j]] <- mcmc(samples.loglike.list[[j]])   
      fitted.temp[[j]] <- mcmc(samples.fitted.list[[j]])   
      Y.temp[[j]] <- mcmc(samples.Y.list[[j]])   
      }
   beta.mcmc <- as.mcmc.list(beta.temp)
   psi.mcmc <- as.mcmc.list(psi.temp)
   sigma2.mcmc <- as.mcmc.list(sigma2.temp)
   tau2.mcmc <- as.mcmc.list(tau2.temp)
   fitted.mcmc <- as.mcmc.list(fitted.temp)
   Y.mcmc <- as.mcmc.list(Y.temp)
   samples <- list(beta=beta.mcmc, psi=psi.mcmc, tau2=tau2.mcmc, sigma2=sigma2.mcmc, fitted=fitted.mcmc, Y=Y.mcmc)

   ## Create a summary object
   n.keep <- floor((n.sample - burnin)/thin)
   summary.beta <- t(rbind(apply(samples.beta.matrix, 2, mean), apply(samples.beta.matrix, 2, quantile, c(0.025, 0.975)))) 
   summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta,p), effectiveSize(beta.mcmc), gelman.diag(beta.mcmc)$psrf[ ,2])
   rownames(summary.beta) <- colnames(X)
   colnames(summary.beta) <- c("Mean", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "PSRF (upper 95% CI)")
   summary.hyper <- array(NA, c(2 ,7))
   summary.hyper[1, 1:3] <- c(mean(samples.tau2.matrix), quantile(samples.tau2.matrix, c(0.025, 0.975)))
   summary.hyper[1, 4:7] <- c(n.keep, accept.tau2, effectiveSize(samples.tau2.matrix),  gelman.diag(tau2.mcmc)$psrf[ ,2])
   summary.hyper[2, 1:3] <- c(mean(samples.sigma2.matrix), quantile(samples.sigma2.matrix, c(0.025, 0.975)))
   summary.hyper[2, 4:7] <- c(n.keep, accept.sigma2, effectiveSize(samples.sigma2.matrix),  gelman.diag(sigma2.mcmc)$psrf[ ,2])
   summary.results <- rbind(summary.beta, summary.hyper)
   rownames(summary.results)[(nrow(summary.results)-1):nrow(summary.results)] <- c("tau2", "sigma2")
   summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
   summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
   }



###################################
#### Compile and return the results
###################################
model.string <- c("Likelihood model - Binomial (logit link function)", "\nRandom effects model - BYM CAR\n")
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