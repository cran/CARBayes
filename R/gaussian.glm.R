gaussian.glm <- function(formula, data=NULL, burnin, n.sample, thin=1, n.chains=1, n.cores=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.nu2=NULL, verbose=TRUE)
{
##############################################
#### Format the arguments and check for errors
##############################################
#### Verbose
a <- common.verbose(verbose)
    
    
#### Frame object
frame.results <- common.frame(formula, data, "gaussian")
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

    
#### Priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(100000, p)
    if(is.null(prior.nu2)) prior.nu2 <- c(1, 0.01)
common.prior.beta.check(prior.mean.beta, prior.var.beta, p)
common.prior.var.check(prior.nu2)  

    
#### MCMC quantities - burnin, n.sample, thin
common.burnin.nsample.thin.check(burnin, n.sample, thin)      
    
    

########################
#### Run the MCMC chains
########################
   if(n.chains==1)
   {
   #### Only 1 chain
   results <- gaussian.glmMCMC(Y=Y, offset=offset, X.standardised=X.standardised, K=K, p=p, which.miss=which.miss, n.miss=n.miss, burnin=burnin, n.sample=n.sample, thin=thin, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.nu2=prior.nu2, verbose=verbose, chain=1)
   }else if(n.chains > 1 & ceiling(n.chains)==floor(n.chains) & n.cores==1)
   {
   #### Multiple chains in  series
   results <- as.list(rep(NA, n.chains))
         for(i in 1:n.chains)
         {
         results[[i]] <- gaussian.glmMCMC(Y=Y, offset=offset, X.standardised=X.standardised, K=K, p=p, which.miss=which.miss, n.miss=n.miss, burnin=burnin, n.sample=n.sample, thin=thin, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.nu2=prior.nu2, verbose=verbose, chain=i)
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
   results <- clusterCall(compclust, fun=gaussian.glmMCMC, Y=Y, offset=offset, X.standardised=X.standardised, K=K, p=p, which.miss=which.miss, n.miss=n.miss, burnin=burnin, n.sample=n.sample, thin=thin, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.nu2=prior.nu2, verbose=verbose, chain="all")
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
   accept.final <- rep(100, 2)
   names(accept.final) <- c("beta", "nu2")

   ## Compute the model fit criterion
   mean.beta <- apply(results$samples.beta, 2, mean)
   fitted.mean <- X.standardised %*% mean.beta + offset
   nu2.mean <- mean(results$samples.nu2)
   deviance.fitted <- -2 * sum(dnorm(Y, mean = fitted.mean, sd = rep(sqrt(nu2.mean),K), log = TRUE), na.rm=TRUE)
   modelfit <- common.modelfit(results$samples.loglike, deviance.fitted)      

   ## Create the Fitted values and residuals
   fitted.values <- apply(results$samples.fitted, 2, mean)
   response.residuals <- as.numeric(Y) - fitted.values
   pearson.residuals <- response.residuals /sqrt(nu2.mean)
   residuals <- data.frame(response=response.residuals, pearson=pearson.residuals)

   ## Create MCMC objects and back transform the regression parameters
   samples.beta.orig <- common.betatransform(results$samples.beta, X.indicator, X.mean, X.sd, p, FALSE)
   samples <- list(beta=mcmc(samples.beta.orig), nu2=mcmc(results$samples.nu2), fitted=mcmc(results$samples.fitted), Y=mcmc(results$samples.Y))
    
   #### Create a summary object
   n.keep <- floor((n.sample - burnin)/thin)
   summary.beta <- t(rbind(apply(samples$beta, 2, mean), apply(samples$beta, 2, quantile, c(0.025, 0.975)))) 
   summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(100,p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
   rownames(summary.beta) <- colnames(X)
   colnames(summary.beta) <- c("Mean", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
    summary.hyper <- array(NA, c(1 ,7))
    summary.hyper[1, 1:3] <- c(mean(samples$nu2), quantile(samples$nu2, c(0.025, 0.975)))
    summary.hyper[1, 4:7] <- c(n.keep, 100, effectiveSize(samples$nu2), geweke.diag(samples$nu2)$z)
    summary.results <- rbind(summary.beta, summary.hyper)
    rownames(summary.results)[nrow(summary.results)] <- c("nu2")
    summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
    summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
   }else
   {
   ## Compute the acceptance rates
   accept.final <- rep(100, 2)
   names(accept.final) <- c("beta", "nu2")

   ## Extract the samples into separate lists
   samples.beta.list <- lapply(results, function(l) l[["samples.beta"]])
   samples.nu2.list <- lapply(results, function(l) l[["samples.nu2"]])
   samples.loglike.list <- lapply(results, function(l) l[["samples.loglike"]])
   samples.fitted.list <- lapply(results, function(l) l[["samples.fitted"]])
   samples.Y.list <- lapply(results, function(l) l[["samples.Y"]])

   ## Convert the samples into separate matrix objects
   samples.beta.matrix <- do.call(what=rbind, args=samples.beta.list)
   samples.nu2.matrix <- do.call(what=rbind, args=samples.nu2.list)
   samples.loglike.matrix <- do.call(what=rbind, args=samples.loglike.list)
   samples.fitted.matrix <- do.call(what=rbind, args=samples.fitted.list)

   ## Compute the model fit criteria
   mean.beta <- apply(samples.beta.matrix, 2, mean)
   fitted.mean <- X.standardised %*% mean.beta + offset
   nu2.mean <- mean(samples.nu2.matrix)
   deviance.fitted <- -2 * sum(dnorm(Y, mean = fitted.mean, sd = rep(sqrt(nu2.mean),K), log = TRUE), na.rm=TRUE)
   modelfit <- common.modelfit(samples.loglike.matrix, deviance.fitted)      

   ## Create the Fitted values and residuals
   fitted.values <- apply(samples.fitted.matrix, 2, mean)
   response.residuals <- as.numeric(Y) - fitted.values
   pearson.residuals <- response.residuals /sqrt(nu2.mean)
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
   nu2.temp <- samples.nu2.list
   loglike.temp <- samples.loglike.list
   fitted.temp <- samples.fitted.list
   Y.temp <- samples.Y.list
      for(j in 1:n.chains)
      {
      beta.temp[[j]] <- mcmc(samples.beta.list[[j]])   
      nu2.temp[[j]] <- mcmc(samples.nu2.list[[j]])   
      loglike.temp[[j]] <- mcmc(samples.loglike.list[[j]])   
      fitted.temp[[j]] <- mcmc(samples.fitted.list[[j]])   
      Y.temp[[j]] <- mcmc(samples.Y.list[[j]])   
      }
   beta.mcmc <- as.mcmc.list(beta.temp)
   nu2.mcmc <- as.mcmc.list(nu2.temp)
   fitted.mcmc <- as.mcmc.list(fitted.temp)
   Y.mcmc <- as.mcmc.list(Y.temp)
   samples <- list(beta=beta.mcmc, nu2=nu2.mcmc, fitted=fitted.mcmc, Y=Y.mcmc)

   ## Create a summary object
   n.keep <- floor((n.sample - burnin)/thin)
   summary.beta <- t(rbind(apply(samples.beta.matrix, 2, mean), apply(samples.beta.matrix, 2, quantile, c(0.025, 0.975)))) 
   summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(100,p), effectiveSize(beta.mcmc), gelman.diag(beta.mcmc)$psrf[ ,2])
   rownames(summary.beta) <- colnames(X)
   colnames(summary.beta) <- c("Mean", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "PSRF (upper 95% CI)")
   summary.hyper <- array(NA, c(1 ,7))
   summary.hyper[1, 1:3] <- c(mean(samples.nu2.matrix), quantile(samples.nu2.matrix, c(0.025, 0.975)))
   summary.hyper[1, 4:7] <- c(n.keep, 100, effectiveSize(nu2.mcmc), gelman.diag(nu2.mcmc)$psrf[ ,2])
   summary.results <- rbind(summary.beta, summary.hyper)
   rownames(summary.results)[nrow(summary.results)] <- c("nu2")
   summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
   summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
   }



###################################
#### Compile and return the results
###################################
model.string <- c("Likelihood model - Gaussian (identity link function)", "\nRandom effects model - None\n")
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



