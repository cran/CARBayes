poisson.localisedCAR <- function(formula, data=NULL, G, W, burnin, n.sample, thin=1, n.chains=1, n.cores=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.tau2=NULL, prior.delta = NULL, MALA=TRUE, verbose=TRUE)
{
##############################################
#### Format the arguments and check for errors
##############################################
#### Verbose
a <- common.verbose(verbose)
    
    
#### Frame object
frame.results <- common.frame.localised(formula, data, "poisson", trials=NA)
K <- frame.results$n
p <- frame.results$p
X <- frame.results$X
X.standardised <- frame.results$X.standardised
X.sd <- frame.results$X.sd
X.mean <- frame.results$X.mean
X.indicator <- frame.results$X.indicator 
offset <- frame.results$offset
Y <- frame.results$Y
log.Y <- log(Y)
log.Y[Y==0] <- -0.1 


#### Check on MALA argument
    if(length(MALA)!=1) stop("MALA is not length 1.", call.=FALSE)
    if(!is.logical(MALA)) stop("MALA is not logical.", call.=FALSE)  


#### Format and check the number of clusters G     
    if(length(G)!=1) stop("G is the wrong length.", call.=FALSE)    
    if(!is.numeric(G)) stop("G is not numeric.", call.=FALSE)    
    if(G<=1) stop("G is less than 2.", call.=FALSE)    
    if(G!=round(G)) stop("G is not an integer.", call.=FALSE) 
    if(floor(G/2)==ceiling(G/2))
    {
    Gstar <- G/2    
    }else
    {
    Gstar <- (G+1)/2          
    }
  
     
#### Priors
    if(p>0)
    {
        if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
        if(is.null(prior.var.beta)) prior.var.beta <- rep(100000, p)
        common.prior.beta.check(prior.mean.beta, prior.var.beta, p)
    }else
    {
    prior.mean.beta <- NULL
    prior.var.beta <- NULL
    }

    if(is.null(prior.tau2)) prior.tau2 <- c(1, 0.01)
common.prior.var.check(prior.tau2)

    if(is.null(prior.delta)) prior.delta <- 10
    if(length(prior.delta)!=1) stop("the prior value for delta is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.delta)) stop("the prior value for delta is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.delta))!=0) stop("the prior value for delta has missing values.", call.=FALSE)    
    if(prior.delta<=0) stop("the prior value for delta is not positive.", call.=FALSE)    


#### Compute the blocking structure for beta     
    if(p>0)
    {
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
    }else
    {
    n.beta.block <- NULL
    list.block <- NULL   
    }


#### MCMC quantities - burnin, n.sample, thin
common.burnin.nsample.thin.check(burnin, n.sample, thin)  



########################
#### Run the MCMC chains
########################
   if(n.chains==1)
   {
   #### Only 1 chain
   results <- poisson.localisedCARMCMC(Y=Y, offset=offset, X.standardised=X.standardised, G=G, Gstar=Gstar, W=W, K=K, p=p, burnin=burnin, n.sample=n.sample, thin=thin, MALA=MALA, n.beta.block=n.beta.block, list.block=list.block, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.tau2=prior.tau2, prior.delta=prior.delta, verbose=verbose, chain=1)
   }else if(n.chains > 1 & ceiling(n.chains)==floor(n.chains) & n.cores==1)
   {
   #### Multiple chains in  series
   results <- as.list(rep(NA, n.chains))
         for(i in 1:n.chains)
         {
         results[[i]] <- poisson.localisedCARMCMC(Y=Y, offset=offset, X.standardised=X.standardised, G=G, Gstar=Gstar, W=W, K=K, p=p, burnin=burnin, n.sample=n.sample, thin=thin, MALA=MALA, n.beta.block=n.beta.block, list.block=list.block, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.tau2=prior.tau2, prior.delta=prior.delta, verbose=verbose, chain=i)
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
   results <- clusterCall(compclust, fun=poisson.localisedCARMCMC, Y=Y, offset=offset, X.standardised=X.standardised, G=G, Gstar=Gstar, W=W, K=K, p=p, burnin=burnin, n.sample=n.sample, thin=thin, MALA=MALA, n.beta.block=n.beta.block, list.block=list.block, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.tau2=prior.tau2, prior.delta=prior.delta, verbose=verbose, chain="all")
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
   accept.phi <- 100 * results$accept[1] / results$accept[2]
   accept.lambda <- 100 * results$accept[3] / results$accept[4]
   accept.delta <- 100 * results$accept[5] / results$accept[6]
   accept.tau2 <- 100
        if(p>0)
        {
        accept.beta <- 100 * results$accept[7] / results$accept[8]   
        accept.final <- c(accept.beta, accept.lambda, accept.delta, accept.phi, accept.tau2)
        names(accept.final) <- c("beta", "lambda", "delta", "phi", "tau2")   
        }else
        {
        accept.final <- c(accept.lambda,  accept.delta, accept.phi, accept.tau2)
        names(accept.final) <- c("lambda", "delta", "phi", "tau2")   
        }

   ## Compute the model fit criterion
   mean.phi <- apply(results$samples.phi, 2, mean)
   mean.Z <- round(apply(results$samples.Z,2,mean),0)
   mean.lambda <- apply(results$samples.lambda,2,mean)
    if(p>0)
    {
    mean.beta <- apply(results$samples.beta, 2, mean)
    regression.vec <- as.numeric(X.standardised %*% mean.beta)   
    }else
    {
    regression.vec <- rep(0,K)
    }
   fitted.mean <- exp(regression.vec  + mean.lambda[mean.Z] + mean.phi + offset)
   deviance.fitted <- -2 * sum(dpois(x=Y, lambda=fitted.mean, log=TRUE))
   modelfit <- common.modelfit(results$samples.loglike, deviance.fitted)

   ## Create the Fitted values and residuals
   fitted.values <- apply(results$samples.fitted, 2, mean)
   response.residuals <- as.numeric(Y) - fitted.values
   pearson.residuals <- response.residuals /sqrt(fitted.values)
   residuals <- data.frame(response=response.residuals, pearson=pearson.residuals)

   ## Create MCMC objects and back transform the regression parameters
       if(p>0)
       {    
       samples.beta.orig <- common.betatransform(results$samples.beta, X.indicator, X.mean, X.sd, p, FALSE)
       samples <- list(beta=mcmc(samples.beta.orig), phi=mcmc(results$samples.phi), lambda=mcmc(results$samples.lambda), Z=mcmc(results$samples.Z), tau2=mcmc(results$samples.tau2), delta=mcmc(results$samples.delta), fitted=mcmc(results$samples.fitted), Y=NA)          
       }else
       {
       samples <- list(phi=mcmc(results$samples.phi), lambda=mcmc(results$samples.lambda), Z=mcmc(results$samples.Z), tau2=mcmc(results$samples.tau2), delta=mcmc(results$samples.delta), fitted=mcmc(results$samples.fitted), Y=NA)          
       }     

   #### Create a summary object
   n.keep <- floor((n.sample - burnin)/thin)
   summary.hyper <- array(NA, c(2 ,7))
   summary.hyper[1, 1:3] <- c(mean(samples$tau2), quantile(samples$tau2, c(0.025, 0.975)))
   summary.hyper[1, 4:7] <- c(n.keep, accept.tau2, effectiveSize(samples$tau2), geweke.diag(samples$tau2)$z)
   summary.hyper[2, 1:3] <- c(mean(samples$delta), quantile(samples$delta, c(0.025, 0.975)))
   summary.hyper[2, 4:7] <- c(n.keep, accept.delta, effectiveSize(samples$delta), geweke.diag(samples$delta)$z)
   summary.lambda <- t(rbind(apply(samples$lambda, 2, mean), apply(samples$lambda, 2, quantile, c(0.025, 0.975)))) 
   summary.lambda <- cbind(summary.lambda, rep(n.keep, G), rep(accept.lambda, G), effectiveSize(samples$lambda), geweke.diag(samples$lambda)$z)
   Z.used <- as.numeric(names(table(samples$Z)))
   summary.lambda <- summary.lambda[Z.used, ]
        if(p>0)
        {
        summary.beta <- t(rbind(apply(samples$beta, 2, mean), apply(samples$beta, 2, quantile, c(0.025, 0.975)))) 
        summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta,p), effectiveSize(samples$beta), geweke.diag(samples$beta)$z)
        rownames(summary.beta) <- colnames(X)
        summary.results <- rbind(summary.beta, summary.lambda, summary.hyper)  
        row.names(summary.results)[(p+1):nrow(summary.results)] <- c(paste("lambda", Z.used, sep=""), "tau2", "delta")
        colnames(summary.results) <- c("Mean", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
        summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
        summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
        }else
        {
        summary.results <- rbind(summary.lambda, summary.hyper)
        row.names(summary.results) <- c(paste("lambda", Z.used, sep=""), "tau2", "delta")
        colnames(summary.results) <- c("Mean", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
        summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
        summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
        }
   }else
   {
   ## Compute the acceptance rates
   accept.temp <- lapply(results, function(l) l[["accept"]])
   accept.temp2 <- do.call(what=rbind, args=accept.temp)
   accept.phi <- 100 * sum(accept.temp2[ ,1]) / sum(accept.temp2[ ,2])
   accept.lambda <- 100 * sum(accept.temp2[ ,3]) / sum(accept.temp2[ ,4])
   accept.delta <- 100 * sum(accept.temp2[ ,5]) / sum(accept.temp2[ ,6])
   accept.tau2 <- 100
        if(p>0)
        {
        accept.beta <- 100 * sum(accept.temp2[ ,7]) / sum(accept.temp2[ ,8])   
        accept.final <- c(accept.beta, accept.lambda, accept.delta, accept.phi, accept.tau2)
        names(accept.final) <- c("beta", "lambda", "delta", "phi", "tau2")   
        }else
        {
        accept.final <- c(accept.lambda,  accept.delta, accept.phi, accept.tau2)
        names(accept.final) <- c("lambda", "delta", "phi", "tau2")   
        }

   ## Extract the samples into separate lists
   if(p>0) samples.beta.list <- lapply(results, function(l) l[["samples.beta"]])
   samples.phi.list <- lapply(results, function(l) l[["samples.phi"]])
   samples.delta.list <- lapply(results, function(l) l[["samples.delta"]])
   samples.lambda.list <- lapply(results, function(l) l[["samples.lambda"]])
   samples.Z.list <- lapply(results, function(l) l[["samples.Z"]])
   samples.tau2.list <- lapply(results, function(l) l[["samples.tau2"]])
   samples.loglike.list <- lapply(results, function(l) l[["samples.loglike"]])
   samples.fitted.list <- lapply(results, function(l) l[["samples.fitted"]])

   ## Convert the samples into separate matrix objects
   if(p>0) samples.beta.matrix <- do.call(what=rbind, args=samples.beta.list)
   samples.phi.matrix <- do.call(what=rbind, args=samples.phi.list)
   samples.delta.matrix <- do.call(what=rbind, args=samples.delta.list)
   samples.lambda.matrix <- do.call(what=rbind, args=samples.lambda.list)
   samples.Z.matrix <- do.call(what=rbind, args=samples.Z.list)
   samples.tau2.matrix <- do.call(what=rbind, args=samples.tau2.list)
   samples.loglike.matrix <- do.call(what=rbind, args=samples.loglike.list)
   samples.fitted.matrix <- do.call(what=rbind, args=samples.fitted.list)

   ## Compute the model fit criteria
   mean.phi <- apply(samples.phi.matrix, 2, mean)
   mean.Z <- round(apply(samples.Z.matrix,2,mean),0)
   mean.lambda <- apply(samples.lambda.matrix,2,mean)
    if(p>0)
    {
    mean.beta <- apply(samples.beta.matrix, 2, mean)
    regression.vec <- as.numeric(X.standardised %*% mean.beta)   
    }else
    {
    regression.vec <- rep(0,K)
    }
   fitted.mean <- exp(regression.vec  + mean.lambda[mean.Z] + mean.phi + offset)
   deviance.fitted <- -2 * sum(dpois(x=Y, lambda=fitted.mean, log=TRUE))
   modelfit <- common.modelfit(samples.loglike.matrix, deviance.fitted)

   ## Create the Fitted values and residuals
   fitted.values <- apply(samples.fitted.matrix, 2, mean)
   response.residuals <- as.numeric(Y) - fitted.values
   pearson.residuals <- response.residuals /sqrt(fitted.values)
   residuals <- data.frame(response=response.residuals, pearson=pearson.residuals)

   ## Backtransform the regression parameters
       if(p>0)
       {    
       samples.beta.list <- samples.beta.list
            for(j in 1:n.chains)
            {
            samples.beta.list[[j]] <- common.betatransform(samples.beta.list[[j]], X.indicator, X.mean, X.sd, p, FALSE)  
            }
       samples.beta.matrix <- do.call(what=rbind, args=samples.beta.list)
       }else
       {}
   
   ## Create MCMC objects
   if(p>0) beta.temp <- samples.beta.list
   phi.temp <- samples.phi.list
   delta.temp <- samples.delta.list
   lambda.temp <- samples.lambda.list
   Z.temp <- samples.Z.list
   tau2.temp <- samples.tau2.list
   loglike.temp <- samples.loglike.list
   fitted.temp <- samples.fitted.list
      for(j in 1:n.chains)
      {
      if(p>0) beta.temp[[j]] <- mcmc(samples.beta.list[[j]])   
      phi.temp[[j]] <- mcmc(samples.phi.list[[j]])   
      delta.temp[[j]] <- mcmc(samples.delta.list[[j]])   
      lambda.temp[[j]] <- mcmc(samples.lambda.list[[j]])   
      Z.temp[[j]] <- mcmc(samples.Z.list[[j]])   
      tau2.temp[[j]] <- mcmc(samples.tau2.list[[j]])   
      loglike.temp[[j]] <- mcmc(samples.loglike.list[[j]])   
      fitted.temp[[j]] <- mcmc(samples.fitted.list[[j]])   
      }
   if(p>0) beta.mcmc <- as.mcmc.list(beta.temp)
   phi.mcmc <- as.mcmc.list(phi.temp)
   delta.mcmc <- as.mcmc.list(delta.temp)
   Z.mcmc <- as.mcmc.list(Z.temp)
   lambda.mcmc <- as.mcmc.list(lambda.temp)
   tau2.mcmc <- as.mcmc.list(tau2.temp)
   fitted.mcmc <- as.mcmc.list(fitted.temp)
        if(p>0)
        {
        samples <- list(beta=beta.mcmc, phi=phi.mcmc, lambda=lambda.mcmc, Z=Z.mcmc, tau2=tau2.mcmc, delta=delta.mcmc, fitted=fitted.mcmc, Y=NA)    
        }else
        {
        samples <- list(phi=phi.mcmc, lambda=lambda.mcmc, Z=Z.mcmc, tau2=tau2.mcmc, delta=delta.mcmc, fitted=fitted.mcmc, Y=NA)    
        }

   ## Create a summary object
   n.keep <- floor((n.sample - burnin)/thin)
   summary.hyper <- array(NA, c(2 ,7))
   summary.hyper[1, 1:3] <- c(mean(samples.tau2.matrix), quantile(samples.tau2.matrix, c(0.025, 0.975)))
   summary.hyper[1, 4:7] <- c(n.keep, accept.tau2, effectiveSize(tau2.mcmc), gelman.diag(tau2.mcmc)$psrf[ ,2])
   summary.hyper[2, 1:3] <- c(mean(samples.delta.matrix), quantile(samples.delta.matrix, c(0.025, 0.975)))
   summary.hyper[2, 4:7] <- c(n.keep, accept.delta, effectiveSize(delta.mcmc), gelman.diag(delta.mcmc)$psrf[ ,2])
   summary.lambda <- t(rbind(apply(samples.lambda.matrix, 2, mean), apply(samples.lambda.matrix, 2, quantile, c(0.025, 0.975)))) 
   summary.lambda <- cbind(summary.lambda, rep(n.keep, G), rep(accept.lambda, G), effectiveSize(lambda.mcmc), gelman.diag(lambda.mcmc)$psrf[ ,2])
   Z.used <- as.numeric(names(table(samples.Z.matrix)))
   summary.lambda <- summary.lambda[Z.used, ]
        if(p>0)
        {
        summary.beta <- t(rbind(apply(samples.beta.matrix, 2, mean), apply(samples.beta.matrix, 2, quantile, c(0.025, 0.975)))) 
        summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta,p), effectiveSize(beta.mcmc), gelman.diag(beta.mcmc)$psrf[ ,2])
        rownames(summary.beta) <- colnames(X)
        summary.results <- rbind(summary.beta, summary.lambda, summary.hyper)  
        row.names(summary.results)[(p+1):nrow(summary.results)] <- c(paste("lambda", Z.used, sep=""), "tau2", "delta")
        colnames(summary.results) <- c("Mean", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "PSRF (upper 95% CI)")
        summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
        summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
        }else
        {
        summary.results <- rbind(summary.lambda, summary.hyper)
        row.names(summary.results) <- c(paste("lambda", Z.used, sep=""), "tau2", "delta")
        colnames(summary.results) <- c("Mean", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "PSRF (upper 95% CI)")
        summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
        summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
        }
}



###################################
#### Compile and return the results
###################################
model.string <- c("Likelihood model - Poisson (log link function)", "\nRandom effects  model - Localised CAR model\n")
n.total <- floor((n.sample - burnin) / thin) * n.chains
mcmc.info <- c(n.total, n.sample, burnin, thin, n.chains)
names(mcmc.info) <- c("Total samples", "n.sample", "burnin", "thin", "n.chains")
results <- list(summary.results=summary.results, samples=samples, fitted.values=fitted.values, residuals=residuals, modelfit=modelfit, accept=accept.final, localised.structure=mean.Z,  formula=formula, model=model.string, mcmc.info=mcmc.info, X=X)
class(results) <- "CARBayes"
    if(verbose)
    {
    b<-proc.time()
    cat("Finished in ", round(b[3]-a[3], 1), "seconds.\n")
    }else
    {}
return(results)
}
