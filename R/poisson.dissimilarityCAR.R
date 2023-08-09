poisson.dissimilarityCAR <- function(formula, data=NULL, W, Z, W.binary=TRUE, burnin, n.sample, thin=1, n.chains=1, n.cores=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.tau2=NULL, MALA=TRUE, verbose=TRUE)
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


#### Check on MALA argument
if(length(MALA)!=1) stop("MALA is not length 1.", call.=FALSE)
if(!is.logical(MALA)) stop("MALA is not logical.", call.=FALSE)  


#### Dissimilarity metric matrix
    if(!is.list(Z)) stop("Z is not a list object.", call.=FALSE)
    if(sum(is.na(as.numeric(lapply(Z, sum, na.rm=FALSE))))>0) stop("Z contains missing 'NA' values.", call.=FALSE)
q <- length(Z)
	if(sum(as.numeric(lapply(Z,nrow))==K) <q) stop("Z contains matrices of the wrong size.", call.=FALSE)
	if(sum(as.numeric(lapply(Z,ncol))==K) <q) stop("Z contains matrices of the wrong size.", call.=FALSE)
	if(min(as.numeric(lapply(Z,min)))<0) stop("Z contains negative values.", call.=FALSE)

    if(!is.logical(W.binary)) stop("W.binary is not TRUE or FALSE.", call.=FALSE)
    if(length(W.binary)!=1) stop("W.binary has the wrong length.", call.=FALSE)

    if(W.binary)
    {
    alpha.max <- rep(NA,q)
    alpha.threshold <- rep(NA,q)
        for(k in 1:q)
        {
        Z.crit <- quantile(as.numeric(Z[[k]])[as.numeric(Z[[k]])!=0], 0.5)
        alpha.max[k] <- -log(0.5) / Z.crit
        alpha.threshold[k] <- -log(0.5) / max(Z[[k]])
        }        
    }else
    {
    alpha.max <- rep(50, q)    
    }


#### Priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(100000, p)
    if(is.null(prior.tau2)) prior.tau2 <- c(1, 0.01)
common.prior.beta.check(prior.mean.beta, prior.var.beta, p)
common.prior.var.check(prior.tau2)


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
   results <- poisson.dissimilarityCARMCMC(Y=Y, offset=offset, X.standardised=X.standardised, Z=Z, W.binary=W.binary, W=W, K=K, p=p, q=q, which.miss=which.miss, n.miss=n.miss, burnin=burnin, n.sample=n.sample, thin=thin, MALA=MALA, n.beta.block=n.beta.block, list.block=list.block, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.tau2=prior.tau2, alpha.max=alpha.max, verbose=verbose, chain=1)
   }else if(n.chains > 1 & ceiling(n.chains)==floor(n.chains) & n.cores==1)
   {
   #### Multiple chains in  series
   results <- as.list(rep(NA, n.chains))
         for(i in 1:n.chains)
         {
         results[[i]] <- poisson.dissimilarityCARMCMC(Y=Y, offset=offset, X.standardised=X.standardised, Z=Z, W.binary=W.binary, W=W, K=K, p=p, q=q, which.miss=which.miss, n.miss=n.miss, burnin=burnin, n.sample=n.sample, thin=thin, MALA=MALA, n.beta.block=n.beta.block, list.block=list.block, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.tau2=prior.tau2, alpha.max=alpha.max, verbose=verbose, chain=i)
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
   results <- clusterCall(compclust, fun=poisson.dissimilarityCARMCMC, Y=Y, offset=offset, X.standardised=X.standardised, Z=Z, W.binary=W.binary, W=W, K=K, p=p, q=q, which.miss=which.miss, n.miss=n.miss, burnin=burnin, n.sample=n.sample, thin=thin, MALA=MALA, n.beta.block=n.beta.block, list.block=list.block, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.tau2=prior.tau2, alpha.max=alpha.max, verbose=verbose, chain="all")
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
   accept.alpha <- 100 * results$accept[5] / results$accept[6]
   accept.final <- c(accept.beta, accept.phi, accept.tau2, accept.alpha)
   names(accept.final) <- c("beta", "phi", "tau2", "alpha")          

   ## Compute the model fit criterion
   mean.beta <- apply(results$samples.beta, 2, mean)
   mean.phi <- apply(results$samples.phi, 2, mean)
   fitted.mean <- exp(X.standardised %*% mean.beta + mean.phi + offset)
   deviance.fitted <- -2 * sum(dpois(x=Y, lambda=fitted.mean, log=TRUE), na.rm=TRUE)
   modelfit <- common.modelfit(results$samples.loglike, deviance.fitted)

   ## Create the Fitted values and residuals
   fitted.values <- apply(results$samples.fitted, 2, mean)
   response.residuals <- as.numeric(Y) - fitted.values
   pearson.residuals <- response.residuals /sqrt(fitted.values)
   residuals <- data.frame(response=response.residuals, pearson=pearson.residuals)

   ## Create MCMC objects and back transform the regression parameters
   samples.beta.orig <- common.betatransform(results$samples.beta, X.indicator, X.mean, X.sd, p, FALSE)
   samples <- list(beta=mcmc(samples.beta.orig), phi=mcmc(results$samples.phi), alpha=mcmc(results$samples.alpha), tau2=mcmc(results$samples.tau2), fitted=mcmc(results$samples.fitted), Y=mcmc(results$samples.Y))
    
   #### Create a summary object
   n.keep <- floor((n.sample - burnin)/thin)
   summary.beta <- t(rbind(apply(samples$beta, 2, mean), apply(samples$beta, 2, quantile, c(0.025, 0.975)))) 
   summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta,p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
   rownames(summary.beta) <- colnames(X)
   colnames(summary.beta) <- c("Mean", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
   summary.alpha <- t(rbind(apply(samples$alpha, 2, mean), apply(samples$alpha, 2, quantile, c(0.025, 0.975)))) 
   summary.alpha <- cbind(summary.alpha, rep(n.keep, q), rep(accept.alpha,q), effectiveSize(samples$alpha), geweke.diag(samples$alpha)$z)
	    if(!is.null(names(Z)))
 	    {
 	    rownames(summary.alpha) <- names(Z)
	    }else
	    {
	    names.Z <- rep(NA,q)
		    for(j in 1:q)
		    {
		    names.Z[j] <- paste("Z[[",j, "]]", sep="")
		    }	
	    rownames(summary.alpha) <- names.Z	
	    }

   summary.hyper <- array(NA, c(1 ,7))
   summary.hyper[1, 1:3] <- c(mean(samples$tau2), quantile(samples$tau2, c(0.025, 0.975)))
   summary.hyper[1, 4:7] <- c(n.keep, accept.tau2, effectiveSize(samples$tau2), geweke.diag(samples$tau2)$z)
   summary.results <- rbind(summary.beta, summary.hyper, summary.alpha)
        if(W.binary)
        {
        alpha.min <- c(rep(NA, (p+1)), alpha.threshold)
        summary.results <- cbind(summary.results, alpha.min)    
        }else
        {}
   rownames(summary.results)[(p+1)] <- c("tau2")
   summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
   summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
        if(W.binary) summary.results[ , 8] <- round(summary.results[ , 8], 4)
   
   #### Create the posterior medians for the neighbourhood matrix W
   W.posterior <- array(NA, c(K,K))
    if(W.binary)
    {
    W.border.prob <- array(NA, c(K,K))    
    }else
    {
    W.border.prob <- NA
    }

	for(i in 1:K)
	{
		for(j in 1:K)
		{
			if(W[i,j]==1)
			{
			z.temp <- NA
				for(k in 1:q)
				{
				z.temp <- c(z.temp, Z[[k]][i,j])
				}	
			z.temp <- z.temp[-1]
			w.temp <- exp(-samples$alpha %*% z.temp)
			    if(W.binary)
			    {
			    w.posterior <- as.numeric(w.temp>=0.5)
			    W.posterior[i,j] <- ceiling(median(w.posterior))
			    W.border.prob[i,j] <- (1 - sum(w.posterior) / length(w.posterior))    
			    }else
			    {
		        W.posterior[i,j] <- median(w.temp)
			    }
			}else
			{
			}	
		}	
	}
   }else
   {
   ## Compute the acceptance rates
   accept.temp <- lapply(results, function(l) l[["accept"]])
   accept.temp2 <- do.call(what=rbind, args=accept.temp)
   accept.beta <- 100 * sum(accept.temp2[ ,1]) / sum(accept.temp2[ ,2])
   accept.phi <- 100 * sum(accept.temp2[ ,3]) / sum(accept.temp2[ ,4])
   accept.tau2 <- 100
   accept.alpha <- 100 * sum(accept.temp2[ ,5]) / sum(accept.temp2[ ,6])
   accept.final <- c(accept.beta, accept.phi, accept.tau2, accept.alpha)
   names(accept.final) <- c("beta", "phi", "tau2", "alpha")
   
   ## Extract the samples into separate lists
   samples.beta.list <- lapply(results, function(l) l[["samples.beta"]])
   samples.phi.list <- lapply(results, function(l) l[["samples.phi"]])
   samples.alpha.list <- lapply(results, function(l) l[["samples.alpha"]])
   samples.tau2.list <- lapply(results, function(l) l[["samples.tau2"]])
   samples.loglike.list <- lapply(results, function(l) l[["samples.loglike"]])
   samples.fitted.list <- lapply(results, function(l) l[["samples.fitted"]])
   samples.Y.list <- lapply(results, function(l) l[["samples.Y"]])

   ## Convert the samples into separate matrix objects
   samples.beta.matrix <- do.call(what=rbind, args=samples.beta.list)
   samples.phi.matrix <- do.call(what=rbind, args=samples.phi.list)
   samples.alpha.matrix <- do.call(what=rbind, args=samples.alpha.list)
   samples.tau2.matrix <- do.call(what=rbind, args=samples.tau2.list)
   samples.loglike.matrix <- do.call(what=rbind, args=samples.loglike.list)
   samples.fitted.matrix <- do.call(what=rbind, args=samples.fitted.list)

   ## Compute the model fit criteria
   mean.beta <- apply(samples.beta.matrix, 2, mean)
   mean.phi <- apply(samples.phi.matrix, 2, mean)
   fitted.mean <- exp(X.standardised %*% mean.beta + mean.phi + offset)   
   deviance.fitted <- -2 * sum(dpois(x=Y, lambda=fitted.mean, log=TRUE), na.rm=TRUE)
   modelfit <- common.modelfit(samples.loglike.matrix, deviance.fitted)

   ## Create the Fitted values and residuals
   fitted.values <- apply(samples.fitted.matrix, 2, mean)
   response.residuals <- as.numeric(Y) - fitted.values
   pearson.residuals <- response.residuals /sqrt(fitted.values)
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
   alpha.temp <- samples.alpha.list
   tau2.temp <- samples.tau2.list
   loglike.temp <- samples.loglike.list
   fitted.temp <- samples.fitted.list
   Y.temp <- samples.Y.list
      for(j in 1:n.chains)
      {
      beta.temp[[j]] <- mcmc(samples.beta.list[[j]])   
      phi.temp[[j]] <- mcmc(samples.phi.list[[j]])   
      alpha.temp[[j]] <- mcmc(samples.alpha.list[[j]])   
      tau2.temp[[j]] <- mcmc(samples.tau2.list[[j]])   
      loglike.temp[[j]] <- mcmc(samples.loglike.list[[j]])   
      fitted.temp[[j]] <- mcmc(samples.fitted.list[[j]])   
      Y.temp[[j]] <- mcmc(samples.Y.list[[j]])   
      }
   beta.mcmc <- as.mcmc.list(beta.temp)
   phi.mcmc <- as.mcmc.list(phi.temp)
   alpha.mcmc <- as.mcmc.list(alpha.temp)
   tau2.mcmc <- as.mcmc.list(tau2.temp)
   fitted.mcmc <- as.mcmc.list(fitted.temp)
   Y.mcmc <- as.mcmc.list(Y.temp)
   samples <- list(beta=beta.mcmc, phi=phi.mcmc, alpha=alpha.mcmc, tau2=tau2.mcmc, fitted=fitted.mcmc, Y=Y.mcmc)

   ## Create a summary object
   n.keep <- floor((n.sample - burnin)/thin)
   summary.beta <- t(rbind(apply(samples.beta.matrix, 2, mean), apply(samples.beta.matrix, 2, quantile, c(0.025, 0.975)))) 
   summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta,p), effectiveSize(beta.mcmc), gelman.diag(beta.mcmc)$psrf[ ,2])
   rownames(summary.beta) <- colnames(X)
   colnames(summary.beta) <- c("Mean", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "PSRF (upper 95% CI)")
   summary.alpha <- t(rbind(apply(samples.alpha.matrix, 2, mean), apply(samples.alpha.matrix, 2, quantile, c(0.025, 0.975)))) 
   summary.alpha <- cbind(summary.alpha, rep(n.keep, q), rep(accept.alpha,q), effectiveSize(alpha.mcmc), gelman.diag(alpha.mcmc)$psrf[ ,2])
	    if(!is.null(names(Z)))
 	    {
 	    rownames(summary.alpha) <- names(Z)
	    }else
	    {
	    names.Z <- rep(NA,q)
		    for(j in 1:q)
		    {
		    names.Z[j] <- paste("Z[[",j, "]]", sep="")
		    }	
	    rownames(summary.alpha) <- names.Z	
	    }
   summary.hyper <- array(NA, c(1 ,7))
   summary.hyper[1, 1:3] <- c(mean(samples.tau2.matrix), quantile(samples.tau2.matrix, c(0.025, 0.975)))
   summary.hyper[1, 4:7] <- c(n.keep, accept.tau2, effectiveSize(tau2.mcmc),  gelman.diag(tau2.mcmc)$psrf[ ,2])
   summary.results <- rbind(summary.beta, summary.hyper, summary.alpha)
        if(W.binary)
        {
        alpha.min <- c(rep(NA, (p+1)), alpha.threshold)
        summary.results <- cbind(summary.results, alpha.min)    
        }else
        {}
   rownames(summary.results)[(p+1)] <- c("tau2")
   summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
   summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
        if(W.binary) summary.results[ , 8] <- round(summary.results[ , 8], 4)
   
   #### Create the posterior medians for the neighbourhood matrix W
   W.posterior <- array(NA, c(K,K))
    if(W.binary)
    {
    W.border.prob <- array(NA, c(K,K))    
    }else
    {
    W.border.prob <- NA
    }

	for(i in 1:K)
	{
		for(j in 1:K)
		{
			if(W[i,j]==1)
			{
			z.temp <- NA
				for(k in 1:q)
				{
				z.temp <- c(z.temp, Z[[k]][i,j])
				}	
			z.temp <- z.temp[-1]
			w.temp <- exp(-samples.alpha.matrix %*% z.temp)
			    if(W.binary)
			    {
			    w.posterior <- as.numeric(w.temp>=0.5)
			    W.posterior[i,j] <- ceiling(median(w.posterior))
			    W.border.prob[i,j] <- (1 - sum(w.posterior) / length(w.posterior))    
			    }else
			    {
		        W.posterior[i,j] <- median(w.temp)
			    }
			}else
			{
			}	
		}	
	}
   }



###################################
#### Compile and return the results
###################################
## Generate the dissimilarity equation
    if(q==1)
    {
    dis.eq <- rownames(summary.results)[nrow(summary.results)]    
    }else
    {
    dis.eq <- paste(rownames(summary.alpha), "+")
    len <- length(dis.eq)
    dis.eq[len] <- substr(dis.eq[len],1,nchar(dis.eq[2])-1)    
    }

    if(W.binary)
    {
    model.string <- c("Likelihood model - Poisson (log link function)", "\nRandom effects model - Binary dissimilarity CAR", "\nDissimilarity metrics - ", dis.eq, "\n")     
    }else
    {
    model.string <- c("Likelihood model - Poisson (log link function)", "\nRandom effects model - Non-binary dissimilarity CAR", "\nDissimilarity metrics - ", dis.eq, "\n")     
    }

n.total <- floor((n.sample - burnin) / thin) * n.chains
mcmc.info <- c(n.total, n.sample, burnin, thin, n.chains)
names(mcmc.info) <- c("Total samples", "n.sample", "burnin", "thin", "n.chains")
results <- list(summary.results=summary.results, samples=samples, fitted.values=fitted.values, residuals=residuals, modelfit=modelfit, accept=accept.final, localised.structure=list(W.posterior=W.posterior, W.border.prob=W.border.prob),  formula=formula, model=model.string, mcmc.info=mcmc.info, X=X)
class(results) <- "CARBayes"
    if(verbose)
    {
    b<-proc.time()
    cat("Finished in ", round(b[3]-a[3], 1), "seconds.\n")
    }else
    {}
return(results)
}

