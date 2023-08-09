gaussian.multilevelCARMCMC <- function(Y, offset, X.standardised, W, ind.area, rho, fix.rho, n, K, p, which.miss, n.miss, burnin, n.sample, thin, prior.mean.beta, prior.var.beta, prior.nu2, prior.tau2, verbose, chain)
{
    # Rcpp::sourceCpp("src/CARBayes.cpp")   
    # source("R/common.functions.R")
    # library(spdep)
    # library(truncnorm)    
    
    
    ##########################################
    #### Generate the initial parameter values
    ##########################################
    #### Generate initial values for each chain
    mod.glm <- lm(Y~X.standardised-1, offset=offset)
    beta.mean <- mod.glm$coefficients
    beta.sd <- sqrt(diag(summary(mod.glm)$cov.unscaled)) * summary(mod.glm)$sigma
    beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)
    res.temp <- Y - X.standardised %*% beta.mean - offset
    res.sd <- sd(res.temp, na.rm=TRUE)/5
    phi <- rnorm(n=K, mean=rep(0,K), sd=res.sd)
    phi.extend <- phi[ind.area]
    tau2 <- var(phi) / 10
    nu2 <- tau2

    
    
    ###################################################################
    #### Compute the fitted values based on the current parameter values
    ####################################################################   
    fitted <- as.numeric(X.standardised %*% beta) + phi.extend + offset
    Y.DA <- Y
    
    
    
    ########################################    
    #### Set up the MCMC model run quantities    
    #########################################
    #### Ind.area parts
    ind.area.list <- as.list(rep(0,K))
    n.individual <- rep(0,K)
    n.individual.miss <- rep(0,K)
    for(r in 1:K)
    {
        ind.area.list[[r]] <- which(ind.area==r)
        n.individual[r] <- length(ind.area.list[[r]])
        n.individual.miss[r] <- sum(which.miss[ind.area.list[[r]]])
    }
    
    
    #### Matrices to store samples
    n.keep <- floor((n.sample - burnin)/thin)
    samples.beta <- array(NA, c(n.keep, p))
    samples.phi <- array(NA, c(n.keep, K))
    samples.nu2 <- array(NA, c(n.keep, 1))
    samples.tau2 <- array(NA, c(n.keep, 1))
    if(!fix.rho) samples.rho <- array(NA, c(n.keep, 1))
    samples.loglike <- array(NA, c(n.keep, n))
    samples.fitted <- array(NA, c(n.keep, n))
    if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))
    
    
    #### Metropolis quantities
    accept <- rep(0,2)
    proposal.sd.rho <- 0.02
    tau2.posterior.shape <- prior.tau2[1] + 0.5*K
    nu2.posterior.shape <- prior.nu2[1] + 0.5*n

    
    #### CAR quantities
    W.quants <- common.Wcheckformat(W)
    W <- W.quants$W
    W.triplet <- W.quants$W.triplet
    n.triplet <- W.quants$n.triplet
    W.triplet.sum <- W.quants$W.triplet.sum
    n.neighbours <- W.quants$n.neighbours 
    W.begfin <- W.quants$W.begfin
    
    
    #### Create the spatial determinant     
    if(!fix.rho)
    {
        Wstar <- diag(apply(W,1,sum)) - W
        Wstar.eigen <- eigen(Wstar)
        Wstar.val <- Wstar.eigen$values
        det.Q <- 0.5 * sum(log((rho * Wstar.val + (1-rho))))    
    }else
    {}
    
    
    #### Check for islands
    W.list<- mat2listw(W, style = "B")
    W.nb <- W.list$neighbours
    W.islands <- n.comp.nb(W.nb)
    islands <- W.islands$comp.id
    n.islands <- max(W.islands$nc)
    if(rho==1) tau2.posterior.shape <- prior.tau2[1] + 0.5 * (K-n.islands)   
    
    
    #### Beta update quantities
    data.precision.beta <- t(X.standardised) %*% X.standardised
    if(length(prior.var.beta)==1)
    {
        prior.precision.beta <- 1 / prior.var.beta
    }else
    {
        prior.precision.beta <- solve(diag(prior.var.beta))
    }
    
    
    #### Start timer
    if(verbose)
    {
        cat("\nMarkov chain", chain,  "- generating", n.keep, "post burnin and thinned samples.\n", sep = " ")
        progressBar <- txtProgressBar(style = 3)
        percentage.points<-round((1:100/100)*n.sample)
    }else
    {
        percentage.points<-round((1:100/100)*n.sample)     
    }
    
    
    
    ######################
    #### Run an MCMC chain
    ######################
    for(j in 1:n.sample)
    {
        ####################################
        ## Sample from Y - data augmentation
        ####################################
        if(n.miss>0)
        {
            Y.DA[which.miss==0] <- rnorm(n=n.miss, mean=fitted[which.miss==0], sd=sqrt(nu2))    
        }else
        {}
        
        
        
        ####################
        ## Sample from beta
        ####################
        fc.precision <- prior.precision.beta + data.precision.beta / nu2
        fc.var <- solve(fc.precision)
        beta.offset <- as.numeric(Y.DA - offset - phi.extend)
        beta.offset2 <- t(X.standardised) %*% beta.offset / nu2 + prior.precision.beta %*% prior.mean.beta
        fc.mean <- fc.var %*% beta.offset2
        chol.var <- t(chol(fc.var))
        beta <- fc.mean + chol.var %*% rnorm(p)        
        
        
        
        ##################
        ## Sample from nu2
        ##################
        fitted.current <-  as.numeric(X.standardised %*% beta) + phi.extend + offset
        nu2.posterior.scale <- prior.nu2[2] + 0.5 * sum((Y.DA - fitted.current)^2)
        nu2 <- 1 / rgamma(1, nu2.posterior.shape, scale=(1/nu2.posterior.scale))    
        
        
        
        ####################
        ## Sample from phi
        ####################
        offset.phi <- (Y.DA - as.numeric(X.standardised %*% beta) - offset) / nu2
        offset.phi2 <- tapply(offset.phi, ind.area, sum, na.rm=T)
        phi <- gaussiancarmultilevelupdate(Wtriplet=W.triplet, Wbegfin=W.begfin, W.triplet.sum, n_individual=n.individual, nsites=K, phi=phi, tau2=tau2, rho=rho, nu2=nu2, offset=offset.phi2)
        if(rho<1)
        {
            phi <- phi - mean(phi)
        }else
        {
            phi[which(islands==1)] <- phi[which(islands==1)] - mean(phi[which(islands==1)])   
        }
        phi.extend <- phi[ind.area]
        
        
        
        ##################
        ## Sample from tau2
        ##################
        temp2 <- quadform(W.triplet, W.triplet.sum, n.triplet, K, phi, phi, rho)
        tau2.posterior.scale <- temp2 + prior.tau2[2] 
        tau2 <- 1 / rgamma(1, tau2.posterior.shape, scale=(1/tau2.posterior.scale))
        
        
        
        ##################
        ## Sample from rho
        ##################
        if(!fix.rho)
        {
            proposal.rho <- rtruncnorm(n=1, a=0, b=1, mean=rho, sd=proposal.sd.rho)  
            temp3 <- quadform(W.triplet, W.triplet.sum, n.triplet, K, phi, phi, proposal.rho)
            det.Q.proposal <- 0.5 * sum(log((proposal.rho * Wstar.val + (1-proposal.rho))))              
            logprob.current <- det.Q - temp2 / tau2
            logprob.proposal <- det.Q.proposal - temp3 / tau2
            hastings <- log(dtruncnorm(x=rho, a=0, b=1, mean=proposal.rho, sd=proposal.sd.rho)) - log(dtruncnorm(x=proposal.rho, a=0, b=1, mean=rho, sd=proposal.sd.rho)) 
            prob <- exp(logprob.proposal - logprob.current + hastings)
            
            #### Accept or reject the proposal
            if(prob > runif(1))
            {
                rho <- proposal.rho
                det.Q <- det.Q.proposal
                accept[1] <- accept[1] + 1           
            }else
            {}              
            accept[2] <- accept[2] + 1           
        }else
        {}
        
        
        
        #########################
        ## Calculate the deviance
        #########################
        fitted <- as.numeric(X.standardised %*% beta) + phi.extend + offset 
        loglike <- dnorm(Y, mean = fitted, sd = rep(sqrt(nu2),n), log=TRUE)
        
        
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
            if(!fix.rho) samples.rho[ele, ] <- rho
            samples.loglike[ele, ] <- loglike
            samples.fitted[ele, ] <- fitted
            if(n.miss>0) samples.Y[ele, ] <- Y.DA[which.miss==0]
        }else
        {}
        
        
        
        #######################################
        #### Update the acceptance rate for rho
        #######################################
        if(ceiling(j/100)==floor(j/100) & j < burnin)
        {
            #### Update the proposal sds
            if(!fix.rho)
            {
                proposal.sd.rho <- common.accceptrates2(accept[1:2], proposal.sd.rho, 40, 50, 0.5)
            }
            accept <- c(0,0)
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
        close(progressBar)
    }else
    {}
    
    
    
    ############################################
    #### Return the results to the main function
    ############################################
    #### Compile the results
    if(n.miss==0) samples.Y = NA
    if(fix.rho) samples.rho=NA
    chain.results <- list(samples.beta=samples.beta, samples.phi=samples.phi, samples.tau2=samples.tau2, samples.nu2=samples.nu2, samples.rho=samples.rho, samples.loglike=samples.loglike, samples.fitted=samples.fitted,
                          samples.Y=samples.Y, accept=accept)
    
    #### Return the results
    return(chain.results)
}