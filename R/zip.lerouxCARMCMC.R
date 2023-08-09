zip.lerouxCARMCMC <- function(Y, offset, offset.omega, X.standardised, V.standardised, W, rho, fix.rho, K, p, q, which.miss, n.miss, which.zero, n.zero, burnin, n.sample, thin, MALA, n.beta.block, list.block, n.delta.block, list.block.delta, prior.mean.beta, prior.var.beta, prior.mean.delta, prior.var.delta, prior.tau2, verbose, chain)
{
    # Rcpp::sourceCpp("src/CARBayes.cpp")   
    # source("R/common.functions.R")
    # library(spdep)
    # library(truncnorm)    
    # 
    # 
    ##########################################
    #### Generate the initial parameter values
    ##########################################
    #### Initial parameter values
    mod.glm <- glm(Y[Y>0]~X.standardised[Y>0, ]-1, offset=offset[Y>0], family="quasipoisson")
    beta.mean <- mod.glm$coefficients
    beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
    beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)
    
    log.Y <- log(Y)
    log.Y[Y==0] <- -0.1  
    res.temp <- log.Y - X.standardised %*% beta.mean - offset
    res.sd <- sd(res.temp, na.rm=TRUE)/5
    phi <- rnorm(n=K, mean=rep(0,K), sd=res.sd)
    tau2 <- var(phi) / 10

    Y.zero <- rep(0,K)
    Y.zero[which.zero] <- 1
    mod.glm2 <- glm(Y.zero~V.standardised-1, offset=offset.omega, family="binomial")
    delta.mean <- mod.glm2$coefficients
    delta.sd <- sqrt(diag(summary(mod.glm2)$cov.scaled))
    delta <- rnorm(n=length(delta.mean), mean=delta.mean, sd=delta.sd)
    
    omega <- exp(V.standardised %*% delta+offset.omega) / (1+exp(V.standardised %*% delta+offset.omega))
    prob.pointmass <- omega[which.zero] / (omega[which.zero]+(1-omega[which.zero])*exp(-exp(as.matrix(X.standardised[which.zero, ]) %*% beta + offset[which.zero])))
    Z <-  rep(0, K)
    Z[which.zero] <- rbinom(n=n.zero, size=1, prob=prob.pointmass)    

    
    
    ###################################################################
    #### Compute the fitted values based on the current parameter values
    ####################################################################   
    fitted <- exp(as.numeric(X.standardised %*% beta) + offset + phi)     
    Y.DA <- Y
    
    
    
    ########################################    
    #### Set up the MCMC model run quantities    
    #########################################
    #### Matrices to store samples   
    n.keep <- floor((n.sample - burnin)/thin)
    samples.beta <- array(NA, c(n.keep, p))
    samples.phi <- array(NA, c(n.keep, K))
    samples.tau2 <- array(NA, c(n.keep, 1))
    if(!fix.rho) samples.rho <- array(NA, c(n.keep, 1))
    samples.delta <- array(NA, c(n.keep, q))    
    samples.Z <- array(NA, c(n.keep, K))    
    samples.loglike <- array(NA, c(n.keep, K))
    samples.fitted <- array(NA, c(n.keep, K))
    if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))
    
    
    #### Metropolis quantities
    accept <- rep(0,8)
    proposal.sd.beta <- 0.01
    proposal.sd.delta <- 0.01
    proposal.sd.phi <- 0.1
    proposal.sd.rho <- 0.02
    tau2.posterior.shape <- prior.tau2[1] + 0.5 * K
    
    
    ##################################
    #### Set up the spatial quantities
    ##################################
    #### CAR quantities
    W.quants <- common.Wcheckformat(W)
    W <- W.quants$W
    W.triplet <- W.quants$W.triplet
    n.triplet <- W.quants$n.triplet
    W.triplet.sum <- W.quants$W.triplet.sum
    n.neighbours <- W.quants$n.neighbours 
    W.begfin <- W.quants$W.begfin
    
    
    #### Create the determinant     
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
    #### Create the MCMC samples      
    for(j in 1:n.sample)
    {
        ####################################
        ## Sample from Y - data augmentation
        ####################################
        if(n.miss>0)
        {
            Y.DA[which.miss==0] <- rpois(n=n.miss, lambda=fitted[which.miss==0]) * (1-Z[which.miss==0]) 
        }else
        {}
        which.zero <- which(Y.DA==0)
        n.zero <- length(which.zero)   
        
        
        
        ###################################
        #### Update Z via data augmentation
        ###################################
        prob.pointmass <- omega[which.zero] / (omega[which.zero] + (1 - omega[which.zero]) * exp(-exp(as.matrix(X.standardised[which.zero, ]) %*% beta + offset[which.zero] + phi[which.zero])))
        Z <-  rep(0, K)
        Z[which.zero] <- rbinom(n=n.zero, size=1, prob=prob.pointmass)    
        
        
        
        ####################
        ## Sample from beta
        ####################
        Z.zero <- which(Z==0)
        offset.temp <- offset[Z.zero] + phi[Z.zero]
        
        if(MALA)
        {
            temp <- poissonbetaupdateMALA(as.matrix(X.standardised[Z.zero, ]), length(Z.zero), p, beta, offset.temp, Y.DA[Z.zero], prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta, list.block)
        }else
        {
            temp <- poissonbetaupdateRW(as.matrix(X.standardised[Z.zero, ]), length(Z.zero), p, beta, offset.temp, Y.DA[Z.zero], prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta, list.block)
        }
        beta <- temp[[1]]
        regression <- X.standardised %*% beta
        accept[1] <- accept[1] + temp[[2]]
        accept[2] <- accept[2] + n.beta.block  
        
        
        
        ####################
        ## Sample from phi
        ####################
        beta.offset <- regression + offset
        temp1 <- zipcarupdateRW(Wtriplet=W.triplet, Wbegfin=W.begfin, W.triplet.sum, nsites=K, phi=phi, tau2=tau2, y=Y.DA, phi_tune=proposal.sd.phi, rho=rho, offset=beta.offset, 1-Z)
        phi <- temp1[[1]]
        if(rho<1)
        {
            phi <- phi - mean(phi)
        }else
        {
            phi[which(islands==1)] <- phi[which(islands==1)] - mean(phi[which(islands==1)])   
        }
        accept[3] <- accept[3] + temp1[[2]]
        accept[4] <- accept[4] + sum(Z==0)    
        
        
        
        
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
                accept[5] <- accept[5] + 1           
            }else
            {
            }              
            accept[6] <- accept[6] + 1           
        }else
        {}
        
        
        
        ######################    
        #### Sample from delta
        ######################
        offset.temp <- offset.omega
        if(MALA)
        {
            temp <- binomialbetaupdateMALA(V.standardised, K, q, delta, offset.temp, Z, 1-Z, rep(1,K), prior.mean.delta, prior.var.delta, n.delta.block, proposal.sd.delta, list.block.delta)
        }else
        {
            temp <- binomialbetaupdateRW(V.standardised, K, q, delta, offset.temp, Z, 1-Z, prior.mean.delta, prior.var.delta, n.delta.block, proposal.sd.delta, list.block.delta)
        }
        delta <- temp[[1]]
        accept[7] <- accept[7] + temp[[2]]
        accept[8] <- accept[8] + n.delta.block  
        omega <- exp(V.standardised %*% delta+offset.omega) / (1+exp(V.standardised %*% delta+offset.omega))
        
        
        
        #########################
        ## Calculate the deviance
        #########################
        lp <- as.numeric(regression) + phi + offset
        fitted <- exp(lp)
        fitted.zip <- fitted * (1-omega)    
        temp <- rep(0,K)
        temp[Z==1] <- log(omega[Z==1])
        loglike <- temp + (1-Z) * (log(1-omega) + dpois(x=as.numeric(Y), lambda=fitted, log=T))
        
        
        
        ###################
        ## Save the results
        ###################
        if(j > burnin & (j-burnin)%%thin==0)
        {
            ele <- (j - burnin) / thin
            samples.beta[ele, ] <- beta
            samples.phi[ele, ] <- phi
            samples.tau2[ele, ] <- tau2
            if(!fix.rho) samples.rho[ele, ] <- rho
            samples.delta[ele, ] <- delta
            samples.Z[ele, ] <- Z
            samples.loglike[ele, ] <- loglike
            samples.fitted[ele, ] <- fitted.zip
            if(n.miss>0) samples.Y[ele, ] <- Y.DA[which.miss==0]
        }else
        {}
        
        
        
        ########################################
        ## Self tune the acceptance probabilties
        ########################################
        if(ceiling(j/100)==floor(j/100) & j < burnin)
        {
            #### Update the proposal sds
            ## beta
            if(p>2)
            {
                proposal.sd.beta <- common.accceptrates1(accept[1:2], proposal.sd.beta, 40, 50)
            }else
            {
                proposal.sd.beta <- common.accceptrates1(accept[1:2], proposal.sd.beta, 30, 40)    
            }
            
            ## delta
            if(q>2)
            {
                proposal.sd.delta <- common.accceptrates1(accept[7:8], proposal.sd.delta, 40, 50)
            }else
            {
                proposal.sd.delta <- common.accceptrates1(accept[7:8], proposal.sd.delta, 30, 40)    
            }
            
            proposal.sd.phi <- common.accceptrates1(accept[3:4], proposal.sd.phi, 40, 50)
            if(!fix.rho) proposal.sd.rho <- common.accceptrates2(accept[5:6], proposal.sd.rho, 40, 50, 0.5)
            accept <- rep(0,8)
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
    chain.results <- list(samples.beta=samples.beta, samples.phi=samples.phi, samples.tau2=samples.tau2, samples.rho=samples.rho, samples.delta=samples.delta, samples.Z=samples.Z, samples.loglike=samples.loglike, samples.fitted=samples.fitted,
                          samples.Y=samples.Y, accept=accept)
    
    #### Return the results
    return(chain.results)
}