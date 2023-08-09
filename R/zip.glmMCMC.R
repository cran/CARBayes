zip.glmMCMC <- function(Y, offset, offset.omega, X.standardised, V.standardised, K, p, q, which.miss, n.miss, which.zero, n.zero, burnin, n.sample, thin, MALA, n.beta.block, list.block, n.delta.block, list.block.delta, prior.mean.beta, prior.var.beta, prior.mean.delta, prior.var.delta, verbose, chain)
{
# Rcpp::sourceCpp("src/CARBayes.cpp")   
# source("R/common.functions.R")
##########################################
#### Generate the initial parameter values
##########################################
#### Generate initial values for each chain
mod.glm <- glm(Y[Y>0]~X.standardised[Y>0, ]-1, offset=offset[Y>0], family="quasipoisson")
beta.mean <- mod.glm$coefficients
beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)

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
fitted <- exp(as.numeric(X.standardised %*% beta) + offset)   
Y.DA <- Y

    
   
########################################    
#### Set up the MCMC model run quantities    
#########################################
#### Matrices to store samples   
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, p))
samples.delta <- array(NA, c(n.keep, q))    
samples.Z <- array(NA, c(n.keep, K))    
samples.loglike <- array(NA, c(n.keep, K))
samples.fitted <- array(NA, c(n.keep, K))
    if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))
    
    
#### Metropolis quantities
accept <- rep(0,4)
proposal.sd.beta <- 0.01
proposal.sd.delta <- 0.01


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
        Y.DA[which.miss==0] <- rpois(n=n.miss, lambda=fitted[which.miss==0]) * (1-Z[which.miss==0]) 
        }else
        {}
    which.zero <- which(Y.DA==0)
    n.zero <- length(which.zero)    
        
    
    
    ###################################
    #### Update Z via data augmentation
    ###################################
    prob.pointmass <- omega[which.zero] / (omega[which.zero] + (1 - omega[which.zero]) * exp(-exp(as.matrix(X.standardised[which.zero, ]) %*% beta + offset[which.zero])))
    Z <-  rep(0, K)
    Z[which.zero] <- rbinom(n=n.zero, size=1, prob=prob.pointmass)    
    
    
    
    ####################
    ## Sample from beta
    ####################
    Z.zero <- which(Z==0)
    offset.temp <- offset[Z.zero]
    
        if(MALA)
        {
        temp <- poissonbetaupdateMALA(as.matrix(X.standardised[Z.zero, ]), length(Z.zero), p, beta, offset.temp, Y.DA[Z.zero], prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta, list.block)
        }else
        {
        temp <- poissonbetaupdateRW(as.matrix(X.standardised[Z.zero, ]), length(Z.zero), p, beta, offset.temp, Y.DA[Z.zero], prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta, list.block)
        }
    beta <- temp[[1]]
    accept[1] <- accept[1] + temp[[2]]
    accept[2] <- accept[2] + n.beta.block  

    
        
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
    accept[3] <- accept[3] + temp[[2]]
    accept[4] <- accept[4] + n.delta.block  
    omega <- exp(V.standardised %*% delta+offset.omega) / (1+exp(V.standardised %*% delta+offset.omega))
        
        
        
    #########################
    ## Calculate the deviance
    #########################
    lp <- as.numeric(X.standardised %*% beta) + offset
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
            proposal.sd.delta <- common.accceptrates1(accept[3:4], proposal.sd.delta, 40, 50)
            }else
            {
            proposal.sd.delta <- common.accceptrates1(accept[3:4], proposal.sd.delta, 30, 40)    
            }
            
        accept <- rep(0,4)
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
    close(progressBar)
    }else
    {}



############################################
#### Return the results to the main function
############################################
#### Compile the results
    if(n.miss==0) samples.Y = NA
chain.results <- list(samples.beta=samples.beta, samples.delta=samples.delta, samples.Z=samples.Z, samples.loglike=samples.loglike, samples.fitted=samples.fitted,
                    samples.Y=samples.Y, accept=accept)

#### Return the results
return(chain.results)
}