binomial.lerouxCARMCMC <- function(Y, failures, trials, offset, X.standardised, W, rho, fix.rho, K, p, which.miss, n.miss, burnin, n.sample, thin, MALA, n.beta.block, list.block, prior.mean.beta, prior.var.beta, prior.tau2, verbose, chain)
{
# Rcpp::sourceCpp("src/CARBayes.cpp")   
# source("R/common.functions.R")
# library(spdep)
# library(truncnorm)    
    
    
##########################################
#### Generate the initial parameter values
##########################################
#### Generate initial values for each chain
dat <- cbind(Y, failures)
mod.glm <- glm(dat~X.standardised-1, offset=offset, family="quasibinomial")
beta.mean <- mod.glm$coefficients
beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)
    
theta.hat <- Y / trials
theta.hat[theta.hat==0] <- 0.01
theta.hat[theta.hat==1] <- 0.99
res.temp <- log(theta.hat / (1 - theta.hat)) - X.standardised %*% beta.mean - offset
res.sd <- sd(res.temp, na.rm=TRUE)/5
phi <- rnorm(n=K, mean=rep(0,K), sd=res.sd)
tau2 <- var(phi) / 10


   
###################################################################
#### Compute the fitted values based on the current parameter values
####################################################################   
lp <- as.numeric(X.standardised %*% beta) + phi + offset
prob <- exp(lp)  / (1 + exp(lp))
fitted <- trials * prob
Y.DA <- Y
failures.DA <- trials - Y.DA

    
   
########################################    
#### Set up the MCMC model run quantities    
#########################################
#### Matrices to store samples
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, p))
samples.phi <- array(NA, c(n.keep, K))
samples.tau2 <- array(NA, c(n.keep, 1))
    if(!fix.rho) samples.rho <- array(NA, c(n.keep, 1))
samples.loglike <- array(NA, c(n.keep, K))
samples.fitted <- array(NA, c(n.keep, K))
    if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))
    

#### Metropolis quantities
accept <- rep(0,6)
proposal.sd.beta <- 0.01
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
        Y.DA[which.miss==0] <- rbinom(n=n.miss, size=trials[which.miss==0], prob=prob[which.miss==0])
        failures.DA <- trials - Y.DA
        }else
        {}
    
            
        
    ####################
    ## Sample from beta
    ####################
    offset.temp <- phi + offset
        if(MALA)
        {
        temp <- binomialbetaupdateMALA(X.standardised, K, p, beta, offset.temp, Y.DA, failures.DA, trials, prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta, list.block)
        }else
        {
        temp <- binomialbetaupdateRW(X.standardised, K, p, beta, offset.temp, Y.DA, failures.DA, prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta, list.block)
        }
    beta <- temp[[1]]
    accept[1] <- accept[1] + temp[[2]]
    accept[2] <- accept[2] + n.beta.block  
        
    
        
    ####################
    ## Sample from phi
    ####################
    beta.offset <- X.standardised %*% beta + offset
    temp1 <- binomialcarupdateRW(Wtriplet=W.triplet, Wbegfin=W.begfin, Wtripletsum=W.triplet.sum, nsites=K, phi=phi, tau2=tau2, y=Y.DA, failures=failures.DA, phi_tune=proposal.sd.phi, rho=rho, offset=beta.offset)
    phi <- temp1[[1]]
        if(rho<1)
        {
        phi <- phi - mean(phi)
        }else
        {
        phi[which(islands==1)] <- phi[which(islands==1)] - mean(phi[which(islands==1)])   
        }
    accept[3] <- accept[3] + temp1[[2]]
    accept[4] <- accept[4] + K
        
        
    
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
        
        
    
        #########################
        ## Calculate the deviance
        #########################
        lp <- as.numeric(X.standardised %*% beta) + phi + offset
        prob <- exp(lp)  / (1 + exp(lp))
        fitted <- trials * prob
        loglike <- dbinom(x=Y, size=trials, prob=prob, log=TRUE)
        
        
        
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
            samples.loglike[ele, ] <- loglike
            samples.fitted[ele, ] <- fitted
                if(n.miss>0) samples.Y[ele, ] <- Y.DA[which.miss==0]
            }else
            {
            }
        
        
        
        ########################################
        ## Self tune the acceptance probabilties
        ########################################
            if(ceiling(j/100)==floor(j/100) & j < burnin)
            {
            #### Update the proposal sds
                if(p>2)
                {
                proposal.sd.beta <- common.accceptrates1(accept[1:2], proposal.sd.beta, 40, 50)
                }else
                {
                proposal.sd.beta <- common.accceptrates1(accept[1:2], proposal.sd.beta, 30, 40)    
                }
            proposal.sd.phi <- common.accceptrates1(accept[3:4], proposal.sd.phi, 40, 50)
                if(!fix.rho)
                {
                proposal.sd.rho <- common.accceptrates2(accept[5:6], proposal.sd.rho, 40, 50, 0.5)
                }
            accept <- c(0,0,0,0,0,0)
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
chain.results <- list(samples.beta=samples.beta, samples.phi=samples.phi, samples.tau2=samples.tau2, samples.rho=samples.rho, samples.loglike=samples.loglike, samples.fitted=samples.fitted,
                    samples.Y=samples.Y, accept=accept)

#### Return the results
return(chain.results)
}