binomial.dissimilarityCARMCMC <- function(Y, failures, trials, offset, X.standardised, Z, W.binary, W, K, p, q, which.miss, n.miss, burnin, n.sample, thin, MALA, n.beta.block, list.block, prior.mean.beta, prior.var.beta, prior.tau2, alpha.max, verbose, chain)
{
# Rcpp::sourceCpp("src/CARBayes.cpp")   
# source("R/common.functions.R")
# library(spdep)
# library(truncnorm)    
# library(spam)
#     
#     
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
alpha <- runif(n=q, min=rep(0,q), max=rep(alpha.max/(2+q)))
lp <- as.numeric(X.standardised %*% beta) + phi + offset
prob <- exp(lp)  / (1 + exp(lp))



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
samples.alpha <- array(NA, c(n.keep, q))
samples.loglike <- array(NA, c(n.keep, K))
samples.fitted <- array(NA, c(n.keep, K))
    if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))


## Metropolis quantities
accept <- rep(0,6)
proposal.sd.alpha <- 0.02 * alpha.max
proposal.sd.beta <- 0.01
proposal.sd.phi <- 0.1
tau2.posterior.shape <- prior.tau2[1] + 0.5 * K



##################################
#### Set up the spatial quantities
##################################
#### CAR quantities
W.quants <- common.Wcheckformat.disimilarity(W)
W <- W.quants$W
W.triplet <- W.quants$W.triplet
n.triplet <- W.quants$n.triplet
W.triplet.sum <- W.quants$W.triplet.sum
n.neighbours <- W.quants$n.neighbours 
W.begfin <- W.quants$W.begfin
spam.W <- W.quants$spam.W   


#### Create the Z triplet form
Z.triplet <- array(NA, c(n.triplet, q))
     for(i in 1:n.triplet)
     {
     row <- W.triplet[i,1]
     col <- W.triplet[i,2]
          for(j in 1:q)
          {
          Z.triplet[i,j] <- Z[[j]][row, col]     
          }     
     }

    if(W.binary)
    {
    W.triplet[ ,3] <- as.numeric(exp(-Z.triplet %*% alpha)>=0.5)        
    }else
    {
    W.triplet[ ,3] <- as.numeric(exp(-Z.triplet %*% alpha))    
    }
W.triplet.sum <- tapply(W.triplet[ ,3], W.triplet[ ,1], sum)
spam.W@entries <- W.triplet[ ,3]      
spam.Wprop <- spam.W     
W.tripletprop <- W.triplet

 
#### Create the matrix form of Q
rho <- 0.99
Q <- -rho * spam.W     
diag(Q) <- rho * rowSums(spam.W) + 1-rho
det.Q <- sum(log(diag(chol.spam(Q))))     


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
    phi <- phi - mean(phi)
    accept[3] <- accept[3] + temp1[[2]]
    accept[4] <- accept[4] + K

 		

    ##################
	## Sample from tau2
	##################
    temp2 <- quadform(W.triplet, W.triplet.sum, n.triplet, K, phi, phi, rho)
    tau2.posterior.scale <- temp2 + prior.tau2[2] 
    tau2 <- 1 / rgamma(1, tau2.posterior.shape, scale=(1/tau2.posterior.scale))		
    
 
               
    ######################
	#### Sample from alpha
	######################
    ## Propose a value
    proposal.alpha <- alpha
        for(r in 1:q)
    	{
    	proposal.alpha[r] <- rtruncnorm(n=1, a=0, b=alpha.max[r],  mean=alpha[r], sd=proposal.sd.alpha[r])
    	}
               
    ## Create the proposal values for W and Q
        if(W.binary)
        {
        W.tripletprop[ ,3] <- as.numeric(exp(-Z.triplet %*% proposal.alpha)>=0.5)        
        }else
        {
        W.tripletprop[ ,3] <- as.numeric(exp(-Z.triplet %*% proposal.alpha))    
        }
    W.triplet.sum.prop <- tapply(W.tripletprop[ ,3], W.tripletprop[ ,1], sum)
    spam.Wprop@entries <- W.tripletprop[ ,3]     
    Qprop <- -rho * spam.Wprop 
    diag(Qprop) <- rho * rowSums(spam.Wprop) + 1-rho
    det.Qprop <- sum(log(diag(chol.spam(Qprop))))     
    temp3 <- quadform(W.tripletprop, W.triplet.sum.prop, n.triplet, K, phi, phi, rho)              
               
    #### Calculate the acceptance probability
    logprob.current <- det.Q - temp2 / tau2
    logprob.proposal <- det.Qprop - temp3 / tau2
    hastings <- sum(log(dtruncnorm(x=alpha, a=rep(0,q), b=alpha.max, mean=proposal.alpha, sd=proposal.sd.alpha)) - log(dtruncnorm(x=proposal.alpha, a=rep(0,q), b=alpha.max, mean=alpha, sd=proposal.sd.alpha)))
    prob <- exp(logprob.proposal - logprob.current + hastings)

    
    #### Accept or reject the proposed value
        if(prob > runif(1))
        {
    	alpha <- proposal.alpha
    	det.Q <- det.Qprop 
        W.triplet[ ,3] <- W.tripletprop[ ,3]
        W.triplet.sum <- W.triplet.sum.prop
        accept[5] <- accept[5] + 1
    	}else
    	{}
    accept[6] <- accept[6] + 1    
               
               
               
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
        samples.alpha[ele, ] <- alpha
        samples.loglike[ele, ] <- loglike
        samples.fitted[ele, ] <- fitted
            if(n.miss>0) samples.Y[ele, ] <- Y.DA[which.miss==0]               
        }else
        {}


    
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
    	proposal.sd.alpha <- common.accceptrates2(accept[5:6], proposal.sd.alpha, 40, 50, alpha.max/4)
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
chain.results <- list(samples.beta=samples.beta, samples.phi=samples.phi, samples.tau2=samples.tau2, samples.alpha=samples.alpha, samples.loglike=samples.loglike, samples.fitted=samples.fitted,
                    samples.Y=samples.Y, accept=accept)

#### Return the results
return(chain.results)
}