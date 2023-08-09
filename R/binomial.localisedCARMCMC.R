binomial.localisedCARMCMC <- function(Y, failures, trials, offset, X.standardised, G, Gstar, W, K, p, burnin, n.sample, thin, MALA, n.beta.block, list.block, prior.mean.beta, prior.var.beta, prior.tau2, prior.delta, verbose, chain)
{
# Rcpp::sourceCpp("src/CARBayes.cpp")   
# source("R/common.functions.R")
# library(spdep)
# library(truncnorm)    
    
    
##########################################
#### Generate the initial parameter values
##########################################
#### Generate initial values for each chain
    if(p==0)
    {
    regression.vec <- rep(0, K)
    beta <- NA    
    }else
    {
    mod.glm <- glm(cbind(Y, failures)~X.standardised, offset=offset, family="quasibinomial")
    beta.mean <- mod.glm$coefficients[-1]
    beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))[-1]
    beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)
    regression.vec <- X.standardised %*% beta    
    }
    
theta.hat <- Y / trials
theta.hat[theta.hat==0] <- 0.01
theta.hat[theta.hat==1] <- 0.99
res.temp <- log(theta.hat / (1 - theta.hat)) - regression.vec - offset
clust <- kmeans(res.temp,G)
lambda <- clust$centers[order(clust$centers)]
Z <- rep(1, K)
    for(j in 2:G)
    {
    Z[clust$cluster==order(clust$centers)[j]] <- j    
    }
delta <- runif(1,1, min(2, prior.delta))
res.sd <- sd(res.temp, na.rm=TRUE)/5
phi <- rnorm(n=K, mean=rep(0,K), sd = res.sd)
    for(i in 1:G)
    {
    phi[which(Z==i)] <- phi[which(Z==i)] - mean(phi[which(Z==i)])
    }
tau2 <- var(phi) / 10



###################################################################
#### Compute the fitted values based on the current parameter values
####################################################################   
lp <- lambda[Z] + phi + regression.vec + offset
prob <- exp(lp)  / (1 + exp(lp))
fitted <- trials * prob   



########################################    
#### Set up the MCMC model run quantities    
#########################################
#### Matrices to store samples
n.keep <- floor((n.sample - burnin)/thin)
samples.phi <- array(NA, c(n.keep, K))
samples.tau2 <- array(NA, c(n.keep, 1))
samples.Z <- array(NA, c(n.keep, K))
samples.lambda <- array(NA, c(n.keep, G))
samples.delta <- array(NA, c(n.keep, 1))
samples.loglike <- array(NA, c(n.keep, K))
samples.fitted <- array(NA, c(n.keep, K))

#### Metropolis quantities    
    if(p>0)
    {
    samples.beta <- array(NA, c(n.keep, p))
    accept <- rep(0,8)
    proposal.sd.beta <- 0.01
    }else
    {
    accept <- rep(0,6)    
    }

proposal.sd.phi <- 0.1
proposal.sd.delta <- 0.1
proposal.sd.lambda <- 0.01
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
    ###################
    ## Sample from beta
    ###################
        if(p>0)
        {
        offset.temp <- phi + offset + lambda[Z]
            if(MALA)
            {
            temp <- binomialbetaupdateMALA(X.standardised, K, p, beta, offset.temp, Y, failures, trials, prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta, list.block)
            }else
            {
            temp <- binomialbetaupdateRW(X.standardised, K, p, beta, offset.temp, Y, failures, prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta, list.block)
            }
        beta <- temp[[1]]
        accept[7] <- accept[7] + temp[[2]]
        accept[8] <- accept[8] + n.beta.block  
        regression.vec <- X.standardised %*% beta    
        }else
        {}

     
     
     ##################
     ## Sample from phi
     ##################
     phi.offset <- regression.vec + offset  + lambda[Z]
     temp1 <- binomialcarupdateRW(Wtriplet=W.triplet, Wbegfin=W.begfin, Wtripletsum=W.triplet.sum, nsites=K, phi=phi, tau2=tau2, y=Y, failures=failures, phi_tune=proposal.sd.phi, rho=1, offset=phi.offset)
     phi <- temp1[[1]]
          for(i in 1:G)
          {
          phi[which(Z==i)] <- phi[which(Z==i)] - mean(phi[which(Z==i)])
          }
     accept[1] <- accept[1] + temp1[[2]]
     accept[2] <- accept[2] + K    

         
     
     ##################
     ## Sample from tau2
     ##################
     temp2 <- quadform(W.triplet, W.triplet.sum, n.triplet, K, phi, phi, 1)
     tau2.posterior.scale <- temp2 + prior.tau2[2] 
     tau2 <- 1 / rgamma(1, tau2.posterior.shape, scale=(1/tau2.posterior.scale))
          
    
         
    #####################
    ## Sample from lambda
    #####################
     proposal.extend <- c(-1000, lambda, 1000)
     lambda.extend <- c(-1000, lambda, 1000)
          for(i in 1:G)
          {
           proposal.extend[(i+1)] <- rtruncnorm(n=1, a=proposal.extend[i], b=proposal.extend[(i+2)], mean=lambda[i], sd=proposal.sd.lambda)    
          }
     proposal <- proposal.extend[2:(G+1)]
     lp.current <- lambda[Z] + phi + regression.vec + offset
     lp.proposal <- proposal[Z] + phi + regression.vec + offset
     prob.current <- exp(lp.current)  / (1 + exp(lp.current))
     prob.proposal <- exp(lp.proposal)  / (1 + exp(lp.proposal))
     prob1 <- sum(Y * (log(prob.proposal) - log(prob.current)) + failures * (log(1-prob.proposal) - log(1-prob.current)))          
          prob <- exp(prob1)
          if(prob > runif(1))
          {
          lambda <- proposal
          accept[3] <- accept[3] + 1  
          }else
          {}
     accept[4] <- accept[4] + 1       

         
         
     ################
     ## Sample from Z
     ################
     Z.proposal <- sample(1:G, size=K, replace=TRUE)
     prior <- delta * ((Z - Gstar)^2 - (Z.proposal-Gstar)^2)    
     lp.current <- lambda[Z] + phi + regression.vec + offset
     lp.proposal <- lambda[Z.proposal] + phi + regression.vec + offset
     prob.current <- exp(lp.current)  / (1 + exp(lp.current))
     prob.proposal <- exp(lp.proposal)  / (1 + exp(lp.proposal))
     like <- Y * (log(prob.proposal) - log(prob.current)) + failures * (log(1-prob.proposal) - log(1-prob.current))        
     prob <- exp(like + prior)   
     test <- prob> runif(K)         
     Z[test] <- Z.proposal[test]         
       
 
         
    ####################
    ## Sample from delta
    ####################
    proposal.delta <-  rtruncnorm(n=1, a=1, b=prior.delta, mean=delta, sd=proposal.sd.delta)    
    prob1 <- sum((Z-Gstar)^2) * (delta - proposal.delta)        
    prob2 <- K * log(sum(exp(-delta *(1:G - Gstar)^2))) - K * log(sum(exp(-proposal.delta *(1:G - Gstar)^2)))
    hastings <- log(dtruncnorm(x=delta, a=1, b=prior.delta, mean=proposal.delta, sd=proposal.sd.delta)) - log(dtruncnorm(x=proposal.delta, a=1, b=prior.delta, mean=delta, sd=proposal.sd.delta)) 
    prob <- exp(prob1 + prob2 + hastings)    
          if(prob > runif(1))
          {
          delta <- proposal.delta
          accept[5] <- accept[5] + 1  
          }else
          {
          }
     accept[6] <- accept[6] + 1       

         
 
    #########################
    ## Calculate the deviance
    #########################
    lp <- lambda[Z] + phi + regression.vec + offset
    prob <- exp(lp)  / (1 + exp(lp))
    fitted <- trials * prob
    loglike <- dbinom(x=Y, size=trials, prob=prob, log=TRUE)

         
         
    ###################
    ## Save the results
    ###################
        if(j > burnin & (j-burnin)%%thin==0)
        {
        ele <- (j - burnin) / thin
        samples.phi[ele, ] <- phi
        samples.lambda[ele, ] <- lambda
        samples.tau2[ele, ] <- tau2
        samples.Z[ele, ] <- Z
        samples.delta[ele, ] <- delta
        samples.loglike[ele, ] <- loglike
        samples.fitted[ele, ] <- fitted
        if(p>0) samples.beta[ele, ] <- beta    
        }else
        {
        }

         

     ########################################
     ## Self tune the acceptance probabilties
     ########################################
        if(ceiling(j/100)==floor(j/100) & j < burnin)
        {
            if(p>0)
            {
                if(p>2)
                {
                proposal.sd.beta <- common.accceptrates1(accept[7:8], proposal.sd.beta, 40, 50)
                }else
                {
                proposal.sd.beta <- common.accceptrates1(accept[7:8], proposal.sd.beta, 30, 40)    
                }
            proposal.sd.phi <- common.accceptrates1(accept[1:2], proposal.sd.phi, 40, 50)
            proposal.sd.lambda <- common.accceptrates1(accept[3:4], proposal.sd.lambda, 20, 40)
            proposal.sd.delta <- common.accceptrates2(accept[5:6], proposal.sd.delta, 40, 50, prior.delta/6)
            accept <- rep(0,8)
            }else
            {
            proposal.sd.phi <- common.accceptrates1(accept[1:2], proposal.sd.phi, 40, 50)
            proposal.sd.lambda <- common.accceptrates1(accept[3:4], proposal.sd.lambda, 20, 40)
            proposal.sd.delta <- common.accceptrates2(accept[5:6], proposal.sd.delta, 40, 50, prior.delta/6)
            accept <- rep(0,6)     
            }
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
    if(p>0)
    {
    chain.results <- list(samples.beta=samples.beta, samples.phi=samples.phi, samples.Z=samples.Z, samples.lambda=samples.lambda, samples.tau2=samples.tau2, samples.delta=samples.delta, samples.loglike=samples.loglike, samples.fitted=samples.fitted,
                    accept=accept)
    }else
    {
    chain.results <- list(samples.phi=samples.phi, samples.Z=samples.Z, samples.lambda=samples.lambda, samples.tau2=samples.tau2, samples.delta=samples.delta, samples.loglike=samples.loglike, samples.fitted=samples.fitted,
                    accept=accept)
    }

#### Return the results
return(chain.results)
}