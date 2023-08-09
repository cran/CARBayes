gaussian.MVlerouxCARMCMC <- function(Y, offset, X.standardised, W, rho, fix.rho, K, p, J, N.all, which.miss, n.miss, miss.locator, burnin, n.sample, thin, prior.mean.beta, prior.var.beta, prior.nu2, prior.Sigma.df, prior.Sigma.scale, verbose, chain)
{
    # Rcpp::sourceCpp("src/CARBayes.cpp")   
    # source("R/common.functions.R")
    # library(spdep)
    # library(truncnorm)  
    # library(MCMCpack)
##########################################
#### Generate the initial parameter values
##########################################
beta <- array(NA, c(p, J))
nu2 <- rep(NA, J)
    for(i in 1:J)
    {
    mod.glm <- lm(Y[ ,i]~X.standardised-1, offset=offset[ ,i])
    beta.mean <- mod.glm$coefficients
    beta.sd <- sqrt(diag(summary(mod.glm)$cov.unscaled)) * summary(mod.glm)$sigma
    beta[ ,i] <- rnorm(n=p, mean=beta.mean, sd=beta.sd)
    nu2[i] <- runif(1, var(mod.glm$residuals)*0.5, var(mod.glm$residuals)*2)
    }
    
res.temp <- Y - X.standardised %*% beta - offset
res.sd <- sd(res.temp, na.rm=TRUE)/5
phi.vec <- rnorm(n=N.all, mean=0, sd=res.sd)
phi <- matrix(phi.vec, nrow=K, byrow=TRUE)
Sigma <- cov(phi)
Sigma.inv <- solve(Sigma)
Sigma.a <- rep(1, J)



    
####################################################################
#### Compute the fitted values based on the current parameter values
####################################################################   
regression <- X.standardised %*% beta
fitted <- regression + phi + offset
Y.DA <- Y
    
    
    
###############################    
#### Set up the MCMC quantities    
###############################
#### Matrices to store samples   
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, J*p))
samples.nu2 <- array(NA, c(n.keep, J))
samples.phi <- array(NA, c(n.keep, N.all))
samples.Sigma <- array(NA, c(n.keep, J, J))
samples.Sigma.a <- array(NA, c(n.keep, J))
    if(!fix.rho) samples.rho <- array(NA, c(n.keep, 1))
samples.loglike <- array(NA, c(n.keep, N.all))
samples.fitted <- array(NA, c(n.keep, N.all))
    if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))
    
    
#### Metropolis quantities
accept <- rep(0,4)
proposal.sd.phi <- 0.1
proposal.sd.rho <- 0.02
nu2.posterior.shape <- prior.nu2[1] + 0.5 * K
Sigma.post.df <- prior.Sigma.df + K  + J - 1
Sigma.a.post.shape <- (prior.Sigma.df + J) / 2



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
Wstar <- diag(apply(W,1,sum)) - W
Q <- rho * Wstar + diag(rep(1-rho,K))
    
    
#### Create the determinant     
    if(!fix.rho)
    {
    Wstar.eigen <- eigen(Wstar)
    Wstar.val <- Wstar.eigen$values
    det.Q <- sum(log((rho * Wstar.val + (1-rho))))    
    }else
    {} 
    
    
#### Check for islands
W.list<- mat2listw(W, style = "B")
W.nb <- W.list$neighbours
W.islands <- n.comp.nb(W.nb)
islands <- W.islands$comp.id
islands.all <- rep(islands,J)
n.islands <- max(W.islands$nc)
    if(rho==1) Sigma.post.df <- prior.Sigma.df + K  + J - 1 - n.islands       
    

#### Specify vector variants
Y.vec <- as.numeric(t(Y))


#### Beta update quantities
data.precision <- t(X.standardised) %*% X.standardised 
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
#### Create the MCMC samples      
   for(j in 1:n.sample)
    {
    ####################################
    ## Sample from Y - data augmentation
    ####################################
        if(n.miss>0)
        {
        nu.mat <- matrix(rep(sqrt(nu2), K), nrow=K, byrow=T)
        Y.DA[miss.locator] <- rnorm(n=n.miss, mean=fitted[miss.locator], sd=nu.mat[miss.locator])    
        }else
        {}
        
        
        
    ###################
    ## Sample from beta
    ###################
        for(r in 1:J)
        {
        fc.precision <- prior.precision.beta + data.precision / nu2[r]
        fc.var <- solve(fc.precision)    
        fc.temp1 <- t(((Y.DA[, r] - phi[ , r] - offset[ , r]) %*%  X.standardised) / nu2[r]) + prior.precision.beta %*% prior.mean.beta  
        fc.mean <- fc.var %*% fc.temp1
        chol.var <- t(chol(fc.var))
        beta[ ,r] <- fc.mean + chol.var %*% rnorm(p)   
        }
    regression <- X.standardised %*% beta    

        
    
        ##################
        ## Sample from nu2
        ##################
        fitted.current <- regression + phi + offset
        nu2.posterior.scale <- prior.nu2[2] + 0.5 * apply((Y.DA - fitted.current)^2, 2, sum)
        nu2 <- 1 / rgamma(J, nu2.posterior.shape, scale=(1/nu2.posterior.scale))

        
        
        ##################
        ## Sample from phi
        ##################
        den.offset <- rho * W.triplet.sum + 1 - rho
        phi.offset <- Y.DA - regression - offset
        Chol.Sigma <- t(chol(proposal.sd.phi*Sigma))
        z.mat <- matrix(rnorm(n=N.all, mean=0, sd=1), nrow=J, ncol=K)
        innovations <- t(Chol.Sigma %*% z.mat)
        temp1 <- gaussianmcarupdateRW(W.triplet, W.begfin, K, J, phi, phi.offset, den.offset, Sigma.inv, rho, nu2, proposal.sd.phi, innovations)      
        phi <- temp1[[1]]
            for(r in 1:J)
            {
            phi[ ,r] <- phi[ ,r] - mean(phi[ ,r])    
            }
        accept[1] <- accept[1] + temp1[[2]]
        accept[2] <- accept[2] + K    
        
        
        
        ####################
        ## Sample from Sigma
        ####################
        Sigma.post.scale <- 2 * prior.Sigma.df * diag(1 / Sigma.a) + t(phi) %*% Q %*% phi
        Sigma <- riwish(Sigma.post.df, Sigma.post.scale)
        Sigma.inv <- solve(Sigma)
        

        
        ######################
        ## Sample from Sigma.a
        ######################
        Sigma.a.posterior.scale <- prior.Sigma.df * diag(Sigma.inv) + 1 / prior.Sigma.scale^2
        Sigma.a <- 1 / rgamma(J, Sigma.a.post.shape, scale=(1/Sigma.a.posterior.scale))   
    
    
    
        ##################
        ## Sample from rho
        ##################
        if(!fix.rho)
        {
            ## Propose a new value
            proposal.rho <- rtruncnorm(n=1, a=0, b=1, mean=rho, sd=proposal.sd.rho)
            Q.prop <- proposal.rho * Wstar + diag(rep(1-proposal.rho), K)
            det.Q.prop <-  sum(log((proposal.rho * Wstar.val + (1-proposal.rho))))    
             
            ## Compute the acceptance rate
            logprob.current <- 0.5 * J * det.Q - 0.5 * sum(diag(t(phi) %*% Q %*% phi %*% Sigma.inv))
            logprob.proposal <- 0.5 * J * det.Q.prop - 0.5 * sum(diag(t(phi) %*% Q.prop %*% phi %*% Sigma.inv))
            hastings <- log(dtruncnorm(x=rho, a=0, b=1, mean=proposal.rho, sd=proposal.sd.rho)) - log(dtruncnorm(x=proposal.rho, a=0, b=1, mean=rho, sd=proposal.sd.rho)) 
            prob <- exp(logprob.proposal - logprob.current + hastings)
            if(prob > runif(1))
            {
                rho <- proposal.rho
                det.Q <- det.Q.prop
                Q <- Q.prop
                accept[3] <- accept[3] + 1           
            }else
            {}              
            accept[4] <- accept[4] + 1       
        }else
        {}
        
        
        
        #########################
        ## Calculate the deviance
        #########################
        fitted <- regression + phi + offset
        loglike <- dnorm(x=as.numeric(t(Y)), mean=as.numeric(t(fitted)), sd=rep(sqrt(nu2), K), log=TRUE)

        
        ###################
        ## Save the results
        ###################
        if(j > burnin & (j-burnin)%%thin==0)
        {
            ele <- (j - burnin) / thin
            samples.beta[ele, ] <- as.numeric(beta)
            samples.nu2[ele, ] <- nu2
            samples.phi[ele, ] <- as.numeric(t(phi))
            samples.Sigma[ele, , ] <- Sigma
            samples.Sigma.a[ele, ] <- Sigma.a
            if(!fix.rho) samples.rho[ele, ] <- rho
            samples.loglike[ele, ] <- loglike
            samples.fitted[ele, ] <- as.numeric(t(fitted))
            if(n.miss>0) samples.Y[ele, ] <- Y.DA[miss.locator]
        }else
        {}
        
        
        
        ########################################
        ## Self tune the acceptance probabilties
        ########################################
        if(ceiling(j/100)==floor(j/100) & j < burnin)
        {
        proposal.sd.phi <- common.accceptrates1(accept[1:2], proposal.sd.phi, 40, 50)
            if(!fix.rho)
            {
                proposal.sd.rho <- common.accceptrates2(accept[3:4], proposal.sd.rho, 40, 50, 0.5)
            }
            accept <- c(0,0,0,0)
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
chain.results <- list(samples.beta=samples.beta, samples.phi=samples.phi, samples.nu2=samples.nu2, samples.Sigma=samples.Sigma, samples.Sigma.a=samples.Sigma.a, samples.rho=samples.rho, samples.loglike=samples.loglike, samples.fitted=samples.fitted,
                      samples.Y=samples.Y, accept=accept)

#### Return the results
return(chain.results)
}            