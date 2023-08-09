multinomial.MVlerouxCARMCMC <- function(Y, trials, offset, X.standardised, W, rho, fix.rho, K, p, J, N.all, N.re, which.miss, n.miss, burnin, n.sample, thin, n.beta.block, list.block, prior.mean.beta, prior.var.beta, prior.Sigma.df, prior.Sigma.scale, verbose, chain)
{
    # Rcpp::sourceCpp("src/CARBayes.cpp")   
    # source("R/common.functions.R")
    # library(spdep)
    # library(truncnorm)  
    # library(MCMCpack)
    ##########################################
    #### Generate the initial parameter values
    ##########################################
    #### Generate initial values for each chain
    beta <- array(NA, c(p, (J-1)))
    for(i in 2:J)
    {
        mod.glm <- glm(cbind(Y[ ,i], trials - Y[ ,i])~X.standardised-1, offset=offset[ ,(i-1)], family="quasibinomial")
        beta.mean <- mod.glm$coefficients
        beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
        beta[ ,(i-1)] <- rnorm(n=p, mean=beta.mean, sd=beta.sd)
    }
    regression <- X.standardised %*% beta
    
    theta.hat <- Y / trials
    theta.hat[theta.hat==0] <- 0.01
    theta.hat[theta.hat==1] <- 0.99
    res.temp <- log(theta.hat[ ,-1]  / theta.hat[ ,1]) - offset - regression
    res.sd <- sd(res.temp, na.rm=TRUE)/5
    phi.vec <- rnorm(n=N.re, mean=0, sd=res.sd)
    phi <- matrix(phi.vec, nrow=K, byrow=TRUE)
    Sigma <- cov(phi)
    Sigma.inv <- solve(Sigma)    
    Sigma.a <- rep(1, (J-1))

    
    
    ####################################################################
    #### Compute the fitted values based on the current parameter values
    ####################################################################   
    regression <- X.standardised %*% beta
    Y.DA <- Y
    
    
    
    #### If only one element in Y is missing then fix it as we know the total number of trials
    which.miss.row <- J-apply(which.miss,1,sum)
    which.miss.1 <- which(which.miss.row==1)
    if(length(length(which.miss.1))>0)
    {
        for(r in 1:length(which.miss.1))
        {
            which.miss[which.miss.1[r], is.na(Y[which.miss.1[r], ])] <- 1
            Y[which.miss.1[r], is.na(Y[which.miss.1[r], ])] <- trials[which.miss.1[r]] - sum(Y[which.miss.1[r], ], na.rm=T)    
        }
        n.miss <- sum(is.na(Y))
        which.miss.row <- J-apply(which.miss,1,sum)
    }else
    {}
    const.like <- lfactorial(trials[which.miss.row==0]) - apply(lfactorial(Y[which.miss.row==0, ]),1,sum)
    K.present <- sum(which.miss.row==0)
    
    #### Determine which rows have missing values
    if(n.miss>0)    which.miss.row2 <- which(which.miss.row>0)   
    
    
    
    #### Matrices to store samples    
    n.keep <- floor((n.sample - burnin)/thin)
    samples.beta <- array(NA, c(n.keep, (J-1)*p))
    samples.phi <- array(NA, c(n.keep, N.re))
    samples.Sigma <- array(NA, c(n.keep, (J-1), (J-1)))
    samples.Sigma.a <- array(NA, c(n.keep, (J-1)))
    if(!fix.rho) samples.rho <- array(NA, c(n.keep, 1))
    samples.loglike <- array(NA, c(n.keep, K.present))
    samples.fitted <- array(NA, c(n.keep, N.all))
    if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))
    
    
    #### Metropolis quantities
    accept.beta <- rep(0,2*(J-1))
    proposal.sd.beta <- rep(0.01, (J-1))
    accept <- rep(0,4)
    proposal.sd.phi <- 0.1
    proposal.sd.rho <- 0.02
    Sigma.post.df <- prior.Sigma.df + K  + J - 2
    Sigma.a.post.shape <- (prior.Sigma.df + J-1) / 2
    
    
    
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
    if(rho==1) Sigma.post.df <- prior.Sigma.df + K  + J - 2 - n.islands   

    
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
            for(g in 1:length(which.miss.row2))   
            {
                ## Determine which row (area) of Y to update
                row <- which.miss.row2[g]
                
                ## Compute the vector of probabilities for that row
                lp <- c(0, regression[row, ] + phi[row, ] + offset[row, ])
                prob <- exp(lp)  / sum(exp(lp))
                
                ## Do the multinomial data augmentation
                if(which.miss.row[row]==J)
                {
                    ## All the Ys are missing
                    Y.DA[row, ] <- as.numeric(rmultinom(n=1, size=trials[row], prob=prob))    
                }else
                {
                    ## Not all the Ys are missing
                    ## Re-normalise the probabilities
                    prob[!is.na(Y[row, ])] <- 0
                    prob <- prob / sum(prob)
                    temp <- as.numeric(rmultinom(n=1, size=trials[row]-sum(Y[row, ], na.rm=T), prob=prob))    
                    Y.DA[row, which.miss[row, ]==0]  <- temp[which.miss[row, ]==0]  
                }
            }
        }else
        {}
        
        
        
        ###################
        ## Sample from beta
        ###################
        offset.temp <- phi + offset
        for(r in 1:(J-1))
        {
            temp <- multinomialbetaupdateRW(X.standardised, K, J, p, r, beta, offset.temp, Y.DA, prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta[r], list.block, rep(0, K))
            beta[ ,r] <- temp[[1]][ ,r]
            accept.beta[r] <- accept.beta[r] + temp[[2]]
            accept.beta[(r+J-1)] <- accept.beta[(r+J-1)] + n.beta.block  
        }
        regression <- X.standardised %*% beta       
        
        
        
        ##################    
        ## Sample from phi
        ##################
        den.offset <- rho * W.triplet.sum + 1 - rho
        phi.offset <- regression + offset
        Chol.Sigma <- t(chol(proposal.sd.phi*Sigma))
        z.mat <- matrix(rnorm(n=N.all-K, mean=0, sd=1), nrow=J-1, ncol=K)
        innovations <- t(Chol.Sigma %*% z.mat)
        temp1 <- multinomialmcarupdateRW(W.triplet, W.begfin, K, J, phi, Y.DA, phi.offset, den.offset, Sigma.inv, rho,  proposal.sd.phi, innovations)      
        phi <- temp1[[1]]
        for(r in 1:(J-1))
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
        Sigma.a <- 1 / rgamma((J-1), Sigma.a.post.shape, scale=(1/Sigma.a.posterior.scale))   
        
        
        
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
            logprob.current <- 0.5 * (J-1) * det.Q - 0.5 * sum(diag(t(phi) %*% Q %*% phi %*% Sigma.inv))
            logprob.proposal <- 0.5 * (J-1) * det.Q.prop - 0.5 * sum(diag(t(phi) %*% Q.prop %*% phi %*% Sigma.inv))
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
        lp <- regression + phi + offset
        lp <- cbind(rep(0,K), lp)
        prob <- exp(lp)  / apply(exp(lp),1,sum)
        fitted <- prob * trials
        loglike <-  const.like + apply(Y[which.miss.row==0, ] * log(prob[which.miss.row==0, ]),1,sum)
        
        
        
        ###################
        ## Save the results
        ###################
        if(j > burnin & (j-burnin)%%thin==0)
        {
            ele <- (j - burnin) / thin
            samples.beta[ele, ] <- as.numeric(beta)
            samples.phi[ele, ] <- as.numeric(t(phi))
            samples.Sigma[ele, , ] <- Sigma
            samples.Sigma.a[ele, ] <- Sigma.a
            if(!fix.rho) samples.rho[ele, ] <- rho
            samples.loglike[ele, ] <- loglike
            samples.fitted[ele, ] <- as.numeric(t(fitted))
            if(n.miss>0) samples.Y[ele, ] <- t(Y.DA)[is.na(t(Y))]
        }else
        {}
        
        
        
        ########################################
        ## Self tune the acceptance probabilties
        ########################################
        if(ceiling(j/100)==floor(j/100) & j < burnin)
        {
            #### Update the proposal sds
            for(r in 1:(J-1))
            {
                if(p>2)
                {
                    proposal.sd.beta[r] <- common.accceptrates1(accept.beta[c(r, (r+J-1))], proposal.sd.beta[r], 40, 50)
                }else
                {
                    proposal.sd.beta[r] <- common.accceptrates1(accept.beta[c(r, (r+J-1))], proposal.sd.beta[r], 30, 40)    
                }
            }
            
            proposal.sd.phi <- common.accceptrates1(accept[1:2], proposal.sd.phi, 40, 50)
            if(!fix.rho)
            {
                proposal.sd.rho <- common.accceptrates2(accept[3:4], proposal.sd.rho, 40, 50, 0.5)
            }
            accept <- c(0,0,0,0)
            accept.beta <- rep(0,2*(J-1))
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
    
    
    #### Close the progress bar if used
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
    chain.results <- list(samples.beta=samples.beta, samples.phi=samples.phi, samples.Sigma=samples.Sigma, samples.Sigma.a=samples.Sigma.a, samples.rho=samples.rho, samples.loglike=samples.loglike, samples.fitted=samples.fitted,
                          samples.Y=samples.Y, accept=accept, accept.beta=accept.beta)
    
    #### Return the results
    return(chain.results)
}            