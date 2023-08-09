multinomial.glmMCMC <- function(Y, trials, offset, X.standardised, K, p, J, N.all, which.miss, n.miss, burnin, n.sample, thin, n.beta.block, list.block, prior.mean.beta, prior.var.beta, verbose, chain)
{
    # Rcpp::sourceCpp("src/CARBayes.cpp")   
    # source("R/common.functions.R")
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

    
    
    ########################################    
    #### Set up the MCMC model run quantities    
    #########################################
    #### Matrices to store samples    
    n.keep <- floor((n.sample - burnin)/thin)
    samples.beta <- array(NA, c(n.keep, (J-1)*p))
    samples.loglike <- array(NA, c(n.keep, K.present))
    samples.fitted <- array(NA, c(n.keep, N.all))
    if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))
    
    
    #### Metropolis quantities
    accept.beta <- rep(0,2*(J-1))
    proposal.sd.beta <- rep(0.01, (J-1))

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
                lp <- c(0, regression[row, ] + offset[row, ])
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
        for(r in 1:(J-1))
        {
            temp <- multinomialbetaupdateRW(X.standardised, K, J, p, r, beta, offset, Y.DA, prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta[r], list.block, rep(0, K))
            beta[ ,r] <- temp[[1]][ ,r]
            accept.beta[r] <- accept.beta[r] + temp[[2]]
            accept.beta[(r+J-1)] <- accept.beta[(r+J-1)] + n.beta.block  
        }
        regression <- X.standardised %*% beta       
        
        
        
        #########################
        ## Calculate the deviance
        #########################
        lp <- regression + offset
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
    chain.results <- list(samples.beta=samples.beta, samples.loglike=samples.loglike, samples.fitted=samples.fitted,
                          samples.Y=samples.Y, accept.beta=accept.beta)

    #### Return the results
    return(chain.results)
}