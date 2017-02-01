poisson.MVlerouxCAR <- function(formula, data=NULL,  W, burnin, n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.Sigma.df=NULL, prior.Sigma.scale=NULL, fix.rho=FALSE, rho=NULL, verbose=TRUE)
{
##############################################
#### Format the arguments and check for errors
##############################################
#### Verbose
a <- common.verbose(verbose)
    
    
#### Frame object
frame.results <- common.frame(formula, data, "poisson")
N.all <- frame.results$n
p <- frame.results$p
X <- frame.results$X
X.standardised <- frame.results$X.standardised
X.sd <- frame.results$X.sd
X.mean <- frame.results$X.mean
X.indicator <- frame.results$X.indicator 
offset <- frame.results$offset
Y <- frame.results$Y
Y.miss <- frame.results$Y.miss
which.miss <- frame.results$which.miss
n.miss <- frame.results$n.miss  

    
#### W matrix
    if(!is.matrix(W)) stop("W is not a matrix.", call.=FALSE)
K <- nrow(W)
    if(ceiling(N.all/K)!= floor(N.all/K)) stop("The number of data points divided by the number of rows in W is not a whole number.", call.=FALSE)
J <- N.all / K


#### Check for errors on rho and fix.rho
    if(!is.logical(fix.rho)) stop("fix.rho is not logical.", call.=FALSE)   
    if(fix.rho & is.null(rho)) stop("rho is fixed but an initial value was not set.", call.=FALSE)   
    if(fix.rho & !is.numeric(rho) ) stop("rho is fixed but is not numeric.", call.=FALSE)  
    if(!fix.rho) rho <- runif(1)
    if(rho<0 ) stop("rho is outside the range [0, 1].", call.=FALSE)  
    if(rho>1 ) stop("rho is outside the range [0, 1].", call.=FALSE)    
    

#### Priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(1000, p)
    if(is.null(prior.Sigma.df)) prior.Sigma.df <- J+1
    if(is.null(prior.Sigma.scale)) prior.Sigma.scale <- diag(rep(1,J))
prior.beta.check(prior.mean.beta, prior.var.beta, p)
prior.varmat.check(prior.Sigma.scale, J)  


#### Compute the blocking structure for beta     
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



#############################
#### Initial parameter values
#############################
mod.glm <- glm(Y~X.standardised-1, offset=offset, family="quasipoisson")
beta.mean <- mod.glm$coefficients
beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
beta <- rnorm(n=p, mean=beta.mean, sd=beta.sd)
    
log.Y <- log(Y)
log.Y[Y==0] <- -0.1  
res.temp <- log.Y - X.standardised %*% beta.mean - offset
res.sd <- sd(res.temp, na.rm=TRUE)/5
phi <- rnorm(n=N.all, mean=0, sd=res.sd)
phi.mat <- matrix(phi, nrow=K, byrow=TRUE)
Sigma <- cov(phi.mat)
Sigma.inv <- solve(Sigma)

    
    
###############################    
#### Set up the MCMC quantities    
###############################
#### Matrices to store samples    
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, p))
samples.phi <- array(NA, c(n.keep, N.all))
samples.Sigma <- array(NA, c(n.keep, J, J))
    if(!fix.rho) samples.rho <- array(NA, c(n.keep, 1))
samples.deviance <- array(NA, c(n.keep, 1))
samples.like <- array(NA, c(n.keep, N.all))
samples.fitted <- array(NA, c(n.keep, N.all))
    if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))

    
#### Metropolis quantities
accept.all <- rep(0,6)
accept <- accept.all
proposal.sd.beta <- 0.01
proposal.sd.phi <- 0.1
proposal.sd.rho <- 0.02
proposal.corr.beta <- solve(t(X.standardised) %*% X.standardised)
chol.proposal.corr.beta <- chol(proposal.corr.beta) 
Sigma.post.df <- prior.Sigma.df + K  
    
    
##################################
#### Set up the spatial quantities
##################################
#### CAR quantities
W.quants <- common.Wcheckformat.leroux(W, K, fix.rho, rho)
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
W.list<- mat2listw(W)
W.nb <- W.list$neighbours
W.islands <- n.comp.nb(W.nb)
islands <- W.islands$comp.id
islands.all <- rep(islands,J)
n.islands <- max(W.islands$nc)
    if(rho==1) Sigma.post.df <- prior.Sigma.df + K - n.islands   

        
#### Make matrix versions of the variables
offset.mat <- matrix(offset, nrow=K, ncol=J, byrow=TRUE) 
regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=J, byrow=TRUE)
Y.mat <- matrix(Y, nrow=K, ncol=J, byrow=TRUE)
Y.mat.miss <- matrix(Y.miss, nrow=K, ncol=J, byrow=TRUE)
which.miss.mat <- matrix(which.miss, nrow=K, ncol=J, byrow=TRUE)
    


###########################
#### Run the Bayesian model
###########################
#### Start timer
    if(verbose)
    {
    cat("Generating", n.keep, "post burnin and thinned (if requested) samples\n", sep = " ")
    progressBar <- txtProgressBar(style = 3)
    percentage.points<-round((1:100/100)*n.sample)
    }else
    {
    percentage.points<-round((1:100/100)*n.sample)     
    }
    
    
#### Create the MCMC samples      
    for(j in 1:n.sample)
    {
    ###################
    ## Sample from beta
    ###################
    offset.temp <- phi + offset
        if(p>2)
        {
        temp <- poissonbetaupdateMALA(X.standardised, N.all, p, beta, offset.temp, Y.miss, prior.mean.beta, prior.var.beta, which.miss, n.beta.block, proposal.sd.beta, list.block)
        }else
        {
        temp <- poissonbetaupdateRW(X.standardised, N.all, p, beta, offset.temp, Y.miss, prior.mean.beta, prior.var.beta, which.miss, proposal.sd.beta)
        }
    beta <- temp[[1]]
    accept[1] <- accept[1] + temp[[2]]
    accept[2] <- accept[2] + n.beta.block  
    regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=J, byrow=TRUE)          

    

    
    
    
    ##################
    ## Sample from phi
    ##################
    den.offset <- rho * W.triplet.sum + 1 - rho
    phi.offset <- regression.mat + offset.mat
    temp1 <- poissonmcarupdate(W.triplet, W.begfin, K, J, phi.mat, Y.mat.miss,  phi.offset, den.offset, Sigma.inv, rho, proposal.sd.phi, which.miss.mat)      
        if(rho<1)
        {
        phi.mat <- temp1[[1]] - mean(temp1[[1]])
        phi <- as.numeric(t(phi.mat))
        }else
        {
        phi.mat <- temp1[[1]]
        phi <- as.numeric(t(phi.mat))
        phi[which(islands.all==1)] <- phi[which(islands.all==1)] - mean(phi[which(islands.all==1)]) 
        phi.mat <- matrix(phi, nrow=K, byrow=TRUE)
        }
    accept[3] <- accept[3] + temp1[[2]]
    accept[4] <- accept[4] + K    

        
        
    ####################
    ## Sample from Sigma
    ####################
    Sigma.post.scale <- t(phi.mat) %*% Q %*% phi.mat + prior.Sigma.scale
    Sigma <- riwish(Sigma.post.df, Sigma.post.scale)
    Sigma.inv <- solve(Sigma)
        
    
        
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
        logprob.current <- 0.5 * J * det.Q - 0.5 * sum(diag(t(phi.mat) %*% Q %*% phi.mat %*% Sigma.inv))
        logprob.proposal <- 0.5 * J * det.Q.prop - 0.5 * sum(diag(t(phi.mat) %*% Q.prop %*% phi.mat %*% Sigma.inv))
        prob <- exp(logprob.proposal - logprob.current)
            if(prob > runif(1))
            {
            rho <- proposal.rho
            det.Q <- det.Q.prop
            Q <- Q.prop
            accept[5] <- accept[5] + 1           
            }else
            {}              
        accept[6] <- accept[6] + 1       
        }else
        {}
        
        
        
    #########################
    ## Calculate the deviance
    #########################
    lp <- as.numeric(X.standardised %*% beta) + phi + offset
    fitted <- exp(lp)
    deviance.all <- dpois(x=as.numeric(Y), lambda=fitted, log=TRUE)
    like <- exp(deviance.all)
    deviance <- -2 * sum(deviance.all, na.rm=TRUE)  
        
        
    ###################
    ## Save the results
    ###################
        if(j > burnin & (j-burnin)%%thin==0)
        {
        ele <- (j - burnin) / thin
        samples.beta[ele, ] <- beta
        samples.phi[ele, ] <- phi
        samples.Sigma[ele, , ] <- Sigma
            if(!fix.rho) samples.rho[ele, ] <- rho
        samples.deviance[ele, ] <- deviance
        samples.like[ele, ] <- like
        samples.fitted[ele, ] <- fitted
            if(n.miss>0) samples.Y[ele, ] <- rpois(n=n.miss, lambda=fitted[which.miss==0])
        }else
        {}
        
        
    
    ########################################
    ## Self tune the acceptance probabilties
    ########################################
    k <- j/100
        if(ceiling(k)==floor(k))
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
        accept.all <- accept.all + accept
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
    cat("\nSummarising results")
    close(progressBar)
    }else
    {}



###################################
#### Summarise and save the results 
###################################
#### Compute the acceptance rates
accept.beta <- 100 * accept.all[1] / accept.all[2]
accept.phi <- 100 * accept.all[3] / accept.all[4]
    if(!fix.rho)
    {
    accept.rho <- 100 * accept.all[5] / accept.all[6]
    }else
    {
    accept.rho <- NA    
    }
accept.Sigma <- 100
accept.final <- c(accept.beta, accept.phi, accept.rho, accept.Sigma)
names(accept.final) <- c("beta", "phi", "rho", "Sigma")
    

#### Deviance information criterion (DIC)
median.beta <- apply(samples.beta, 2, median)
median.phi <- apply(samples.phi, 2, median)
fitted.median <- exp(X.standardised %*% median.beta + median.phi + offset)
deviance.fitted <- -2 * sum(dpois(x=Y, lambda=fitted.median, log=TRUE), na.rm=TRUE)
p.d <- median(samples.deviance) - deviance.fitted
DIC <- 2 * median(samples.deviance) - deviance.fitted     
    
    
#### Watanabe-Akaike Information Criterion (WAIC)
LPPD <- sum(log(apply(samples.like,2,mean)), na.rm=TRUE)
p.w <- sum(apply(log(samples.like),2,var), na.rm=TRUE)
WAIC <- -2 * (LPPD - p.w)
    
#### Compute the Conditional Predictive Ordinate
CPO <- rep(NA, N.all)
    for(j in 1:N.all)
    {
    CPO[j] <- 1/median((1 / dpois(x=Y[j], lambda=samples.fitted[ ,j])))    
    }
LMPL <- sum(log(CPO), na.rm=TRUE)  
    
    
#### transform the parameters back to the origianl covariate scale.
samples.beta.orig <- common.betatransform(samples.beta, X.indicator, X.mean, X.sd, p, FALSE)


#### Create a summary object
samples.beta.orig <- mcmc(samples.beta.orig)
summary.beta <- t(apply(samples.beta.orig, 2, quantile, c(0.5, 0.025, 0.975))) 
summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta,p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
rownames(summary.beta) <- colnames(X)
colnames(summary.beta) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
    
summary.hyper <- array(NA, c((J+1) ,7))
summary.hyper[1:J, 1] <- diag(apply(samples.Sigma, c(2,3), quantile, c(0.5)))
summary.hyper[1:J, 2] <- diag(apply(samples.Sigma, c(2,3), quantile, c(0.025)))
summary.hyper[1:J, 3] <- diag(apply(samples.Sigma, c(2,3), quantile, c(0.975)))
summary.hyper[1:J, 4] <- n.keep
summary.hyper[1:J, 5] <- accept.Sigma
summary.hyper[1:J, 6] <- diag(apply(samples.Sigma, c(2,3), effectiveSize))
    for(r in 1:J)
    {
    summary.hyper[r, 7] <- geweke.diag(samples.Sigma[ ,r,r])$z    
    }
    
    if(!fix.rho)
    {
    summary.hyper[(J+1), 1:3] <- quantile(samples.rho, c(0.5, 0.025, 0.975))
    summary.hyper[(J+1), 4:7] <- c(n.keep, accept.rho, effectiveSize(samples.rho), geweke.diag(samples.rho)$z)
    }else
    {
    summary.hyper[(J+1), 1:3] <- c(rho, rho, rho)
    summary.hyper[(J+1), 4:7] <- rep(NA, 4)
    }
    
summary.results <- rbind(summary.beta, summary.hyper)
rownames(summary.results)[(p+1): nrow(summary.results)] <- c(paste(rep("Sigma",J), 1:J, 1:J, sep=""), "rho")
summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
    
    
#### Create the fitted values and residuals
fitted.values <- apply(samples.fitted, 2, median)
response.residuals <- as.numeric(Y) - fitted.values
pearson.residuals <- response.residuals /sqrt(fitted.values)
deviance.residuals <- sign(response.residuals) * sqrt(2 * (Y * log(Y/fitted.values) + fitted.values - Y))
residuals <- data.frame(response=response.residuals, pearson=pearson.residuals, deviance=deviance.residuals)
    
    
#### Compile and return the results
loglike <- (-0.5 * deviance.fitted)
modelfit <- c(DIC, p.d, WAIC, p.w, LMPL, loglike)
names(modelfit) <- c("DIC", "p.d", "WAIC", "p.w", "LMPL", "loglikelihood")
model.string <- c("Likelihood model - Poisson (log link function)", "\nRandom effects model - Leroux MCAR\n")
    if(fix.rho) samples.rho=NA
    if(n.miss==0) samples.Y = NA
samples <- list(beta=samples.beta.orig, phi=mcmc(samples.phi), Sigma=samples.Sigma, rho=mcmc(samples.rho), fitted=mcmc(samples.fitted), Y=mcmc(samples.Y))
results <- list(summary.results=summary.results, samples=samples, fitted.values=fitted.values, residuals=residuals, modelfit=modelfit, accept=accept.final, localised.structure=NULL,  formula=formula, model=model.string, X=X)
class(results) <- "CARBayes"
    

#### Finish by stating the time taken    
    if(verbose)
    {
    b<-proc.time()
    cat(" finished in ", round(b[3]-a[3], 1), "seconds")
    }else
    {}
return(results)
}
