multinomial.glm <- function(formula, data=NULL, trials, burnin, n.sample, thin=1, n.chains=1, n.cores=1, prior.mean.beta=NULL, prior.var.beta=NULL, verbose=TRUE)
{
##############################################
#### Format the arguments and check for errors
##############################################
#### Verbose
a <- common.verbose(verbose)
    
    
#### Frame object
frame.results <- common.frame(formula, data, "multinomial")
K <- frame.results$n
p <- frame.results$p
X <- frame.results$X
X.standardised <- frame.results$X.standardised
X.sd <- frame.results$X.sd
X.mean <- frame.results$X.mean
X.indicator <- frame.results$X.indicator 
offset <- frame.results$offset
Y <- frame.results$Y
which.miss <- frame.results$which.miss
n.miss <- frame.results$n.miss  
J <- ncol(Y)
N.all <- K * J


#### Check and format the trials argument
    if(sum(is.na(trials))>0) stop("the numbers of trials has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(trials)) stop("the numbers of trials has non-numeric values.", call.=FALSE)
int.check <- K-sum(ceiling(trials)==floor(trials))
    if(int.check > 0) stop("the numbers of trials has non-integer values.", call.=FALSE)
    if(min(trials)<=0) stop("the numbers of trials has zero or negative values.", call.=FALSE)
diffs <- apply(Y, 1, sum, na.rm=T) - trials
    if(max(diffs)>0) stop("the response variable has larger values that the numbers of trials.", call.=FALSE)


    
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
if(n.miss>0)    which.miss.row2 <- which(which.miss.row>0)   


#### Priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(100000, p)
common.prior.beta.check(prior.mean.beta, prior.var.beta, p)

    
#### Compute the blocking structure for beta   
block.temp <- common.betablock(p, 5)
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
    
    

########################
#### Run the MCMC chains
########################
if(n.chains==1)
{
    #### Only 1 chain
    results <- multinomial.glmMCMC(Y=Y, trials=trials, offset=offset, X.standardised=X.standardised, K=K, p=p, J=J, N.all=N.all, which.miss=which.miss, n.miss=n.miss, burnin=burnin, n.sample=n.sample, thin=thin, n.beta.block=n.beta.block, list.block=list.block, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, verbose=verbose, chain=1)
}else if(n.chains > 1 & ceiling(n.chains)==floor(n.chains) & n.cores==1)
{
    #### Multiple chains in  series
    results <- as.list(rep(NA, n.chains))
    for(i in 1:n.chains)
    {
        results[[i]] <- multinomial.glmMCMC(Y=Y, trials=trials, offset=offset, X.standardised=X.standardised, K=K, p=p, J=J, N.all=N.all, which.miss=which.miss, n.miss=n.miss, burnin=burnin, n.sample=n.sample, thin=thin, n.beta.block=n.beta.block, list.block=list.block, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, verbose=verbose, chain=i)
    }
}else if(n.chains > 1 & ceiling(n.chains)==floor(n.chains) & n.cores>1 & ceiling(n.cores)==floor(n.cores))   
{   
    #### Multiple chains in parallel
    results <- as.list(rep(NA, n.chains))
    if(verbose)
    {
        compclust <- makeCluster(n.cores, outfile="CARBayesprogress.txt")
        cat("The current progress of the model fitting algorithm has been output to CARBayesprogress.txt in the working directory")
    }else
    {
        compclust <- makeCluster(n.cores)
    }
    results <- clusterCall(compclust, fun=multinomial.glmMCMC, Y=Y, trials=trials, offset=offset, X.standardised=X.standardised, K=K, p=p, J=J, N.all=N.all, which.miss=which.miss, n.miss=n.miss, burnin=burnin, n.sample=n.sample, thin=thin, n.beta.block=n.beta.block, list.block=list.block, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, verbose=verbose, chain="all")
    stopCluster(compclust)
}else
{
    stop("n.chains or n.cores are not positive integers.", call.=FALSE)  
}


#### end timer
if(verbose)
{
    cat("\nSummarising results.\n")
}else
{}



###################################
#### Summarise and save the results 
###################################
if(n.chains==1)
{
    ## Compute the acceptance rates
    accept.beta <- 100 * sum(results$accept.beta[1:(J-1)]) / sum(results$accept.beta[(J:(2*(J-1)))])
    accept.final <- accept.beta
    names(accept.final) <- c("beta")

    ## Compute the model fit criterion
    mean.beta <- matrix(apply(results$samples.beta, 2, mean), nrow=p, ncol=(J-1), byrow=F)
    mean.logit <- X.standardised %*% mean.beta + offset
    mean.logit <- cbind(rep(0,K), mean.logit)
    mean.prob <- exp(mean.logit)  / apply(exp(mean.logit),1,sum)
    deviance.fitted <- -2* sum(const.like + apply(Y[which.miss.row==0, ] * log(mean.prob[which.miss.row==0, ]),1,sum))
    modelfit <- common.modelfit(results$samples.loglike, deviance.fitted)

    ## Create the Fitted values and residuals
    fitted.values <- matrix(apply(results$samples.fitted, 2, mean), nrow=K, ncol=J, byrow=T)
    response.residuals <- Y - fitted.values
    var.y <- fitted.values * (1-fitted.values /  trials)
    pearson.residuals <- response.residuals / sqrt(var.y)
    residuals <- list(response=response.residuals, pearson=pearson.residuals)
    
    ## Create MCMC objects and back transform the regression parameters
    samples.beta.orig <- results$samples.beta
        for(r in 1:(J-1))
        {
        samples.beta.orig[ ,((r-1)*p+1):(r*p)] <- common.betatransform(results$samples.beta[ ,((r-1)*p+1):(r*p)], X.indicator, X.mean, X.sd, p, FALSE)
        }
     samples <- list(beta=mcmc(samples.beta.orig), fitted=mcmc(results$samples.fitted), Y=mcmc(results$samples.Y))
    
    #### Create a summary object
    n.keep <- floor((n.sample - burnin)/thin)
    summary.beta <- t(rbind(apply(samples$beta, 2, mean), apply(samples$beta, 2, quantile, c(0.025, 0.975)))) 
    summary.beta <- cbind(summary.beta, rep(n.keep, (J-1)*p), rep(accept.beta,(J-1)*p), effectiveSize(samples$beta), geweke.diag(samples$beta)$z)
    col.name <- rep(NA, p*(J-1))
    
    if(is.null(colnames(Y)))
    {
        for(r in 1:(J-1))
        {
            col.name[((r-1)*p+1):(r*p)] <- paste("Category ", r+1,  " - ", colnames(X), sep="")   
        }
    }else
    {
        for(r in 1:(J-1))
        {
            col.name[((r-1)*p+1):(r*p)] <- paste(colnames(Y)[(r+1)],  " - ", colnames(X), sep="")   
        }
    }
    
    rownames(summary.beta) <- col.name
    colnames(summary.beta) <- c("Mean", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
    summary.results <- summary.beta
    summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
    summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
}else
{
    ## Compute the acceptance rates
    accept.temp <- lapply(results, function(l) l[["accept.beta"]])
    accept.temp2 <- do.call(what=rbind, args=accept.temp)
    accept.beta <- 100 * sum(accept.temp2[ ,1:(J-1)]) / sum(accept.temp2[ ,(J:(2*(J-1)))])
    accept.final <- c(accept.beta)
    names(accept.final) <- c("beta")

    ## Extract the samples into separate lists
    samples.beta.list <- lapply(results, function(l) l[["samples.beta"]])
    samples.loglike.list <- lapply(results, function(l) l[["samples.loglike"]])
    samples.fitted.list <- lapply(results, function(l) l[["samples.fitted"]])
    samples.Y.list <- lapply(results, function(l) l[["samples.Y"]])
    
    ## Convert the samples into separate matrix objects
    samples.beta.matrix <- do.call(what=rbind, args=samples.beta.list)
    samples.loglike.matrix <- do.call(what=rbind, args=samples.loglike.list)
    samples.fitted.matrix <- do.call(what=rbind, args=samples.fitted.list)
    
    ## Compute the model fit criteria
    mean.beta <- matrix(apply(samples.beta.matrix, 2, mean), nrow=p, ncol=(J-1), byrow=F)
    mean.logit <- X.standardised %*% mean.beta + offset
    mean.logit <- cbind(rep(0,K), mean.logit)
    mean.prob <- exp(mean.logit)  / apply(exp(mean.logit),1,sum)
    deviance.fitted <- -2* sum(const.like + apply(Y[which.miss.row==0, ] * log(mean.prob[which.miss.row==0, ]),1,sum))
    modelfit <- common.modelfit(samples.loglike.matrix, deviance.fitted)

    ## Create the Fitted values and residuals
    fitted.values <- matrix(apply(samples.fitted.matrix, 2, mean), nrow=K, ncol=J, byrow=T)
    response.residuals <- Y - fitted.values
    var.y <- fitted.values * (1-fitted.values /  trials)
    pearson.residuals <- response.residuals / sqrt(var.y)
    residuals <- list(response=response.residuals, pearson=pearson.residuals)

    ## Backtransform the regression parameters
    samples.beta.list <- samples.beta.list
    for(j in 1:n.chains)
    {
        for(r in 1:(J-1))
        {
        samples.beta.list[[j]][ ,((r-1)*p+1):(r*p)] <- common.betatransform(samples.beta.list[[j]][ ,((r-1)*p+1):(r*p)], X.indicator, X.mean, X.sd, p, FALSE)  
        }
    }
    samples.beta.matrix <- do.call(what=rbind, args=samples.beta.list)

    ## Create MCMC objects
    beta.temp <- samples.beta.list
    loglike.temp <- samples.loglike.list
    fitted.temp <- samples.fitted.list
    Y.temp <- samples.Y.list
    for(j in 1:n.chains)
    {
        beta.temp[[j]] <- mcmc(samples.beta.list[[j]])   
        loglike.temp[[j]] <- mcmc(samples.loglike.list[[j]])   
        fitted.temp[[j]] <- mcmc(samples.fitted.list[[j]])   
        Y.temp[[j]] <- mcmc(samples.Y.list[[j]])   
    }
    beta.mcmc <- as.mcmc.list(beta.temp)
    fitted.mcmc <- as.mcmc.list(fitted.temp)
    Y.mcmc <- as.mcmc.list(Y.temp)
    samples <- list(beta=beta.mcmc, fitted=fitted.mcmc, Y=Y.mcmc)
    
    ## Create a summary object
    n.keep <- floor((n.sample - burnin)/thin)
    summary.beta <- t(rbind(apply(samples.beta.matrix, 2, mean), apply(samples.beta.matrix, 2, quantile, c(0.025, 0.975)))) 
    summary.beta <- cbind(summary.beta, rep(n.keep, (J-1)*p), rep(accept.beta,(J-1)*p), effectiveSize(beta.mcmc), gelman.diag(beta.mcmc)$psrf[ ,2])
    col.name <- rep(NA, p*(J-1))
    
    if(is.null(colnames(Y)))
    {
        for(r in 1:(J-1))
        {
            col.name[((r-1)*p+1):(r*p)] <- paste("Category ", r+1,  " - ", colnames(X), sep="")   
        }
    }else
    {
        for(r in 1:(J-1))
        {
            col.name[((r-1)*p+1):(r*p)] <- paste(colnames(Y)[(r+1)],  " - ", colnames(X), sep="")   
        }
    }
    
    rownames(summary.beta) <- col.name
    colnames(summary.beta) <- c("Mean", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "PSRF (upper 95% CI)")
    summary.results <- summary.beta
    summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
    summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
}



###################################
#### Compile and return the results
###################################
model.string <- c("Likelihood model - Multinomial (logit link function)", "\nRandom effects model - None\n")
n.total <- floor((n.sample - burnin) / thin) * n.chains
mcmc.info <- c(n.total, n.sample, burnin, thin, n.chains)
names(mcmc.info) <- c("Total samples", "n.sample", "burnin", "thin", "n.chains")
results <- list(summary.results=summary.results, samples=samples, fitted.values=fitted.values, residuals=residuals, modelfit=modelfit, accept=accept.final, localised.structure=NULL,  formula=formula, model=model.string, mcmc.info=mcmc.info, X=X)
class(results) <- "CARBayes"
if(verbose)
{
    b<-proc.time()
    cat("Finished in ", round(b[3]-a[3], 1), "seconds.\n")
}else
{}
return(results)
}




