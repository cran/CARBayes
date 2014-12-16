gaussian.clusterCAR <- function(formula, data=NULL, G, W, burnin=0, n.sample=1000, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.tau2=NULL, prior.nu2=NULL, prior.delta = NULL, verbose=TRUE)
{
#### Check on the verbose option
     if(is.null(verbose)) verbose=TRUE     
     if(!is.logical(verbose)) stop("the verbose option is not logical.", call.=FALSE)

     if(verbose)
     {
     cat("Setting up the model\n")
     a<-proc.time()
     }else{}

##############################################
#### Format the arguments and check for errors
##############################################
#### Overall formula object
frame <- try(suppressWarnings(model.frame(formula, data=data, na.action=na.pass)), silent=TRUE)
     if(class(frame)=="try-error") stop("the formula inputted contains an error, e.g the variables may be different lengths.", call.=FALSE)
                   
 

#### Design matrix
## Create the matrix
X <- try(suppressWarnings(model.matrix(object=attr(frame, "terms"), data=frame)), silent=TRUE)
     if(class(X)=="try-error") stop("the covariate matrix contains inappropriate values.", call.=FALSE)
     if(sum(is.na(X))>0) stop("the covariate matrix contains missing 'NA' values.", call.=FALSE)

n <- nrow(X)
p <- ncol(X)
   

## Check for linearly related columns
cor.X <- suppressWarnings(cor(X))
diag(cor.X) <- 0
          
     if(max(cor.X, na.rm=TRUE)==1) stop("the covariate matrix has two exactly linearly related columns.", call.=FALSE)
     if(min(cor.X, na.rm=TRUE)==-1) stop("the covariate matrix has two exactly linearly related columns.", call.=FALSE)
          
     if(p>1)
     {
          if(sort(apply(X, 2, sd))[2]==0) stop("the covariate matrix has two intercept terms.", call.=FALSE)
     }else
     {
     }
          
     ## Standardise the matrix
     X.standardised <- X
     X.sd <- apply(X, 2, sd)
     X.mean <- apply(X, 2, mean)
     X.indicator <- rep(NA, p)       # To determine which parameter estimates to transform back
          
     for(j in 1:p)
     {
          if(length(table(X[ ,j]))>2)
          {
          X.indicator[j] <- 1
          X.standardised[ ,j] <- (X[ ,j] - mean(X[ ,j])) / sd(X[ ,j])
          }else if(length(table(X[ ,j]))==1)
          {
          X.indicator[j] <- 2
          }else
          {
          X.indicator[j] <- 0
          }
     }
          
          
          
#### Response variable
## Create the response
Y <- model.response(frame)
n <- length(Y)

    if(sum(is.na(Y))>0) stop("the response has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(Y)) stop("the response variable has non-numeric values.", call.=FALSE)


     
#### Offset variable
## Create the offset
offset <- try(model.offset(frame), silent=TRUE)
          
## Check for errors
     if(class(offset)=="try-error")   stop("the offset is not numeric.", call.=FALSE)
     if(is.null(offset))  offset <- rep(0,n)
     if(sum(is.na(offset))>0) stop("the offset has missing 'NA' values.", call.=FALSE)
     if(!is.numeric(offset)) stop("the offset variable has non-numeric values.", call.=FALSE)
 
     
#### Number of clusters G
  if(length(G)!=1) stop("G is the wrong length.", call.=FALSE)    
  if(!is.numeric(G)) stop("G is not numeric.", call.=FALSE)    
  if(G<=0) stop("G is not positive.", call.=FALSE)    
  if(G!=round(G)) stop("G is not an integer.", call.=FALSE) 
Gbar <- (G+1)/2 
  
     
#### Priors
     if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
     if(is.null(prior.var.beta)) prior.var.beta <- rep(1000, p)
     if(is.null(prior.tau2)) prior.tau2 <- c(0.001, 0.001)
     if(is.null(prior.nu2)) prior.nu2 <- c(0.001, 0.001)
     if(is.null(prior.delta)) prior.delta <- 10

     
     if(length(prior.mean.beta)!=p) stop("the vector of prior means for beta is the wrong length.", call.=FALSE)    
     if(!is.numeric(prior.mean.beta)) stop("the vector of prior means for beta is not numeric.", call.=FALSE)    
     if(sum(is.na(prior.mean.beta))!=0) stop("the vector of prior means for beta has missing values.", call.=FALSE)    
          
     if(length(prior.var.beta)!=p) stop("the vector of prior variances for beta is the wrong length.", call.=FALSE)    
     if(!is.numeric(prior.var.beta)) stop("the vector of prior variances for beta is not numeric.", call.=FALSE)    
     if(sum(is.na(prior.var.beta))!=0) stop("the vector of prior variances for beta has missing values.", call.=FALSE)    
     if(min(prior.var.beta) <=0) stop("the vector of prior variances has elements less than zero", call.=FALSE)

     if(length(prior.tau2)!=2) stop("the prior value for tau2 is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.tau2)) stop("the prior value for tau2 is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.tau2))!=0) stop("the prior value for tau2 has missing values.", call.=FALSE)    
          
    if(length(prior.nu2)!=2) stop("the prior value for nu2 is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.nu2)) stop("the prior value for nu2 is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.nu2))!=0) stop("the prior value for nu2 has missing values.", call.=FALSE)    

    if(length(prior.delta)!=1) stop("the prior value for delta is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.delta)) stop("the prior value for delta is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.delta))!=0) stop("the prior value for delta has missing values.", call.=FALSE)    
    if(prior.delta<=0) stop("the prior value for delta is not positive.", call.=FALSE)    

     
     
     
#### Initial parameter values
## Cluster allocation
mod.glm <- glm(Y~X.standardised-1 + offset(offset), family=gaussian)
res <- Y - offset - (X.standardised %*% mod.glm$coefficients) + mod.glm$coefficients[which(X.indicator==2)]
cluster.factor <- kmeans(x=res, centers=G, nstart=1000)$cluster
Z.temp <- array(NA, c(n,G))
          for(i in 1:G)
          {
          Z.temp[ ,i] <- as.numeric(cluster.factor==i)    
          }
classes.means <- tapply(res, cluster.factor, mean)     
Z <- cluster.factor
     for(i in 1:G)
     {
     Z[cluster.factor==order(classes.means)[i]] <-  i    
     }
     
lambda <- tapply(res, Z, mean)
delta <- runif(1,0, prior.delta)
beta <- mod.glm$coefficients
beta[which(X.indicator==2)] <- 0    
phi <- rnorm(n=n, mean=rep(0,n), sd=rep(0.1, n))
     for(i in 1:G)
     {
     phi[which(Z==i)] <- phi[which(Z==i)] - mean(phi[which(Z==i)])
     }
tau2 <- runif(1)
nu2 <- runif(1)
rho <- runif(1)  
beta.regression <- as.numeric(X.standardised %*% beta)   
     
         
#### MCMC quantities
## Checks
    if(!is.numeric(burnin)) stop("burn-in is not a number", call.=FALSE)
    if(!is.numeric(n.sample)) stop("n.sample is not a number", call.=FALSE)    
    if(n.sample <= 0) stop("n.sample is less than or equal to zero.", call.=FALSE)
    if(burnin < 0) stop("burn-in is less than zero.", call.=FALSE)
    if(n.sample <= burnin)  stop("Burn-in is greater than n.sample.", call.=FALSE)
    if(!is.numeric(thin)) stop("thin is not a number", call.=FALSE)
    if(thin <= 0) stop("thin is less than or equal to zero.", call.=FALSE)

## Compute the blocking structure for beta     
blocksize <- 5
     if(blocksize >= p)
     {
     n.beta.block <- 1
     beta.beg <- 1
     beta.fin <- p
     }else
     {
     n.standard <- 1 + floor((p-blocksize) / blocksize)
     remainder <- p - n.standard * blocksize
     
          if(remainder==0)
          {
          beta.beg <- c(1,seq((blocksize+1), p, blocksize))
          beta.fin <- seq(blocksize, p, blocksize)
          n.beta.block <- length(beta.beg)
          }else
          {
          beta.beg <- c(1, seq((blocksize+1), p, blocksize))
          beta.fin <- c(seq((blocksize), p, blocksize), p)
          n.beta.block <- length(beta.beg)
          }
     }         

     
## Matrices to store samples
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, p))
samples.phi <- array(NA, c(n.keep, n))
samples.tau2 <- array(NA, c(n.keep, 1))
samples.nu2 <- array(NA, c(n.keep, 1))
samples.rho <- array(NA, c(n.keep, 1))
samples.Z <- array(NA, c(n.keep, n))
samples.lambda <- array(NA, c(n.keep, G))
samples.delta <- array(NA, c(n.keep, 1))
samples.deviance <- array(NA, c(n.keep, 1))
samples.fitted <- array(NA, c(n.keep, n))
     
         
## Metropolis quantities
accept.all <- rep(0,6)
accept <- accept.all
proposal.sd.rho <- 0.02
proposal.sd.delta <- 0.1
proposal.sd.lambda <- 0.01
tau2.posterior.shape <- prior.tau2[1] + 0.5 * n
nu2.posterior.shape <- prior.nu2[1] + 0.5*n

      
     
#### CAR quantities
    if(!is.matrix(W)) stop("W is not a matrix.", call.=FALSE)
    if(nrow(W)!= n) stop("W has the wrong number of rows.", call.=FALSE)
    if(ncol(W)!= n) stop("W has the wrong number of columns.", call.=FALSE)
    if(sum(is.na(W))>0) stop("W has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(W)) stop("W has non-numeric values.", call.=FALSE)
    if(min(W)<0) stop("W has negative elements.", call.=FALSE)
    if(sum(W!=t(W))>0) stop("W is not symmetric.", call.=FALSE)


## Create the triplet object
W.triplet <- c(NA, NA, NA)
     for(i in 1:n)
     {
          for(j in 1:n)
          {
               if(W[i,j]==1)
               {
               W.triplet <- rbind(W.triplet, c(i,j, W[i,j]))     
               }else{}
          }
     }
W.triplet <- W.triplet[-1, ]     
n.triplet <- nrow(W.triplet) 
W.triplet.sum <- tapply(W.triplet[ ,3], W.triplet[ ,1], sum)
n.neighbours <- tapply(W.triplet[ ,3], W.triplet[ ,1], length)

     
## Create the start and finish points for W updating
W.begfin <- array(NA, c(n, 2))     
temp <- 1
     for(i in 1:n)
     {
     W.begfin[i, ] <- c(temp, (temp + n.neighbours[i]-1))
     temp <- temp + n.neighbours[i]
     }
     
             
## Create the determinant     
Wstar <- diag(apply(W,1,sum)) - W
Wstar.eigen <- eigen(Wstar)
Wstar.val <- Wstar.eigen$values
det.Q <-  0.5 * sum(log((rho * Wstar.val + (1-rho))))    


#### Beta update quantities
data.precision.beta <- t(X.standardised) %*% X.standardised
data.var.beta <- solve(data.precision.beta)
data.temp.beta <- data.var.beta %*% t(X.standardised)
     if(length(prior.var.beta)==1)
	{
	prior.precision.beta <- 1 / prior.var.beta
	}else
	{
	prior.precision.beta <- solve(diag(prior.var.beta))
	}
     
     
###########################
#### Run the Bayesian model
###########################
## Start timer
     if(verbose)
     {
     cat("Collecting", n.sample, "samples\n", sep = " ")
     progressBar <- txtProgressBar(style = 3)
     percentage.points<-round((1:100/100)*n.sample)
     }else
     {
     percentage.points<-round((1:100/100)*n.sample)     
     }
     

     
    for(j in 1:n.sample)
    {
    ####################
    ## Sample from beta
    ####################
    data.mean.beta <- data.temp.beta %*% (Y - phi - offset - lambda[Z])
    U <- chol((prior.precision.beta + data.precision.beta / nu2))
    Uinv <- backsolve(U, diag(rep(1,p)))
    fc.mean.beta <- Uinv %*% (t(Uinv) %*% (prior.precision.beta %*% prior.mean.beta + (data.precision.beta / nu2) %*% data.mean.beta))
    beta <- fc.mean.beta + Uinv %*% rnorm(p)
    beta[which(X.indicator==2)] <- 0
    beta.regression <- as.numeric(X.standardised %*% beta)         
         
   

    ##################
    ## Sample from nu2
    ##################
    fitted.current <-  beta.regression + phi + offset + lambda[Z]
    nu2.posterior.scale <- prior.nu2[2] + 0.5 * sum((Y - fitted.current)^2)
    nu2 <- 1 / rgamma(1, nu2.posterior.shape, scale=(1/nu2.posterior.scale))    

         
         
    ####################
    ## Sample from phi
    ####################
    offset.phi <- (Y - beta.regression - offset - lambda[Z]) / nu2    
    phi <- gaussiancarupdate(Wtriplet=W.triplet, Wbegfin=W.begfin, W.triplet.sum, nsites=n, phi=phi, tau2=tau2, rho=rho, nu2=nu2, offset=offset.phi)
         for(i in 1:G)
          {
          phi[which(Z==i)] <- phi[which(Z==i)] - mean(phi[which(Z==i)])
          }


     ##################
     ## Sample from tau2
     ##################
     temp2 <- quadform(W.triplet, W.triplet.sum, n.triplet, n, phi, phi, rho)
     tau2.posterior.scale <- temp2 + prior.tau2[2] 
     tau2 <- 1 / rgamma(1, tau2.posterior.shape, scale=(1/tau2.posterior.scale))
              
               

     ##################
     ## Sample from rho
     ##################
     proposal.rho <- rtrunc(n=1, spec="norm", a=0, b=1, mean=rho, sd=proposal.sd.rho)  
     temp3 <- quadform(W.triplet, W.triplet.sum, n.triplet, n, phi, phi, proposal.rho)
     det.Q.proposal <- 0.5 * sum(log((proposal.rho * Wstar.val + (1-proposal.rho))))              
     logprob.current <- det.Q - temp2 / tau2
     logprob.proposal <- det.Q.proposal - temp3 / tau2
     prob <- exp(logprob.proposal - logprob.current)
    
     #### Accept or reject the proposal
         if(prob > runif(1))
         {
         rho <- proposal.rho
         det.Q <- det.Q.proposal
         accept[1] <- accept[1] + 1           
         }else
         {
         }              
     accept[2] <- accept[2] + 1           
         
         
         
    ####################
    ## Sample from lambda
    ####################
    proposal <- c(-1000, lambda, 1000)
          for(i in 1:G)
          {
           proposal[(i+1)] <- rtrunc(n=1, spec="norm", a=proposal[i], b=proposal[(i+2)], mean=proposal[(i+1)], sd=proposal.sd.lambda)    
          }
     proposal <- proposal[2:(G+1)]
     prob <- exp(sum((Y - lambda[Z] - phi - beta.regression - offset)^2 - (Y - proposal[Z] - phi - beta.regression - offset)^2)/  (2*nu2))
          if(prob > runif(1))
          {
          lambda <- proposal
          accept[3] <- accept[3] + 1  
          }else
          {
          }
     accept[4] <- accept[4] + 1       
     
    
         
    ################
    ## Sample from Z
    ################
    Z.offset <- phi + offset + beta.regression
    Z.proposal <- sample(1:G, size=n, replace=TRUE)
    prior <- delta * ((Z - Gbar)^2 - (Z.proposal-Gbar)^2)    
    like <- ((Y - lambda[Z] - Z.offset)^2 - (Y - lambda[Z.proposal] - Z.offset)^2) / (2*nu2)
    prob <- exp(like + prior)   
    test <- prob> runif(n)         
    Z[test] <- Z.proposal[test]         

         
       
    ##################
    ## Sample from delta
    ##################
    proposal.delta <-  rtrunc(n=1, spec="norm", a=0, b=prior.delta, mean=delta, sd=proposal.sd.delta)    
    prob1 <- sum((Z-Gbar)^2) * (delta - proposal.delta)        
    prob2 <- n * log(sum(exp(-delta *(1:G - Gbar)^2))) - n * log(sum(exp(-proposal.delta *(1:G - Gbar)^2)))
    prob <- exp(prob1 + prob2)    
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
    fitted <- lambda[Z] + phi + beta.regression + offset
    deviance.all <- dnorm(Y, mean = fitted, sd = rep(sqrt(nu2),n), log=TRUE)
    deviance <- -2 * sum(deviance.all)  
         
         
    ###################
    ## Save the results
    ###################
        if(j > burnin & (j-burnin)%%thin==0)
        {
        ele <- (j - burnin) / thin
        samples.beta[ele, ] <- beta
        samples.phi[ele, ] <- phi
        samples.lambda[ele, ] <- lambda
        samples.tau2[ele, ] <- tau2
        samples.nu2[ele, ] <- nu2
        samples.rho[ele, ] <- rho
        samples.Z[ele, ] <- Z
        samples.delta[ele, ] <- delta
        samples.deviance[ele, ] <- deviance
        samples.fitted[ele, ] <- fitted
        }else
        {
        }

         

    ########################################
    ## Self tune the acceptance probabilties
    ########################################
    k <- j/100
        if(ceiling(k)==floor(k))
        {
        #### Determine the acceptance probabilities
        accept.rho <- 100 * accept[1] / accept[2]
        accept.lambda <- 100 * accept[3] / accept[4]
        accept.delta <- 100 * accept[5] / accept[6]
        accept.all <- accept.all + accept
        accept <- rep(0,6)
            
               
       #### rho tuning parameter
               if(accept.rho > 70)
               {
               proposal.sd.rho <- min(2 * proposal.sd.rho, 0.5)
               }else if(accept.rho < 50)              
               {
               proposal.sd.rho <- 0.5 * proposal.sd.rho
               }else
               {
               }
             
        #### lambda tuning parameter
            if(accept.lambda > 40)
            {
            proposal.sd.lambda <- 2 * proposal.sd.lambda
            }else if(accept.lambda < 30)              
            {
            proposal.sd.lambda <- 0.5 * proposal.sd.lambda
            }else
            {
            }              

        #### delta tuning parameter
            if(accept.delta > 40)
            {
            proposal.sd.delta <- min(2 * proposal.sd.delta, prior.delta/6)
            }else if(accept.delta < 30)              
            {
            proposal.sd.delta <- 0.5 * proposal.sd.delta
            }else
            {
            }              
        }else
        {   
        }

        
    ################################       
    ## print progress to the console
    ################################
          if(j %in% percentage.points & verbose)
          {
          setTxtProgressBar(progressBar, j/n.sample)
          }
     }

# end timer
     if(verbose)
     {
     cat("\nSummarising results")
     close(progressBar)
     }else
     {}
     
###################################
#### Summarise and save the results 
###################################
## Compute the acceptance rates
accept.beta <- 100
accept.phi <- 100
accept.rho <- 100 * accept.all[1] / accept.all[2]
accept.lambda <- 100 * accept.all[3] / accept.all[4]
accept.delta <- 100 * accept.all[5] / accept.all[6]
accept.tau2 <- 100
accept.nu2 <- 100
accept.final <- c(accept.beta, accept.phi, accept.rho, accept.lambda, accept.delta, accept.tau2, accept.tau2)
names(accept.final) <- c("beta", "phi", "rho", "lambda", "delta", "tau2", "nu2")


## Deviance information criterion (DIC)
median.phi <- apply(samples.phi, 2, median)
median.beta <- apply(samples.beta,2,median)
median.Z <- round(apply(samples.Z,2,median),0)
median.lambda <- apply(samples.lambda,2,median)
fitted.median <- X.standardised %*% median.beta + median.lambda[median.Z] + median.phi + offset
nu2.median <- median(samples.nu2)
deviance.fitted <- -2 * sum(dnorm(Y, mean = fitted.median, sd = rep(sqrt(nu2.median),n), log = TRUE))
p.d <- median(samples.deviance) - deviance.fitted
DIC <- 2 * median(samples.deviance) - deviance.fitted       


#### Compute the Conditional Predictive Ordinate
CPO <- rep(NA, n)
     for(j in 1:n)
     {
     CPO[j] <- 1/median((1 / dnorm(Y[j], mean=samples.fitted[ ,j], sd=samples.nu2)))    
     }
LMPL <- sum(log(CPO))   

     
#### transform the parameters back to the origianl covariate scale.
samples.beta.orig <- samples.beta
number.cts <- sum(X.indicator==1)     
     if(number.cts>0)
     {
          for(r in 1:p)
          {
               if(X.indicator[r]==1)
               {
               samples.beta.orig[ ,r] <- samples.beta[ ,r] / X.sd[r]
               }else
               {
               }
          }
     }else
     {
     }

     
#### Shrink the lambda object to only include those levels used
Z.used <- as.numeric(names(table(samples.Z)))
Greal <- length(Z.used)     
samples.lambda <- samples.lambda[ ,Z.used]    
colnames(samples.lambda) <- Z.used
     

#### Create a summary object
samples.beta <- mcmc(samples.beta)
summary.beta <- t(apply(samples.beta, 2, quantile, c(0.5, 0.025, 0.975))) 
summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta, p))
rownames(summary.beta) <- colnames(X)
colnames(summary.beta) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept")
     
samples.lambda <- mcmc(samples.lambda)
summary.lambda <- t(apply(samples.lambda, 2, quantile, c(0.5, 0.025, 0.975))) 
summary.lambda <- cbind(summary.lambda, rep(n.keep, length(Z.used)), rep(accept.lambda, length(Z.used)))
colnames(summary.lambda) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept")

summary.hyper <- array(NA, c(4 ,5))
summary.hyper[1, 1:3] <- quantile(samples.nu2, c(0.5, 0.025, 0.975))
summary.hyper[1, 4:5] <- c(n.keep, accept.nu2)
summary.hyper[2, 1:3] <- quantile(samples.tau2, c(0.5, 0.025, 0.975))
summary.hyper[2, 4:5] <- c(n.keep, accept.tau2)
summary.hyper[3, 1:3] <- quantile(samples.rho, c(0.5, 0.025, 0.975))
summary.hyper[3, 4:5] <- c(n.keep, accept.rho)
summary.hyper[4, 1:3] <- quantile(samples.delta, c(0.5, 0.025, 0.975))
summary.hyper[4, 4:5] <- c(n.keep, accept.delta)
         
summary.results <- rbind(summary.beta, summary.lambda, summary.hyper)
rownames(summary.results)[(p+length(Z.used)+1):(p+length(Z.used)+4)] <- c("nu2", "tau2", "rho", "delta")
summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
summary.results[ , 4:5] <- round(summary.results[ , 4:5], 1)
summary.results <- summary.results[-which(X.indicator==2), ]


#### Create the Fitted values and residuals
fitted.values <- apply(samples.fitted, 2, median)
residuals <- as.numeric(Y) - fitted.values

     
## Compile and return the results
modelfit <- c(DIC, p.d, LMPL)
names(modelfit) <- c("DIC", "p.d", "LMPL")     
     
     
model.string <- c("Likelihood model - Gaussian (identity link function)", "\nRandom effects  model - Leroux CAR with clusters\n")
     if(length(which(X.indicator==2))==p)
     {
     samples <- list(phi=mcmc(samples.phi), lambda=mcmc(samples.lambda), Z=mcmc(samples.Z), tau2=mcmc(samples.tau2), nu2=mcmc(samples.nu2), rho=mcmc(samples.rho), delta=mcmc(samples.delta), fitted=mcmc(samples.fitted))
     }else
     {
     samples <- list(beta=mcmc(samples.beta.orig[ ,-which(X.indicator==2)]), phi=mcmc(samples.phi), lambda=mcmc(samples.lambda), Z=mcmc(samples.Z), tau2=mcmc(samples.tau2), nu2=mcmc(samples.nu2), rho=mcmc(samples.rho), delta=mcmc(samples.delta), fitted=mcmc(samples.fitted))          
     }
results <- list(summary.results=summary.results, samples=samples, fitted.values=fitted.values, residuals=residuals, modelfit=modelfit, accept=accept.final, localised.structure=median.Z,  formula=formula, model=model.string, X=X)
class(results) <- "carbayes"

     if(verbose)
     {
     b<-proc.time()
     cat(" finished in ", round(b[3]-a[3], 1), "seconds")
     }else
     {}
return(results)
}
