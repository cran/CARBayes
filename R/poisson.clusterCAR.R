poisson.clusterCAR <- function(Y, expected, q, W, burnin=0, n.sample=1000, thin=1, prior.tau2=NULL, prior.rho=NULL, verbose=TRUE)
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
#### Data
logSMR <- log(Y / expected)
logSMR[Y==0] <- -4.5
n <- length(Y)   

    
    if(sum(is.na(expected))>0) stop("the expected numbers of cases has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(expected)) stop("the numbers of expected cases has non-numeric values.", call.=FALSE)
    if(min(expected)<=0) stop("the numbers of expected cases has zero or negative values.", call.=FALSE)

    if(sum(is.na(Y))>0) stop("the response has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(Y)) stop("the response variable has non-numeric values.", call.=FALSE)
int.check <- n-sum(ceiling(Y)==floor(Y))
    if(int.check > 0) stop("the respons variable has non-integer values.", call.=FALSE)
    if(min(Y)<0) stop("the response variable has negative values.", call.=FALSE)
 

#### Initial parameter values
## Cluster allocation
cluster.factor <- kmeans(x=logSMR, centers=q, nstart=1000)$cluster
X.temp <- array(NA, c(n,q))
          for(i in 1:q)
          {
          X.temp[ ,i] <- as.numeric(cluster.factor==i)    
          }
classes.means <- tapply(logSMR, cluster.factor, mean)     
X <- as.matrix(X.temp[ ,order(classes.means)], nrow=n)
x <- cluster.factor
     for(i in 1:q)
     {
     x[cluster.factor==order(classes.means)[i]] <-  i    
     }

## Other parameters
beta <- lm(logSMR~X-1)$coefficients
phi <- rnorm(n=n, mean=rep(0,n), sd=rep(0.1, n))
     for(i in 1:q)
     {
     phi[which(x==i)] <- phi[which(x==i)] - mean(phi[which(x==i)])
     }
tau2 <- runif(1)

    if(is.null(prior.rho)) prior.rho <- seq(0, 0.99, 0.01)
    if(min(prior.rho) < 0 | max(prior.rho) >=1) stop("prior.rho is outside the interval [0,1).", call.=FALSE)
    if(sum(duplicated(prior.rho))>0) stop("prior.rho has duplicate values.", call.=FALSE)     
rho <- sample(prior.rho, size=1)
which.rho <- which(rho==prior.rho)    

     
#### Priors
    if(is.null(prior.tau2)) prior.tau2 <- c(0.001, 0.001)
    if(length(prior.tau2)!=2) stop("the prior value for tau2 is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.tau2)) stop("the prior value for tau2 is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.tau2))!=0) stop("the prior value for tau2 has missing values.", call.=FALSE)    
    
     
#### MCMC quantities
## Checks
    if(!is.numeric(burnin)) stop("burn-in is not a number", call.=FALSE)
    if(!is.numeric(n.sample)) stop("n.sample is not a number", call.=FALSE)    
    if(n.sample <= 0) stop("n.sample is less than or equal to zero.", call.=FALSE)
    if(burnin < 0) stop("burn-in is less than zero.", call.=FALSE)
    if(n.sample <= burnin)  stop("Burn-in is greater than n.sample.", call.=FALSE)
    if(!is.numeric(thin)) stop("thin is not a number", call.=FALSE)
    if(thin <= 0) stop("thin is less than or equal to zero.", call.=FALSE)


## Matrices to store samples
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, q))
samples.phi <- array(NA, c(n.keep, n))
samples.x <- array(NA, c(n.keep, n))
samples.tau2 <- array(NA, c(n.keep, 1))
samples.rho <- array(NA, c(n.keep, 1))
samples.deviance <- array(NA, c(n.keep, n))
samples.fitted <- array(NA, c(n.keep, n))
samples.WBSS <- array(NA, c(n.keep, 2))  
     
## Metropolis quantities
accept.all <- rep(0,6)
accept <- accept.all
proposal.sd.phi <- 0.1
proposal.sd.beta <- 0.01
tau2.posterior.shape <- prior.tau2[1] + 0.5 * n
rho.step <- 5
      
     
#### CAR quantities
    if(!is.matrix(W)) stop("W is not a matrix.", call.=FALSE)
    if(nrow(W)!= n) stop("W has the wrong number of rows.", call.=FALSE)
    if(ncol(W)!= n) stop("W has the wrong number of columns.", call.=FALSE)
    if(sum(is.na(W))>0) stop("W has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(W)) stop("W has non-numeric values.", call.=FALSE)
    if(min(W)<0) stop("W has negative elements.", call.=FALSE)
    if(sum(W!=t(W))>0) stop("W is not symmetric.", call.=FALSE)


## Create the duplet form
n.neighbours <- as.numeric(apply(W, 1, sum))
W.duplet <- c(NA, NA)
     for(i in 1:n)
     {
          for(j in 1:n)
          {
               if(W[i,j]==1)
               {
               W.duplet <- rbind(W.duplet, c(i,j))     
               }else{}
          }
     }
W.duplet <- W.duplet[-1, ]     
n.duplet <- nrow(W.duplet) 


## Create the list object
Wlist <- as.list(rep(NA,n))     
     for(i in 1:n)
     {
     Wlist[[i]] <- which(W[i, ]==1)     
     }
     
## Create the set of determinants     
Wstar <- diag(n.neighbours) - W
Wstar.eigen <- eigen(Wstar)
Wstar.val <- Wstar.eigen$values
n.rho <- length(prior.rho)
det.Q <- rep(NA, n.rho)
     for(i in 1:n.rho)
     {
     det.Q[i] <-  0.5 * sum(log((prior.rho[i] * Wstar.val + (1-prior.rho[i]))))    
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
    proposal <- c(-1000, beta, 1000)
          for(i in 1:q)
          {
           proposal[(i+1)] <- rtrunc(n=1, spec="norm", a=proposal[i], b=proposal[(i+2)], mean=proposal[(i+1)], sd=proposal.sd.beta)    
          }
     proposal <- proposal[2:(q+1)]
     lp.current <- beta[x] + phi
     lp.proposal <- proposal[x] + phi
     prob1 <- sum(expected * (exp(lp.current) - exp(lp.proposal)))
     prob2 <- sum(Y * (lp.proposal - lp.current))
     prob <- exp(prob1 + prob2)
          if(prob > runif(1))
          {
          beta <- proposal
          accept[1] <- accept[1] + 1  
          }else
          {
          }
     accept[2] <- accept[2] + 1       

     
    
    ################
    ## Sample from x
    ################
    offset <- expected * exp(phi)
    proposal <- sample(1:q, size=n, replace=TRUE)
    temp1 <- poissonxupdate(y=Y, offset=offset, x=x, nsites=n, beta=beta, proposal=proposal)
    x <- temp1
         
         
       
    ####################
    ## Sample from phi
    ####################
    beta.offset <- beta[x] + log(expected)
    temp1 <- poissoncarupdate(W_list=Wlist, nsites=n, phi=phi, nneighbours=n.neighbours, tau2=tau2, y=Y, phi_tune=proposal.sd.phi, rho_num=rho, rho_den=rho, offset=beta.offset)
    phi <- temp1[[1]]
          for(i in 1:q)
          {
          phi[which(x==i)] <- phi[which(x==i)] - mean(phi[which(x==i)])
          }
         
    accept[3] <- accept[3] + temp1[[2]]
    accept[4] <- accept[4] + n


         
    ##################
    ## Sample from tau2
    ##################
    temp2 <- quadform(W_duplet1=W.duplet[ ,1], W_duplet2=W.duplet[ ,2], n_duplet=n.duplet,  nsites=n, phi=phi, nneighbours=n.neighbours, diagonal=rho, offdiagonal=rho)      
    tau2.posterior.scale <- temp2 + prior.tau2[2] 
    tau2 <- 1 / rgamma(1, tau2.posterior.shape, scale=(1/tau2.posterior.scale))
         

    ##################
    ## Sample from rho
    ##################
    which.proposal.rho <- sample(setdiff(max(1, which.rho-rho.step):min(which.rho+rho.step, n.rho), which.rho), size=1)    
    proposal.rho <- prior.rho[which.proposal.rho]
    temp3 <- quadform(W_duplet1=W.duplet[ ,1], W_duplet2=W.duplet[ ,2], n_duplet=n.duplet,  nsites=n, phi=phi, nneighbours=n.neighbours, diagonal=proposal.rho, offdiagonal=proposal.rho)      
    logprob.current <- det.Q[which.rho] - temp2 / tau2
    logprob.proposal <- det.Q[which.proposal.rho] - temp3 / tau2
    prob <- exp(logprob.proposal - logprob.current)
    
    #### Accept or reject the proposal
         if(prob > runif(1))
         {
         rho <- proposal.rho
         which.rho <- which.proposal.rho
         accept[5] <- accept[5] + 1           
         }else
         {
         }              
    accept[6] <- accept[6] + 1           

   
         
    #########################
    ## Calculate the deviance
    #########################
    fitted <- expected * exp(beta[x] + phi)
    deviance <- dpois(x=as.numeric(Y), lambda=fitted)

         

     ####################################     
     ## Compute the between and within SS
     ####################################     
     BSS <- mean((beta[2:q] - beta[1:(q-1)])^2)     
     WSS <- mean((logSMR - beta[x])^2)          
     WBSS <- WSS / BSS
     WBSS2 <- (BSS/(q-1)) /  (WSS/(n-q))  
         
         
    ###################
    ## Save the results
    ###################
        if(j > burnin & (j-burnin)%%thin==0)
        {
        ele <- (j - burnin) / thin
        samples.beta[ele, ] <- beta
        samples.phi[ele, ] <- phi
        samples.tau2[ele, ] <- tau2
        samples.rho[ele, ] <- rho
        samples.x[ele, ] <- x
        samples.deviance[ele, ] <- deviance
        samples.fitted[ele, ] <- fitted
        samples.WBSS[ele, ] <- c(WBSS, WBSS2)
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
        accept.beta <- 100 * accept[1] / accept[2]
        accept.phi <- 100 * accept[3] / accept[4]
        accept.rho <- 100 * accept[5] / accept[6]
        accept.all <- accept.all + accept
        accept <- c(0,0,0,0,0,0)
            
        #### beta tuning parameter
            if(accept.beta > 70)
            {
            proposal.sd.beta <- 2 * proposal.sd.beta
            }else if(accept.beta < 50)              
            {
            proposal.sd.beta <- 0.5 * proposal.sd.beta
            }else
            {
            }
            
        #### phi tuning parameter
            if(accept.phi > 40)
            {
            proposal.sd.phi <- 2 * proposal.sd.phi
            }else if(accept.phi < 30)              
            {
            proposal.sd.phi <- 0.5 * proposal.sd.phi
            }else
            {
            }              
               
       #### rho tuning parameter
               if(accept.rho > 70)
               {
               rho.step <- min(rho.step+1, round(n.rho/4))
               }else if(accept.rho < 50)              
               {
               proposal.sd.rho <- max(rho.step-1, 1)
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
accept.beta <- 100 * accept.all[1] / accept.all[2]
accept.phi <- 100 * accept.all[3] / accept.all[4]
accept.rho <- 100 * accept.all[5] / accept.all[6]
accept.tau2 <- 100
accept.final <- c(accept.beta, accept.phi, accept.rho, accept.tau2)
names(accept.final) <- c("beta", "phi", "rho", "tau2")

## Deviance information criterion (DIC)
median.phi <- apply(samples.phi, 2, median)
median.beta <- apply(samples.beta,2,median)
median.x <- round(apply(samples.x,2,median),0)
fitted.median <- expected * exp(median.beta[median.x] + median.phi)
deviance.fitted <- -2 * sum(dpois(x=Y, lambda=fitted.median, log=TRUE))
deviance.sum <- apply(-2 * log(samples.deviance), 1, sum)
p.d <- mean(deviance.sum) - deviance.fitted
DIC <- 2 * mean(deviance.sum) - deviance.fitted
like.fitted <- apply(samples.deviance, 2, mean)
DIC3 <- 2 * mean(deviance.sum)   + 2 * sum(log(like.fitted))     


#### Compute the Conditional Predictive Ordinate
CPO.temp <- 1 / samples.deviance
CPO <- 1/apply(CPO.temp, 2, mean)
MPL <- sum(log(CPO))    

     
 

#### Create a summary object
samples.beta <- mcmc(samples.beta)
summary.beta <- t(apply(samples.beta, 2, quantile, c(0.5, 0.025, 0.975))) 
summary.beta <- cbind(summary.beta, rep(n.keep, q), rep(accept.beta, q))
rownames(summary.beta) <- 1:q
colnames(summary.beta) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept")

summary.hyper <- array(NA, c(2 ,5))
summary.hyper[1, 1:3] <- quantile(samples.tau2, c(0.5, 0.025, 0.975))
summary.hyper[1, 4:5] <- c(n.keep, accept.tau2)
summary.hyper[2, 1:3] <- quantile(samples.rho, c(0.5, 0.025, 0.975))
summary.hyper[2, 4:5] <- c(n.keep, accept.rho)

summary.results <- rbind(summary.beta, summary.hyper)
rownames(summary.results)[(nrow(summary.results)-1):nrow(summary.results)] <- c("tau2", "rho")
summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
summary.results[ , 4:5] <- round(summary.results[ , 4:5], 1)



#### Create the Fitted values and residuals
fitted.values <- array(NA, c(n, 5))
colnames(fitted.values) <- c("Mean", "Sd", "Median", "2.5%", "97.5%")
fitted.values[ ,1] <- apply(samples.fitted, 2, mean)
fitted.values[ ,2] <- apply(samples.fitted, 2, sd)
fitted.values[ ,3:5] <- t(apply(samples.fitted, 2, quantile, c(0.5, 0.025, 0.975)))
fitted.values <- round(fitted.values, 4)

residuals <- array(NA, c(n, 5))
colnames(residuals) <- c("Mean", "Sd", "Median", "2.5%", "97.5%")
residuals.temp <- array(NA, c(nrow(samples.beta), n))
     for(i in 1:nrow(samples.beta))
     {
     residuals.temp[i, ] <- as.numeric(Y) - samples.fitted[i, ]
     }
residuals[ ,1] <- apply(residuals.temp, 2, mean)
residuals[ ,2] <- apply(residuals.temp, 2, sd)
residuals[ ,3:5] <- t(apply(residuals.temp, 2, quantile, c(0.5, 0.025, 0.975)))
residuals <- round(residuals, 4)


## Compile and return the results
WBSS.final <- median(samples.WBSS[ ,1])
CH.final <- median(samples.WBSS[ ,2])
modelfit <- c(DIC, p.d, DIC3, MPL, WBSS.final, CH.final)
names(modelfit) <- c("DIC", "p.d", "DIC3", "MPL", "WBSS", "CH")     
     
     
model.string <- c("Likelihood model - Poisson (log link function)", "\nRandom effects model - Leroux CAR\n")
formula <- c("Cluster model")
samples <- list(beta=mcmc(samples.beta), x=mcmc(samples.x), phi=mcmc(samples.phi), tau2=mcmc(samples.tau2), rho=mcmc(samples.rho))
results <- list(formula=formula, samples=samples, fitted.values=fitted.values, residuals=residuals, W.summary=W, modelfit=modelfit, summary.results=summary.results, model=model.string, accept=accept.final)
class(results) <- "carbayes"

     if(verbose)
     {
     b<-proc.time()
     cat(" finished in ", round(b[3]-a[3], 1), "seconds")
     }else
     {}
return(results)
}
