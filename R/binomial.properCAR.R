binomial.properCAR <-
function(formula, data=NULL, beta=NULL, phi=NULL, tau2=NULL, rho=NULL, trials, W, burnin=0, n.sample=1000, thin=1, blocksize.beta=5, prior.mean.beta=NULL, prior.var.beta=NULL, prior.tau2=NULL, prior.rho=NULL)
{
cat("Setting up the model\n")
a<-proc.time()

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



#### Response variable and trials
## Create the response
Y <- model.response(frame)
failures <- trials - Y
     
## Check for errors
    if(sum(is.na(trials))>0) stop("the numbers of trials has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(trials)) stop("the numbers of trials has non-numeric values.", call.=FALSE)
int.check <- n-sum(ceiling(trials)==floor(trials))
    if(int.check > 0) stop("the numbers of trials has non-integer values.", call.=FALSE)
    if(min(trials)<=0) stop("the numbers of trials has zero or negative values.", call.=FALSE)

    if(sum(is.na(Y))>0) stop("the response has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(Y)) stop("the response variable has non-numeric values.", call.=FALSE)
int.check <- n-sum(ceiling(Y)==floor(Y))
    if(int.check > 0) stop("the respons variable has non-integer values.", call.=FALSE)
    if(min(Y)<0) stop("the response variable has negative values.", call.=FALSE)
    if(sum(Y>trials)>0) stop("the response variable has larger values that the numbers of trials.", call.=FALSE)



#### Offset variable
## Create the offset
offset <- try(model.offset(frame), silent=TRUE)

## Check for errors
    if(class(offset)=="try-error")   stop("the offset is not numeric.", call.=FALSE)
    if(is.null(offset))  offset <- rep(0,n)
    if(sum(is.na(offset))>0) stop("the offset has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(offset)) stop("the offset variable has non-numeric values.", call.=FALSE)



#### Initial parameter values
## Regression parameters beta
dat <- cbind(Y, failures)

    if(is.null(beta)) beta <- glm(dat~X.standardised-1, offset=offset, family=binomial)$coefficients 
    if(length(beta)!= p) stop("beta is the wrong length.", call.=FALSE)
    if(sum(is.na(beta))>0) stop("beta has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(beta)) stop("beta has non-numeric values.", call.=FALSE)

## Random effects phi
    if(is.null(phi)) phi <- rnorm(n=n, mean=rep(0,n), sd=rep(0.1, n))
    if(length(phi)!= n) stop("phi is the wrong length.", call.=FALSE)
    if(sum(is.na(phi))>0) stop("phi has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(phi)) stop("phi has non-numeric values.", call.=FALSE)

## Random effects variance tau2
    if(is.null(tau2)) tau2 <- runif(1)
    if(length(tau2)!= 1) stop("tau2 is the wrong length.", call.=FALSE)
    if(sum(is.na(tau2))>0) stop("tau2 has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(tau2)) stop("tau2 has non-numeric values.", call.=FALSE)
    if(tau2 <= 0) stop("tau2 is negative or zero.", call.=FALSE)

## Global correlation parameter rho
    if(is.null(prior.rho)) prior.rho <- seq(0, 0.99, 0.01)
    if(min(prior.rho) < 0 | max(prior.rho) >=1) stop("prior.rho is outside the interval [0,1).", call.=FALSE)
    if(sum(duplicated(prior.rho))>0) stop("prior.rho has duplicate values.", call.=FALSE)     
    if(is.null(rho)) rho <- sample(prior.rho, size=1)
    if(sum(rho==prior.rho)!= 1) stop("rho is not one of the possible values.", call.=FALSE)
    if(length(rho)!= 1) stop("rho is the wrong length.", call.=FALSE)
    if(sum(is.na(rho))>0) stop("rho has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(rho)) stop("rho has non-numeric values.", call.=FALSE)
    if(rho < 0 | rho >=1) stop("rho is outside the interval [0,1).", call.=FALSE)
which.rho <- which(rho==prior.rho)    
    
    
#### Priors
## Put in default priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(1000, p)
    if(is.null(prior.tau2)) prior.tau2 <- c(0.001, 0.001)

    
## Checks    
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

     
     
#### MCMC quantities
## Checks
    if(!is.numeric(burnin)) stop("burn-in is not a number", call.=FALSE)
    if(!is.numeric(n.sample)) stop("n.sample is not a number", call.=FALSE)    
    if(n.sample <= 0) stop("n.sample is less than or equal to zero.", call.=FALSE)
    if(burnin < 0) stop("burn-in is less than zero.", call.=FALSE)
    if(n.sample <= burnin)  stop("Burn-in is greater than n.sample.", call.=FALSE)
    if(!is.numeric(thin)) stop("thin is not a number", call.=FALSE)
    if(thin <= 0) stop("thin is less than or equal to zero.", call.=FALSE)

    if(!is.numeric(blocksize.beta)) stop("blocksize.beta is not a number", call.=FALSE)
    if(blocksize.beta <= 0) stop("blocksize.beta is less than or equal to zero", call.=FALSE)
    if(!(floor(blocksize.beta)==ceiling(blocksize.beta))) stop("blocksize.beta has non-integer values.", call.=FALSE)


## Compute the blocking structure for beta     
     if(blocksize.beta >= p)
     {
     n.beta.block <- 1
     beta.beg <- 1
     beta.fin <- p
     }else
     {
     n.standard <- 1 + floor((p-blocksize.beta) / blocksize.beta)
     remainder <- p - n.standard * blocksize.beta
     
          if(remainder==0)
          {
          beta.beg <- c(1,seq((blocksize.beta+1), p, blocksize.beta))
          beta.fin <- seq(blocksize.beta, p, blocksize.beta)
          n.beta.block <- length(beta.beg)
          }else
          {
          beta.beg <- c(1, seq((blocksize.beta+1), p, blocksize.beta))
          beta.fin <- c(seq((blocksize.beta), p, blocksize.beta), p)
          n.beta.block <- length(beta.beg)
          }
     }     



## Matrices to store samples
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, p))
samples.phi <- array(NA, c(n.keep, n))
samples.tau2 <- array(NA, c(n.keep, 1))
samples.rho <- array(NA, c(n.keep, 1))
samples.deviance <- array(NA, c(n.keep, 1))


## Metropolis quantities
accept.all <- rep(0,6)
accept <- accept.all
proposal.sd.beta <- 0.01
proposal.sd.phi <- 0.1
proposal.corr.beta <- solve(t(X.standardised) %*% X.standardised)
chol.proposal.corr.beta <- chol(proposal.corr.beta) 
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
n.neighbours <- as.numeric(apply(W, 1, sum))
n.rho <- length(prior.rho)
det.Q <- rep(NA, n.rho)
     
     for(i in 1:n.rho)
     {
     Q <- diag(n.neighbours) - prior.rho[i] * W
     determinant <- determinant(Q, logarithm=TRUE)
     det.Q[i] <- 0.5 * determinant$modulus * determinant$sign
     }




###########################
#### Run the Bayesian model
###########################
## Start timer
cat("Collecting", n.sample, "samples\n", sep = " ")
progressBar <- txtProgressBar(style = 3)
percentage.points<-round((1:100/100)*n.sample)

     for(j in 1:n.sample)
    {
    ####################
    ## Sample from beta
    ####################
    proposal <- beta + (sqrt(proposal.sd.beta)* t(chol.proposal.corr.beta)) %*% rnorm(p)
    proposal.beta <- beta
    offset.temp <- phi + offset

       for(r in 1:n.beta.block)
       {
       proposal.beta[beta.beg[r]:beta.fin[r]] <- proposal[beta.beg[r]:beta.fin[r]]
       prob <- binomialbetaupdate(X.standardised, n, p, beta, proposal.beta, offset.temp, Y, failures, prior.mean.beta, prior.var.beta)
            if(prob > runif(1))
            {
            beta[beta.beg[r]:beta.fin[r]] <- proposal.beta[beta.beg[r]:beta.fin[r]]
            accept[1] <- accept[1] + 1  
            }else
            {
            proposal.beta[beta.beg[r]:beta.fin[r]] <- beta[beta.beg[r]:beta.fin[r]]
            }
        }

    accept[2] <- accept[2] + n.beta.block    
 


    ####################
    ## Sample from phi
    ####################
    beta.offset <- X.standardised %*% beta + offset
    temp1 <- binomialcarupdate(W_list=Wlist, nsites=n, phi=phi, nneighbours=n.neighbours, tau2=tau2, y=Y, failures=failures, phi_tune=proposal.sd.phi, rho_num=rho, rho_den=1, offset=beta.offset)
    phi <- temp1[[1]]
    phi <- phi - mean(phi)
    accept[3] <- accept[3] + temp1[[2]]
    accept[4] <- accept[4] + n    
    

    ##################
    ## Sample from tau2
    ##################
    temp2 <- quadform(W_duplet1=W.duplet[ ,1], W_duplet2=W.duplet[ ,2], n_duplet=n.duplet,  nsites=n, phi=phi, nneighbours=n.neighbours, diagonal=1, offdiagonal=rho)      
    tau2.posterior.scale <- temp2 + prior.tau2[2] 
    tau2 <- 1 / rgamma(1, tau2.posterior.shape, scale=(1/tau2.posterior.scale))
    
    
    
    ##################
    ## Sample from rho
    ##################
    which.proposal.rho <- sample(setdiff(max(1, which.rho-rho.step):min(which.rho+rho.step, n.rho), which.rho), size=1)    
    proposal.rho <- prior.rho[which.proposal.rho]
    temp3 <- quadform(W_duplet1=W.duplet[ ,1], W_duplet2=W.duplet[ ,2], n_duplet=n.duplet,  nsites=n, phi=phi, nneighbours=n.neighbours, diagonal=1, offdiagonal=proposal.rho)      
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
    logit <- as.numeric(X.standardised %*% beta) + phi + offset    
    prob <- exp(logit)  / (1 + exp(logit))
    deviance <- -2 * sum(dbinom(x=Y, size=trials, prob=prob, log=TRUE))



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
        samples.deviance[ele, ] <- deviance
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
          if(j %in% percentage.points)
          {
          setTxtProgressBar(progressBar, j/n.sample)
          }
     }

# end timer
cat("\nSummarising results")
close(progressBar)
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
median.beta <- apply(samples.beta, 2, median)
median.phi <- apply(samples.phi, 2, median)
median.logit <- as.numeric(X.standardised %*% median.beta) + median.phi + offset    
median.prob <- exp(median.logit)  / (1 + exp(median.logit))
fitted.median <- trials * median.prob
deviance.fitted <- -2 * sum(dbinom(x=Y, size=trials, prob=median.prob, log=TRUE))
p.d <- mean(samples.deviance) - deviance.fitted
DIC <- 2 * mean(samples.deviance) - deviance.fitted


#### Compute the Conditional Predictive Ordinate
CPO.temp <- array(NA, c(nrow(samples.phi), n))
    for(i in 1:nrow(samples.phi))
    {
    temp.logit <- samples.phi[i, ]  + X.standardised %*% samples.beta[i, ] + offset
    temp <- exp(temp.logit)  / (1 + exp(temp.logit))
    CPO.temp[i, ] <- 1 / dbinom(x=Y, size=trials, prob=temp)
    }
CPO <- 1/apply(CPO.temp, 2, mean)
MPL <- sum(log(CPO))     


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
    }else if(X.indicator[r]==2 & p>1)
    {
      X.transformed <- which(X.indicator==1)
      samples.temp <- as.matrix(samples.beta[ ,X.transformed])
      for(s in 1:length(X.transformed))
      {
        samples.temp[ ,s] <- samples.temp[ ,s] * X.mean[X.transformed[s]]  / X.sd[X.transformed[s]]
      }
      intercept.adjustment <- apply(samples.temp, 1,sum) 
      samples.beta.orig[ ,r] <- samples.beta[ ,r] - intercept.adjustment
    }else
    {
    }
  }
}else
{
}



#### Create a summary object
samples.beta.orig <- mcmc(samples.beta.orig)
summary.beta <- t(apply(samples.beta.orig, 2, quantile, c(0.5, 0.025, 0.975))) 
summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta,p))
rownames(summary.beta) <- colnames(X)
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



#### Create the random effects summary
random.effects <- array(NA, c(n, 5))
colnames(random.effects) <- c("Mean", "Sd", "Median", "2.5%", "97.5%")
random.effects[ ,1] <- apply(samples.phi, 2, mean)
random.effects[ ,2] <- apply(samples.phi, 2, sd)
random.effects[ ,3:5] <- t(apply(samples.phi, 2, quantile, c(0.5, 0.025, 0.975)))
random.effects <- round(random.effects, 4)



#### Create the Fitted values
fitted.values <- array(NA, c(n, 5))
residuals <- array(NA, c(n, 5))
colnames(fitted.values) <- c("Mean", "Sd", "Median", "2.5%", "97.5%")
colnames(residuals) <- c("Mean", "Sd", "Median", "2.5%", "97.5%")
fitted.temp <- array(NA, c(nrow(samples.beta), n))
residuals.temp <- array(NA, c(nrow(samples.beta), n))
    for(i in 1:nrow(samples.beta))
    {
    temp.logit <- X.standardised %*% samples.beta[i, ] + samples.phi[i, ] + offset    
    temp <- trials * exp(temp.logit)  / (1 + exp(temp.logit))
    fitted.temp[i, ] <- temp
    residuals.temp[i, ] <- Y - temp    
    }
fitted.values[ ,1] <- apply(fitted.temp, 2, mean)
fitted.values[ ,2] <- apply(fitted.temp, 2, sd)
fitted.values[ ,3:5] <- t(apply(fitted.temp, 2, quantile, c(0.5, 0.025, 0.975)))
fitted.values <- round(fitted.values, 4)
residuals[ ,1] <- apply(residuals.temp, 2, mean)
residuals[ ,2] <- apply(residuals.temp, 2, sd)
residuals[ ,3:5] <- t(apply(residuals.temp, 2, quantile, c(0.5, 0.025, 0.975)))
residuals <- round(residuals, 4)
     
     





## Compile and return the results
model.string <- c("Likelihood model - Binomial (logit link function)", "\nRandom effects model - Proper CAR\n")
samples <- list(beta=samples.beta.orig, phi=mcmc(samples.phi), tau2=mcmc(samples.tau2), rho=mcmc(samples.rho))
results <- list(formula=formula, samples=samples,  fitted.values=fitted.values, random.effects=random.effects, residuals=residuals, W.summary=W, DIC=DIC, p.d=p.d, MPL=MPL, summary.results=summary.results, model=model.string, accept=accept.final)
class(results) <- "carbayes"
b<-proc.time()
cat(" finished in ", round(b[3]-a[3], 1), "seconds")
return(results)
}

