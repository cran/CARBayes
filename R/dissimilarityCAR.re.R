dissimilarityCAR.re <- function(formula, family=NULL, data=NULL,  trials=NULL, rho=0.99, W, Z, burnin=0, n.sample=1000, thin=1, blocksize.beta=5, prior.mean.beta=NULL, prior.var.beta=NULL, prior.nu2=NULL, prior.tau2=NULL, verbose=TRUE)
{
#### This is a wrapper function that calls one of
## binomial.properCAR
## gaussian.properCAR
## poisson.properCAR
      if(is.null(family)) stop("the family argument is missing", call.=FALSE)
     
#### Run the appropriate model according to the family arugment
     if(family=="binomial")
     {
          if(is.null(trials)) stop("a binomial model was specified but the trials arugment was not specified", call.=FALSE)
     model <- binomial.dissimilarityCAR(formula=formula, data=data,  trials=trials, rho=rho, W=W, Z=Z, burnin=burnin, n.sample=n.sample, thin=thin, blocksize.beta=blocksize.beta, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.tau2=prior.tau2, verbose=verbose)
     }else if(family=="gaussian")
     {
          if(!is.null(trials)) stop("you do not need a trials arugment as a binomial model was not specified", call.=FALSE)
     model <- gaussian.dissimilarityCAR(formula=formula, data=data,  rho=rho, W=W, Z=Z, burnin=burnin, n.sample=n.sample, thin=thin, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.nu2=prior.nu2, prior.tau2=prior.tau2, verbose=verbose)          
     }else if(family=="poisson")
     {
          if(!is.null(trials)) stop("you do not need a trials arugment as a binomial model was not specified", call.=FALSE)
     model <- poisson.dissimilarityCAR(formula=formula, data=data,  rho=rho, W=W, Z=Z, burnin=burnin, n.sample=n.sample, thin=thin, blocksize.beta=blocksize.beta, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.tau2=prior.tau2, verbose=verbose)          
     }else
     {
     stop("the family arugment is not one of `binomial', `gaussian' or `poisson'.", call.=FALSE)     
     }
return(model)     
}