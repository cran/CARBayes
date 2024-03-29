S.CARbym <- function(formula, formula.omega=NULL, family, data=NULL,  trials=NULL, W, burnin, n.sample, thin=1, n.chains=1, n.cores=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.tau2=NULL, prior.sigma2=NULL, prior.mean.delta=NULL, prior.var.delta=NULL, MALA=TRUE, verbose=TRUE)
{
#### This is a wrapper function that calls one of
## binomial.bymCAR
## poisson.bymCAR
## zip.bymCAR
      if(is.null(family)) stop("the family argument is missing", call.=FALSE)
     
#### Run the appropriate model according to the family arugment
     if(family=="binomial")
     {
          if(is.null(trials)) stop("a binomial model was specified but the trials arugment was not specified", call.=FALSE)
         if(!is.null(formula.omega)) stop("you do not need a formula.omega argument as the zip model was not specified", call.=FALSE)
    model <- binomial.bymCAR(formula=formula, data=data,  trials=trials, W=W, burnin=burnin, n.sample=n.sample, thin=thin, n.chains=n.chains, n.cores=n.cores, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.tau2=prior.tau2, prior.sigma2=prior.sigma2, MALA=MALA, verbose=verbose)
     }else if(family=="poisson")
     {
          if(!is.null(trials)) stop("you do not need a trials arugment as a binomial model was not specified", call.=FALSE)
         if(!is.null(formula.omega)) stop("you do not need a formula.omega argument as the zip model was not specified", call.=FALSE)
     model <- poisson.bymCAR(formula=formula, data=data,  W=W, burnin=burnin, n.sample=n.sample, thin=thin, n.chains=n.chains, n.cores=n.cores, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.tau2=prior.tau2, prior.sigma2=prior.sigma2, MALA=MALA, verbose=verbose)          
     }else if(family=="zip")
     {
         if(!is.null(trials)) stop("you do not need a trials arugment as a binomial model was not specified", call.=FALSE)
         if(is.null(formula.omega)) stop("a zip model was specified but the formula.omega argument was not specified", call.=FALSE)
         model <- zip.bymCAR(formula=formula, formula.omega=formula.omega, data=data,  W=W, burnin=burnin, n.sample=n.sample, thin=thin, n.chains=n.chains, n.cores=n.cores, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.tau2=prior.tau2, prior.sigma2=prior.sigma2, prior.mean.delta=prior.mean.delta, prior.var.delta=prior.var.delta, MALA=MALA, verbose=verbose)          
     }else
     {
     stop("the family arugment is not one of `binomial', `poisson', or `zip'.", call.=FALSE)     
     }
return(model)
}