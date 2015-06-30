S.CARlocalised <- function(formula, exposure=NULL, family,  data=NULL, G, trials=NULL, W, burnin, n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.mean.alpha=NULL, prior.var.alpha=NULL, prior.tau2=NULL, prior.delta=NULL, verbose=TRUE)
{
#### This is a wrapper function that calls one of
## binomial.localisedCAR
## poisson.localisedCAR
## poisson.localisedCARagg
      if(is.null(family)) stop("the family argument is missing", call.=FALSE)
     
#### Run the appropriate model according to the family arugment
     if(family=="binomial")
     {
          if(is.null(trials)) stop("a binomial model was specified but the trials arugment was not specified", call.=FALSE)
     model <- binomial.localisedCAR(formula=formula, data=data, G=G, trials=trials, W=W, burnin=burnin, n.sample=n.sample, thin=thin, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.tau2=prior.tau2, prior.delta=prior.delta, verbose=verbose)
     }else if(family=="poisson" & is.null(exposure))
     {
          if(!is.null(trials)) stop("you do not need a trials arugment as a binomial model was not specified", call.=FALSE)
      model <- poisson.localisedCAR(formula=formula, data=data, G=G, W=W, burnin=burnin, n.sample=n.sample, thin=thin, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.tau2=prior.tau2, prior.delta=prior.delta, verbose=verbose)          
     }else if(family=="poisson" & !is.null(exposure))
     {
          if(!is.null(trials)) stop("you do not need a trials arugment as a binomial model was not specified", call.=FALSE)
      model <- poisson.localisedCARagg(formula=formula, exposure=exposure, data=data, G=G, W=W, burnin=burnin, n.sample=n.sample, thin=thin, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.mean.alpha=prior.mean.alpha, prior.var.alpha=prior.var.alpha, prior.tau2=prior.tau2, prior.delta=prior.delta, verbose=verbose)                
     }else if(family!="poisson" & !is.null(exposure))
     {
     stop("the within area variation in exposure model is only for a poisson likelihood.", call.=FALSE)     
     }else
     {
     stop("the family arugment is not one of `binomial' or `poisson'.", call.=FALSE)     
     }
return(model)     
}