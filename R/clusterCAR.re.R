clusterCAR.re <- function(Y, q, family=NULL, expected=NULL, trials=NULL, W, burnin=0, n.sample=1000, thin=1, prior.nu2=NULL, prior.tau2=NULL, prior.rho=NULL, verbose=TRUE)
{
#### This is a wrapper function that calls one of
## binomial.clusterCAR
## gaussian.clusterCAR
## poisson.clusterCAR
      if(is.null(family)) stop("the family argument is missing", call.=FALSE)
     
#### Run the appropriate model according to the family arugment
     if(family=="binomial")
     {
          if(is.null(trials)) stop("a binomial model was specified but the trials arugment was not specified", call.=FALSE)
          if(!is.null(expected)) stop("you do not need the expected arugment as a poisson model was not specified", call.=FALSE)
     model <- binomial.clusterCAR(Y=Y, q=q,  trials=trials, W=W, burnin=burnin, n.sample=n.sample, thin=thin, prior.tau2=prior.tau2, prior.rho=prior.rho, verbose=verbose)
     }else if(family=="gaussian")
     {
          if(!is.null(trials)) stop("you do not need a trials arugment as a binomial model was not specified", call.=FALSE)
          if(!is.null(expected)) stop("you do not need the expected arugment as a poisson model was not specified", call.=FALSE)
     model <- gaussian.clusterCAR(Y=Y, q=q,  W=W, burnin=burnin, n.sample=n.sample, thin=thin, prior.nu2=prior.nu2, prior.tau2=prior.tau2, prior.rho=prior.rho, verbose=verbose)          
     }else if(family=="poisson")
     {
          if(is.null(expected)) stop("a poisson model was specified but the expected arugment was not specified", call.=FALSE)
          if(!is.null(trials)) stop("you do not need a trials arugment as a binomial model was not specified", call.=FALSE)
     model <- poisson.clusterCAR(Y=Y, q=q, expected=expected,  W=W, burnin=burnin, n.sample=n.sample, thin=thin, prior.tau2=prior.tau2, prior.rho=prior.rho, verbose=verbose)          
     }else
     {
     stop("the family arugment is not one of `binomial', `gaussian' or `poisson'.", call.=FALSE)     
     }
return(model)     
}