S.CARmultilevel <- function(formula, family, data=NULL,  trials=NULL, W, ind.area,  burnin, n.sample, thin=1, n.chains=1, n.cores=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.nu2=NULL, prior.tau2=NULL, rho=NULL, MALA=TRUE, verbose=TRUE)
{
    #### This is a wrapper function that calls one of
    ## binomial.multilevelCAR
    ## gaussian.multilevelCAR
    ## poisson.multilevelCAR
    if(is.null(family)) stop("the family argument is missing", call.=FALSE)

    #### Run the appropriate model according to the family arugment
    if(family=="binomial")
    {
        if(is.null(trials)) stop("a binomial model was specified but the trials arugment was not specified", call.=FALSE)
        model <- binomial.multilevelCAR(formula=formula, data=data,  trials=trials, W=W, ind.area=ind.area, burnin=burnin, n.sample=n.sample, thin=thin, n.chains=n.chains, n.cores=n.cores, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.tau2=prior.tau2, rho=rho, MALA=MALA, verbose=verbose)
        }else if(family=="gaussian")
    {
        if(!is.null(trials)) stop("you do not need a trials arugment as a binomial model was not specified", call.=FALSE)
        model <- gaussian.multilevelCAR(formula=formula, data=data,  W=W, ind.area=ind.area, burnin=burnin, n.sample=n.sample, thin=thin, n.chains=n.chains, n.cores=n.cores, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.nu2=prior.nu2, prior.tau2=prior.tau2, rho=rho, verbose=verbose)          
    }else if(family=="poisson")
    {
        if(!is.null(trials)) stop("you do not need a trials arugment as a binomial model was not specified", call.=FALSE)
        model <- poisson.multilevelCAR(formula=formula, data=data,  W=W, ind.area=ind.area, burnin=burnin, n.sample=n.sample, thin=thin, n.chains=n.chains, n.cores=n.cores, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.tau2=prior.tau2, rho=rho, MALA=MALA, verbose=verbose)          
    }else
    {
        stop("the family arugment is not one of `binomial', `gaussian' or `poisson'.", call.=FALSE)     
    }
    return(model)     
}