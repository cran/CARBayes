S.glm <- function(formula, family, data=NULL,  trials=NULL, burnin, n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.nu2=NULL, verbose=TRUE)
{
    #### This is a wrapper function that calls one of
    ## binomial.glm
    ## gaussian.glm
    ## poisson.glm
    if(is.null(family)) stop("the family argument is missing", call.=FALSE)
    
    #### Run the appropriate model according to the family arugment
    if(family=="binomial")
    {
        if(is.null(trials)) stop("a binomial model was specified but the trials arugment was not specified", call.=FALSE)
        model <- binomial.glm(formula=formula, data=data,  trials=trials, burnin=burnin, n.sample=n.sample, thin=thin, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, verbose=verbose)
    }else if(family=="gaussian")
    {
        if(!is.null(trials)) stop("you do not need a trials arugment as a binomial model was not specified", call.=FALSE)
        model <- gaussian.glm(formula=formula, data=data,  burnin=burnin, n.sample=n.sample, thin=thin, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.nu2=prior.nu2, verbose=verbose)          
    }else if(family=="poisson")
    {
        if(!is.null(trials)) stop("you do not need a trials arugment as a binomial model was not specified", call.=FALSE)
        model <- poisson.glm(formula=formula, data=data, burnin=burnin, n.sample=n.sample, thin=thin, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, verbose=verbose)          
    }else
    {
        stop("the family arugment is not one of `binomial', `gaussian' or `poisson'.", call.=FALSE)     
    }
    return(model)     
}