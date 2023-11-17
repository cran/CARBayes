S.RAB <- function(formula, family, data=NULL, trials=NULL, W, V, nlambda=100, verbose=TRUE)
{
    #### This is a wrapper function that calls one of
    ## binomial.RAB
    ## gaussian.RAB
    ## poisson.RAB
    if(is.null(family)) stop("the family argument is missing", call.=FALSE)
    
    #### Run the appropriate model according to the family arugment
    if(family=="binomial")
    {
        if(is.null(trials)) stop("a binomial model was specified but the trials arugment was not specified", call.=FALSE)
        model <- binomial.RAB(formula=formula, data=data,  trials=trials, W=W, V=V, nlambda=nlambda, verbose=verbose)
    }else if(family=="gaussian")
    {
        if(!is.null(trials)) stop("you do not need a trials arugment as a binomial model was not specified", call.=FALSE)
        model <- gaussian.RAB(formula=formula, data=data, W=W, V=V, nlambda=nlambda, verbose=verbose)          
    }else if(family=="poisson")
    {
        if(!is.null(trials)) stop("you do not need a trials arugment as a binomial model was not specified", call.=FALSE)
        model <- poisson.RAB(formula=formula, data=data, W=W, V=V, nlambda=nlambda, verbose=verbose)          
    }else
    {
        stop("the family arugment is not one of `binomial', `gaussian', or `poisson'.", call.=FALSE)     
    }
    return(model)     
}