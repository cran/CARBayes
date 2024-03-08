poisson.RAB <- function(formula, data=NULL, W, V, nlambda, verbose=TRUE)
{
##############################################
#### Format the arguments and check for errors
##############################################
#### Verbose
a <- common.verbose(verbose)
    
    
#### Frame object
frame.results <- common.frame(formula, data, "gaussian")
K <- frame.results$n
p <- frame.results$p
X <- frame.results$X
offset <- frame.results$offset
Y <- frame.results$Y
which.miss <- frame.results$which.miss
which.present <- which(!is.na(Y))
n.miss <- frame.results$n.miss  
    if(p==0) stop("The model (via the formula object) must at least have an intercept term.", call.=FALSE)


#### Ancillary data
    if(!is.numeric(V)) stop("The ancillary data V is not a vector.", call.=FALSE)
    if(length(V) != K) stop("The ancillary data V is not the same length as the remaining data.", call.=FALSE)
    if(sum(is.na(V))>0) stop("The ancillary data V has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(V)) stop("The ancillary data V has non-numeric values.", call.=FALSE)


#### Neighbourhood matrix W
    if(!is.matrix(W)) stop("W is not a matrix.", call.=FALSE)
    if(ncol(W)!= nrow(W)) stop("W is not a square matrix.", call.=FALSE)    
    if(sum(is.na(W))>0) stop("W has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(W)) stop("W has non-numeric values.", call.=FALSE)
    if(min(W)<0) stop("W has negative elements.", call.=FALSE)
    if(sum(W!=t(W))>0) stop("W is not symmetric.", call.=FALSE)
    if(min(apply(W, 1, sum))==0) stop("W has some areas with no neighbours (one of the row sums equals zero).", call.=FALSE)    


#### Create the shortest path matrix
graph.W <- graph.adjacency(W, mode="undirected")
graph.dist <- shortest.paths(graph.W)  



#####################################################
#### Create the basis functions and the data elements
#####################################################
#### Create the three sets of basis functions
B.anisotropic.exp <- basiscomputeexponential(D=graph.dist, nrows=K, ncols=K, Z=V, startcol=1) 
B.anisotropic.inv <- basiscomputeinverse(D=graph.dist, nrows=K, ncols=K, Z=V, startcol=1) 
B.anisotropic.linear <- basiscomputelinear(D=graph.dist, nrows=K, ncols=K, Z=V, startcol=1) 


#### Combine with the covariate matrix if needed
X.anisotropic.exp <- cbind(X, B.anisotropic.exp)
X.anisotropic.inv <- cbind(X, B.anisotropic.inv)
X.anisotropic.linear <- cbind(X, B.anisotropic.linear)


#### Remove an intercept term if it present
    if(var(X.anisotropic.exp[ ,1])==0)
    {
    X.anisotropic.exp <- X.anisotropic.exp[ ,-1]
    X.anisotropic.inv <- X.anisotropic.inv[ ,-1]
    X.anisotropic.linear <- X.anisotropic.linear[ ,-1]
    p <- p-1    
    }else
    {}


#### Remove rows with missing values for model fitting
Y.train <- Y[which.present]
offset.train <- offset[which.present]
K.train <- length(Y.train)
X.anisotropic.exp.train <- X.anisotropic.exp[which.present, ]  
X.anisotropic.inv.train <- X.anisotropic.inv[which.present, ]
X.anisotropic.linear.train <- X.anisotropic.linear[which.present, ]
W.train <- W[which.present, which.present]
W.list.train <- mat2listw(W.train, style="B")



########################################
#### Fit the models and make predictions
########################################
#### Update the user on the functions progress
    if(verbose) cat("Fitting the model.")


#### Fit the models with the 3 different types of basis functions
penfac <- c(rep(0, p), rep(1,K))
mod.ridge.exp <- glmnet(x=X.anisotropic.exp.train, y=Y.train, offset=offset.train, alpha=0, nlambda=nlambda, penalty.factor = penfac, family = "poisson", intercept=TRUE, standardize=FALSE)
mod.ridge.inv <- glmnet(x=X.anisotropic.inv.train, y=Y.train, offset=offset.train, alpha=0, nlambda=nlambda, penalty.factor = penfac, family = "poisson", intercept=TRUE, standardize=FALSE)
mod.ridge.linear <- glmnet(x=X.anisotropic.linear.train, y=Y.train, offset=offset.train, alpha=0, nlambda=nlambda, penalty.factor = penfac, family = "poisson", intercept=TRUE, standardize=FALSE)


#### Compute the level of residual spatial autocorrelation for each model and lambda value
## Remove 0s from Y for computing residual Moran's I
Ytemp <- Y.train
Ytemp[Y.train==0] <- 0.1

## Exponential model
fits.lp.exp <- predict(object=mod.ridge.exp, newx=X.anisotropic.exp.train, newoffset=offset.train)
m <- ncol(fits.lp.exp)
results.exp <- data.frame(lambda=mod.ridge.exp$lambda, I=rep(NA, m))
    for(j in 1:m)
    {
    resids <- log(Ytemp) - fits.lp.exp[ ,j]
    results.exp$I[j] <- moran.mc(x=resids, listw=W.list.train, zero.policy = TRUE, nsim=1)$statistic
    }
row.exp <- which(abs(results.exp$I)==min(abs(results.exp$I)))
moran.exp <- results.exp$I[row.exp]      

## Inverse model   
fits.lp.inv <- predict(object=mod.ridge.inv, newx=X.anisotropic.inv.train, newoffset=offset.train)
m <- ncol(fits.lp.inv)
results.inv <- data.frame(lambda=mod.ridge.inv$lambda, I=rep(NA, m))
    for(j in 1:m)
    {
    resids <- log(Ytemp) - fits.lp.inv[ ,j]
    results.inv$I[j] <- moran.mc(x=resids, listw=W.list.train, zero.policy = TRUE, nsim=1)$statistic
    }
row.inv <- which(abs(results.inv$I)==min(abs(results.inv$I)))
moran.inv <- results.inv$I[row.inv]
    
## Linear model   
fits.lp.linear <- predict(object=mod.ridge.linear, newx=X.anisotropic.linear.train, newoffset=offset.train)
m <- ncol(fits.lp.linear)
results.linear <- data.frame(lambda=mod.ridge.linear$lambda, I=rep(NA, m))
    for(j in 1:m)
    {
    resids <- log(Ytemp) - fits.lp.linear[ ,j]
    results.linear$I[j] <- moran.mc(x=resids, listw=W.list.train, zero.policy = TRUE, nsim=1)$statistic
    }
row.linear <- which(abs(results.linear$I)==min(abs(results.linear$I)))
moran.linear <- results.linear$I[row.linear]       


#### Choose the final model
moran.all <- abs(c(moran.exp, moran.inv, moran.linear))
model <- which(moran.all == min(moran.all))[1]    
      if(model==1)
      {
      model.string <- c("Likelihood model - Poisson (log link function)", "Spatial structure model - Anistropic exponential distance-decay basis functions")     
      model <- mod.ridge.exp
      row <- row.exp
      X.final <- X.anisotropic.exp
      X.final.train <- X.final[which.present, ] 
      lambda.hat <- results.exp$lambda[row] * K.train
      I <- results.exp$I[row]
      }else if(model==2)
      {
      model.string <- c("Likelihood model - Poisson (log link function)", "Spatial structure model - Anistropic inverse distance-decay basis functions")     
      model <- mod.ridge.inv
      row <- row.inv
      X.final <- X.anisotropic.inv
      X.final.train <- X.final[which.present, ] 
      lambda.hat <- results.inv$lambda[row] * K.train
      I <- results.inv$I[row]
      }else if(model==3)
      {
      model.string <- c("Likelihood model - Poisson (log link function)", "Spatial structure model - Anistropic linear distance-decay basis functions")     
      model <- mod.ridge.linear
      row <- row.linear
      X.final <- X.anisotropic.linear
      X.final.train <- X.final[which.present, ] 
      lambda.hat <- results.linear$lambda[row] * K.train
      I <- results.linear$I[row]
      }else{}



#### Compute the parameter estimate for beta
beta.hat <- c(model$a0[row], model$beta[ ,row])



#####################################  
#### Summarise and return the results
#####################################
#### Update the user on the progress
    if(verbose) cat("\nSummarising results.\n")


#### Compute the final fitted / predicted values and residuals 
lp.all <- as.numeric(beta.hat[1] + X.final %*% beta.hat[-1] + offset)
fitted.values <- exp(lp.all)
response.residuals <- Y - fitted.values
pearson.residuals <- response.residuals / sqrt(fitted.values)
residuals <- data.frame(response=response.residuals, pearson=pearson.residuals)   


#### Format the final X matrix returned
X.extend <- cbind(rep(1, K), X.final)
colnames(X.extend)[1] <- "(Intercept)"
colnames(X.extend)[(p+2):(p+K+1)] <- paste("Basis function", 1:K, sep=" ")



#######################
#### Return the results
#######################
results <- list(beta.hat=beta.hat, sigma2.hat=NA, lambda.hat=lambda.hat, I=I, fitted.values=fitted.values, residuals=residuals, formula=formula, model.string=model.string, X=X.extend, model=model)
    if(verbose)
    {
    b<-proc.time()
    cat("Finished in ", round(b[3]-a[3], 1), "seconds.\n")
    }else
    {}
return(results)
}