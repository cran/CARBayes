\name{CARBayes-package}
\alias{CARBayes-package}
\alias{CARBayes}
\docType{package}

\title{
Spatial Areal Unit Modelling
}


\description{
Implements Bayesian hierarchical spatial areal unit models. In such models the spatial correlation is modelled by  a set of random effects, which are assigned a conditional autoregressive (CAR) prior distribution. Examples of the models included are the BYM model as well as a recently developed localised spatial smoothing model. In addition, the package includes functions to enable plotting spatial maps of data (by combining the data with an appropriate shapefile) using the functionality of the sp package. There is a vignette accompanying this package. The initial creation of this package was supported by the Economic and Social Research Council (ESRC) grant RES-000-22-4256, and on-going development is supported by the Engineering and Physical Science Research Council (EPSRC) grant EP/J017442/1.
}

\details{
\tabular{ll}{
Package: \tab CARBayes\cr
Type: \tab Package\cr
Version: \tab 4.0\cr
Date: \tab 2014-12-15\cr
License: \tab GPL (>= 2)\cr
}

Implements Bayesian hierarchical models for spatial areal unit data. The models are an extension of generalised linear models, and the response variable can be binomial, Gaussian or Poisson. The mean function for the data comprises a covariate component and a set of random effects. This package implements a number of different random effects models, which represent different types of spatial autocorrelation structure. These models are within the conditional autoregressive (CAR) class, and capture spatial correlation via a neighbourhood matrix W. This matrix typically (although other options are allowed) has a jkth element equal to one and zero, depending on whether areas j and k are spatially close (typically chosen to mean sharing a common border). Inference for all models is via Markov Chain Monte Carlo (MCMC) simulation. The following random effects models are implemented in this package.

Independent and identically distributed random effects.

Spatially correlated random effects modelled by the intrinsic autoregressive (IAR) model proposed by Besag et al (1991).

Spatially correlated random effects modelled by the Besag-York-Mollie (BYM) model proposed by Besag et al (1991).

Spatially correlated random effects modelled by the proper CAR model proposed by Leroux (1999).

Spatially correlated random effects modelled by the localised spatial smoothing approach proposed by Lee and Mitchell (2012).

Spatially correlated random effects modelled by the localised spatial smoothing approach proposed by Lee and Sarran (2014).

More detailed descriptions of the methods are given by the references below. CARBayes has several changes from version 3, one of which is the availability of an up-to-date vignette. Other changes are the removal of the proper CAR model by Stern et al (1999), and a change in the notation of the main functions from, for example, lerouxCAR.re() to S.CARleroux(). Finally, two additional functions have been added, which allow the posterior samples (summarise.samples()) and linear combinations of the covariate component (summarise.lincomb()), for example to summarise non-linear relationships, to be summarised. 
}




\author{
Maintainer: Duncan Lee <Duncan.Lee@glasgow.ac.uk>
}


\references{
Besag, J. and  York, J and Mollie, A (1991). Bayesian image restoration with two applications in spatial statistics. Annals of the Institute of Statistics and Mathematics 43, 1-59.

Eddelbuettel, D and Francois, R (2011). Rcpp: seamless R and C++ Integration. Journal of Statistical Software, 40, 1-18. 

Stern, H and Cressie, N. (1999). Inference for extremes in disease mapping, Chapter
Disease mapping and Risk Assessment for Public Health. Lawson, A and Biggeri, D
and Boehning, E and Lesaffre, E and Viel, J and Bertollini, R (eds). Wiley.

Lee, D. (2011). A comparison of conditional autoregressive model used in Bayesian disease mapping. Spatial and Spatio-temporal Epidemiology 2, 79-89.

Lee, D. and Mitchell, R (2012). Boundary detection in disease mapping studies. Biostatistics, 13, 415-426.

Lee, D (2013). CARBayes: An R Package for Bayesian Spatial Modeling with Conditional Autoregressive Priors. Journal of Statistical Software, 55, 13, 1-24.

Leroux, B., Lei, X and Breslow, N (1999). Estimation of disease rates in small areas: A new mixed model for spatial dependence, Chapter Statistical Models in Epidemiology, the Environment and Clinical Trials, Halloran, M and Berry, D (eds), pp. 135-178. Springer-Verlag, New York.

Lee, D and Sarran, C  (2014). Controlling for unmeasured confounding and spatial misalignment in long-term air pollution and health studies. arXiv:1412.4479.
}



\examples{
library(shapefiles)
library(sp)
library(maptools)
library(splines)
library(spdep)
library(CARBayes)


#### Read in the data
data(spatialhousedata)
housedata <- spatialhousedata@data

#### Transform the price and driveshop variables
housedata$logprice <- log(housedata$price)
housedata$logdriveshop <- log(housedata$driveshop)

#### Fit a covariate only  model
form <- logprice~ns(crime,3)+rooms+sales+factor(type) + logdriveshop
model <- lm(formula=form, data=housedata)


#### Compute a Moran's I test on the residuals
W.nb <- poly2nb(spatialhousedata, row.names = rownames(housedata))
W.list <- nb2listw(W.nb, style="B")
resid.model <- residuals(model)
moran.mc(x=resid.model, listw=W.list, nsim=1000)


#### Run a regression model with spatially correlated random effects
W.mat <- nb2mat(W.nb, style="B")
\dontrun{
model.spatial <- S.CARleroux(formula=form, data=housedata, 
family="gaussian", W=W.mat, burnin=10000, n.sample=40000, 
verbose=FALSE)
print(model.spatial)
}
}
