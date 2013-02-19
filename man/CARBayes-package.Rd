\name{CARBayes-package}
\alias{CARBayes-package}
\alias{CARBayes}
\docType{package}

\title{
Spatial areal unit modelling
}

\description{
This package implements Bayesian hierarchical spatial areal unit models. In such models the spatial correlation is modelled by  a set of random effects, which are assigned a conditional autoregressive (CAR) prior distribution. Examples of the models included are the BYM model as well as a recently developed localised spatial smoothing model. In addition, the package includes functions to enable plotting spatial maps of data (by combining the data with an appropriate shapefile) using the functionality of the sp package.

The creation of this package was supported by the Economic and Social Research Council (ESRC) grant RES-000-22-4256.
}

\details{
\tabular{ll}{
Package: \tab CARBayes\cr
Type: \tab Package\cr
Version: \tab 1.3\cr
Date: \tab 2013-02-19\cr
License: \tab GPL (>= 2)\cr
}

The package implements Bayesian hierarchical models for spatial areal unit data. The models are an extension of generalised linear models, and the response variable can be binomial, Gaussian or Poisson. The mean function for the data comprises a covariate component and a set of random effects. This package implements a number of different random effects models, which represent different types of spatial correlation structure. These models are within the conditional autoregressive (CAR) class, and capture spatial correlation via a neighbourhood matrix W. This matrix typically (although other options are allowed) has a jkth element equal to one and zero, depending on whether areas j and k are spatially close (typically chosen to mean sharing a common border). Inference for all models is via Markov Chain Monte Carlo (MCMC) simulation. The following random effects models are implemented in this package.

Independent and identically distributed random effects.

Spatially correlated random effects modelled by the intrinsic autoregressive (IAR) model proposed by Besag et al (1991).

Spatially correlated random effects modelled by the Besag-York-Mollie (BYM) model proposed by Besag et al (1991).

Spatially correlated random effects modelled by the proper CAR model proposed by Stern and Cressie (1999).

Spatially correlated random effects modelled by the proper CAR model proposed by Leroux (1999).

Spatially correlated random effects modelled by the localised spatial smoothing approach proposed by Lee and Mitchell (2012).

More detailed descriptions of the methods are given by the references below.
}


\author{
Maintainer: Duncan Lee <Duncan.Lee@glasgow.ac.uk>
}


\references{
Besag, J., J. York, and A. Mollie (1991). Bayesian image restoration with two applications in spatial statistics. Annals of the Institute of Statistics and Mathematics 43, 1-59.

Stern, H and Cressie, N. (1999). Inference for extremes in disease mapping, Chapter
Disease mapping and Risk Assessment for Public Health. Lawson, A and Biggeri, D
and Boehning, E and Lesaffre, E and Viel, J and Bertollini, R (eds). Wiley.

Lee, D. (2011). A comparison of conditional autoregressive model used in Bayesian disease mapping. Spatial and Spatio-temporal Epidemiology 2, 79-89.

Lee, D. and R. Mitchell (2012). Boundary detection in disease mapping studies. Biostatistics, 13, 415-426.

Leroux, B., X. Lei, and N. Breslow (1999). Estimation of disease rates in small areas: A new mixed model for spatial dependence, Chapter Statistical Models in Epidemiology, the Environment and Clinical Trials, Halloran, M and Berry, D (eds), pp. 135-178. Springer-Verlag, New York.
}



\examples{
################################################################
#### Load the libraries required to do an extended analysis
################################################################
library(MASS)
library(lattice)
library(coda)
library(MCMCpack)
library(spam)
library(evd)
library(stats4)
library(truncdist)
library(foreign)
library(shapefiles)
library(sp)
library(boot)
library(Matrix)
library(nlme)
library(maptools)
library(grid)
library(deldir)
library(splines)
library(spdep)
library(CARBayes)


#############################
#### Example 1 - house prices
#############################
#### Read in the data
data(housedata) 
data(shp)
data(dbf)

#### Combine the data and shapefile and create the neighbourhood matrix W
data.combined <- combine.data.shapefile(housedata, shp, dbf)
W.nb <- poly2nb(data.combined, row.names = rownames(housedata))
W.mat <- nb2mat(W.nb, style="B")

#### Run a regression model with spatially correlated random effects
housedata$crime.inv <- 1 / housedata$crime
form <- housedata$price~housedata$crime.inv+housedata$rooms+housedata$sales+factor(housedata$type)
\dontrun{
model.spatial <- gaussian.properCAR(formula=form, W=W.mat, burnin=20000, n.sample=100000, 
thin=10, prior.max.nu2=1000000, prior.max.tau2=1000000, prior.var.beta=rep(10000000, 7))

#### Plot the fitted values on a map
fitted.housedata <- model.spatial$fitted.values[ ,3]
data.combined@data <- data.frame(data.combined@data, fitted.housedata)
northarrow <- list("SpatialPolygonsRescale", layout.north.arrow(), offset = c(220000,647000), 
scale = 4000)
scalebar <- list("SpatialPolygonsRescale", layout.scale.bar(), offset = c(225000,647000), 
scale = 10000, fill=c("transparent","black"))
text1 <- list("sp.text", c(225000,649000), "0")
text2 <- list("sp.text", c(230000,649000), "5000 m")
spplot(data.combined, c("fitted.housedata"), sp.layout=list(northarrow, scalebar, text1, text2), 
at=seq(min(fitted.housedata)-1, max(fitted.housedata)+1, length.out=51), 
col.regions=terrain.colors(n=50))
}


################################
#### Example 2 - disease mapping
################################
#### Read in the data
data(respdata) 
data(shp)
data(dbf)

#### Combine the data and shapefile and create the neighbourhood matrix W
data.combined <- combine.data.shapefile(respdata, shp, dbf)
W.nb <- poly2nb(data.combined, row.names = rownames(respdata))
W.mat <- nb2mat(W.nb, style="B")

#### Create the dissimilrity metric
Z.incomedep <- as.matrix(dist(cbind(respdata$incomedep2010, respdata$incomedep2010), 
method="manhattan", diag=TRUE, upper=TRUE)) * W.mat

#### Run the local CAR model
formula <- respdata$observed2010 ~ offset(log(respdata$expected2010))
\dontrun{
model.dissimilarity <- poisson.dissimilarityCAR(formula=formula, W=W.mat, 
Z=list(Z.incomedep=Z.incomedep), rho=0.99, fix.rho=TRUE, burnin=20000, 
n.sample=100000, thin=10)

#### Plot a map with the boundaries overlayed
border.locations <- model.dissimilarity$W.posterior
risk.estimates <- model.dissimilarity$fitted.values[ ,3] / respdata$expected2010
data.combined@data <- data.frame(data.combined@data, risk.estimates)
boundary.final <- highlight.borders(border.locations=border.locations, ID=rownames(respdata), 
shp=shp, dbf=dbf)
boundaries = list("sp.points", boundary.final, col="white", pch=19, cex=0.2)
northarrow <- list("SpatialPolygonsRescale", layout.north.arrow(), offset = c(220000,647000), 
scale = 4000)
scalebar <- list("SpatialPolygonsRescale", layout.scale.bar(), offset = c(225000,647000), 
scale = 10000, fill=c("transparent","black"))
text1 <- list("sp.text", c(225000,649000), "0")
text2 <- list("sp.text", c(230000,649000), "5000 m")
spplot(data.combined, c("risk.estimates"), sp.layout=list(northarrow, scalebar, text1, text2, boundaries),
at=seq(min(risk.estimates)-0.1, max(risk.estimates)+0.1, length.out=51), 
col.regions=terrain.colors(n=50))
}
}
