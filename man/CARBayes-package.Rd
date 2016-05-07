\name{CARBayes-package}
\alias{CARBayes-package}
\alias{CARBayes}
\docType{package}

\title{
Spatial Generalised Linear Mixed Models for Areal Unit Data
}


\description{
Implements a class of univariate and multivariate spatial generalised linear mixed 
models for areal unit data, with inference in a Bayesian setting using Markov chain 
Monte Carlo (McMC) simulation. The response variable can be binomial, Gaussian or 
Poisson. Spatial autocorrelation is modelled by  a set of random effects, which are 
assigned a conditional autoregressive (CAR) prior distribution. A number of different 
CAR priors are available for the random effects, including a multivariate CAR (MCAR) 
model for multivariate spatial data. Full details are given in the vignette accompanying 
this package. The initial creation of this package was supported by the Economic and 
Social Research Council (ESRC) grant RES-000-22-4256, and on-going development was 
supported by the Engineering and Physical Science Research Council (EPSRC) grant 
EP/J017442/1. Version 4.5 has one main change from version 4.4, namely the introduction 
of an MCAR model for multivariate spatial data, for details see the vignette.
}

\details{
\tabular{ll}{
Package: \tab CARBayes\cr
Type: \tab Package\cr
Version: \tab 4.5\cr
Date: \tab 2016-05-05\cr
License: \tab GPL (>= 2)\cr
}
}




\author{
Maintainer: Duncan Lee <Duncan.Lee@glasgow.ac.uk>
}


\references{
Besag, J. and  York, J and Mollie, A (1991). Bayesian image restoration with two 
applications in spatial statistics. Annals of the Institute of Statistics and 
Mathematics 43, 1-59.

Gelfand, A and Vounatsou, P (2003). Proper multivariate conditional autoregressive
models for spatial data analysis, Biostatistics, 4, 11-25.

Lee, D. and Mitchell, R (2012). Boundary detection in disease mapping studies. 
Biostatistics, 13, 415-426.

Lee, D and Sarran, C  (2015). Controlling for unmeasured confounding and spatial 
misalignment in long-term air pollution and health studies, Environmetrics, 26, 
477-487.

Leroux, B., Lei, X and Breslow, N (1999). Estimation of disease rates in small 
areas: A new mixed model for spatial dependence, Chapter Statistical Models in 
Epidemiology, the Environment and Clinical Trials, Halloran, M and Berry, D (eds), 
pp. 135-178. Springer-Verlag, New York.
}



\examples{
## See the examples in the function specific help files and in the vignette
## accompanying this package.
}
