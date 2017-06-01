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
Monte Carlo (MCMC) simulation. The response variable can be binomial, Gaussian or 
Poisson, and spatial autocorrelation is modelled by  a set of random effects that are 
assigned a conditional autoregressive (CAR) prior distribution. A number of different 
models are available for univariate spatial data, including models with no random 
effects as well as random effects modelled by different types of CAR prior. Additionally,  
a multivariate CAR (MCAR) model for multivariate spatial data is available, as is
a two-level hierarchical model for individuals within areas. Full details are given in 
the vignette accompanying this package. The initial creation of this package was 
supported by the Economic and Social Research Council (ESRC) grant RES-000-22-4256, 
and on-going development has / is supported by the Engineering and Physical Science 
Research Council (EPSRC) grant EP/J017442/1, ESRC grant ES/K006460/1, and 
Innovate UK / Natural Environment Research Council (NERC) grant NE/N007352/1. 
Version 5.0 has a  number of changes from version 4.7, including: (i) inclusion of a
simple model with no random effects (S.CARglm()); (ii) inclusion of a two-level 
hierarchical model for modelling individuals within areas (S.CARmultilevel()); 
(iii) an argument in the model fitting functions to allow the random effects to be 
updated using either random walk Metropolis steps or the  Metropolis adjusted Langevin 
algorithm (MALA); and (iv) calculation of the percentage of the deviance explained 
by the model.
}

\details{
\tabular{ll}{
Package: \tab CARBayes\cr
Type: \tab Package\cr
Version: \tab 5.0\cr
Date: \tab 2017-06-01\cr
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

Kavanagh, L., D. Lee, and G. Pryce (2016). Is Poverty Decentralising? Quantifying 
Uncertainty in the Decentralisation of Urban Poverty, Annals of the American 
Association of Geographers, 106, 1286-1298.

Lee, D. and Mitchell, R (2012). Boundary detection in disease mapping studies. 
Biostatistics, 13, 415-426.

Lee, D and Sarran, C  (2015). Controlling for unmeasured confounding and spatial 
misalignment in long-term air pollution and health studies, Environmetrics, 26, 
477-487.

Leroux, B., Lei, X and Breslow, N (2000). Estimation of disease rates in small 
areas: A new mixed model for spatial dependence, Chapter Statistical Models in 
Epidemiology, the Environment and Clinical Trials, Halloran, M and Berry, D (eds), 
pp. 135-178. Springer-Verlag, New York.

Roberts, G and Rosenthal, J (1998). Optimal scaling of discrete approximations to
the Langevin diffusions, Journal of the Royal Statistical Society Series B 60, 
255-268.
}



\examples{
## See the examples in the function specific help files and in the vignette
## accompanying this package.
}
