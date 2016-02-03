\name{CARBayes-package}
\alias{CARBayes-package}
\alias{CARBayes}
\docType{package}

\title{
Spatial Generalised Linear Mixed Models for Areal Unit Data
}


\description{
Implements a class of spatial generalised linear mixed models for areal unit data, 
with inference in a Bayesian setting using Markov chain Monte Carlo (McMC) 
simulation. The response variable can be binomial, Gaussian or Poisson. The 
spatial autocorrelation is modelled by  a set of random effects, which are assigned 
a conditional autoregressive (CAR) prior distribution. A number of different CAR 
priors are available for the random effects, and full details are given in the 
vignette accompanying this package. The initial creation of this package was supported 
by the Economic and Social Research Council (ESRC) grant RES-000-22-4256, and on-going 
development was supported by the Engineering and Physical Science Research Council 
(EPSRC) grant EP/J017442/1.

Version 4.4 has two main changes from version 4.3. Firstly, the three functions 
S.CARiar(), S.independent() and S.CARleroux() have been merged into S.CARleroux(), as 
the latter is a generalisation of the other two. This has been achieved by adding two
additional arguments to S.CARleroux, fix.rho (logical) and rho (numeric). The old
model S.independent() can be obtained by setting (fix.rho=TRUE, rho=0), which
corresponds to independent random effects. Similarly,  the old model S.CARiar() 
corresponding to the intrinsic CAR model can be obtained by setting (fix.rho=TRUE, 
rho=1). The second change is that the modelfit component of the  fitted model object 
now additionally returns the  Watanabe-Akaike Information Criterion (WAIC) and an 
estimate of the effective number of effective parameters (p.w).
}

\details{
\tabular{ll}{
Package: \tab CARBayes\cr
Type: \tab Package\cr
Version: \tab 4.4\cr
Date: \tab 2016-02-03\cr
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
