### R code from vignette source 'CARBayesvignette.Rnw'

###################################################
### code chunk number 1: CARBayesvignette.Rnw:194-195
###################################################
library(CARBayes)


###################################################
### code chunk number 2: CARBayesvignette.Rnw:200-204
###################################################
library(shapefiles)
library(sp)
library(maptools)
library(spdep)


###################################################
### code chunk number 3: CARBayesvignette.Rnw:209-210
###################################################
library(CARBayesdata)


###################################################
### code chunk number 4: CARBayesvignette.Rnw:246-249
###################################################
data(lipdata)
data(lipdbf)
data(lipshp)


###################################################
### code chunk number 5: CARBayesvignette.Rnw:254-256
###################################################
lipdbf$dbf <- lipdbf$dbf[ ,c(2,1)]
data.combined <- combine.data.shapefile(data=lipdata, shp=lipshp, dbf=lipdbf)


###################################################
### code chunk number 6: CARBayesvignette.Rnw:262-264
###################################################
W.nb <- poly2nb(data.combined, row.names = rownames(lipdata))
W.mat <- nb2mat(W.nb, style="B")


###################################################
### code chunk number 7: CARBayesvignette.Rnw:281-283
###################################################
data(propertydata.spatial)
propertydata <- propertydata.spatial@data


###################################################
### code chunk number 8: CARBayesvignette.Rnw:327-336
###################################################
northarrow <- list("SpatialPolygonsRescale", layout.north.arrow(), offset = c(220000,647000), 
scale = 4000)
scalebar <- list("SpatialPolygonsRescale", layout.scale.bar(), offset = c(225000,647000), 
scale = 10000, fill=c("transparent","black"))
text1 <- list("sp.text", c(225000,649000), "0")
text2 <- list("sp.text", c(230000,649000), "5000 m")
spplot(propertydata.spatial, c("price"), sp.layout=list(northarrow, scalebar, text1, text2), 
scales=list(draw = TRUE), at=seq(min(propertydata$price)-1, max(propertydata$price)+1, length.out=8), 
col.regions=hsv(0,seq(0.05,1,length.out=7),1), col="transparent")


###################################################
### code chunk number 9: CARBayesvignette.Rnw:351-353
###################################################
propertydata$logprice <- log(propertydata$price)
propertydata$logdriveshop <- log(propertydata$driveshop)


###################################################
### code chunk number 10: CARBayesvignette.Rnw:358-361
###################################################
library(splines)
form <- logprice~ns(crime,3)+rooms+sales+factor(type) + logdriveshop
model <- lm(formula=form, data=propertydata)


###################################################
### code chunk number 11: CARBayesvignette.Rnw:369-373
###################################################
W.nb <- poly2nb(propertydata.spatial, row.names = rownames(propertydata))
W.list <- nb2listw(W.nb, style="B")
resid.model <- residuals(model)
moran.mc(x=resid.model, listw=W.list, nsim=1000)


###################################################
### code chunk number 12: CARBayesvignette.Rnw:394-397
###################################################
W.mat <- nb2mat(W.nb, style="B")
model.spatial <- S.CARleroux(formula=form, data=propertydata, family="gaussian", 
W=W.mat, burnin=1000, n.sample=2000, verbose=FALSE)


###################################################
### code chunk number 13: CARBayesvignette.Rnw:406-407
###################################################
summary(model.spatial)


###################################################
### code chunk number 14: CARBayesvignette.Rnw:417-418
###################################################
summarise.samples(model.spatial$samples$beta, quantiles=c(0.5, 0.025, 0.975))


###################################################
### code chunk number 15: CARBayesvignette.Rnw:423-425
###################################################
crime.effect <- summarise.lincomb(model=model.spatial, columns=c(2,3,4), 
quantiles=c(0.5, 0.025, 0.975), distribution=FALSE)


###################################################
### code chunk number 16: CARBayesvignette.Rnw:445-450
###################################################
ord <- order(propertydata$crime)
plot(propertydata$crime[ord], crime.effect$quantiles[ord ,1], pch=19, ylim=c(-0.55,0.05), xlab="Crime rate", ylab="Effect of crime", type="l")
rug(propertydata$crime)
lines(propertydata$crime[ord], crime.effect$quantiles[ord ,2], pch=19, col="red")
lines(propertydata$crime[ord], crime.effect$quantiles[ord ,3], pch=19, col="red")


###################################################
### code chunk number 17: CARBayesvignette.Rnw:476-479
###################################################
data(respiratorydata.spatial)
respiratorydata <- respiratorydata.spatial@data
head(respiratorydata)


###################################################
### code chunk number 18: CARBayesvignette.Rnw:485-489
###################################################
respiratorydata$SIR2010 <- respiratorydata$observed2010 / respiratorydata$expected2010
respiratorydata.spatial@data$SIR2010 <- respiratorydata$SIR2010
W.nb <- poly2nb(respiratorydata.spatial, row.names = rownames(respiratorydata))
W.mat <- nb2mat(W.nb, style="B")


###################################################
### code chunk number 19: CARBayesvignette.Rnw:497-504
###################################################
northarrow <- list("SpatialPolygonsRescale", layout.north.arrow(), offset = c(220000,647000), scale = 4000)
scalebar <- list("SpatialPolygonsRescale", layout.scale.bar(), offset = c(225000,647000), scale = 10000, fill=c("transparent","black"))
text1 <- list("sp.text", c(225000,649000), "0")
text2 <- list("sp.text", c(230000,649000), "5000 m")
spplot(respiratorydata.spatial, c("SIR2010"), sp.layout=list(northarrow, scalebar, text1, text2), 
       scales=list(draw = TRUE), at=seq(min(respiratorydata$SIR2010)-0.05, max(respiratorydata$SIR2010)+0.05, length.out=8), 
       col.regions=hsv(0,seq(0.05,1,length.out=7),1), col="transparent")


###################################################
### code chunk number 20: CARBayesvignette.Rnw:523-526
###################################################
Z.incomedep <- as.matrix(dist(cbind(respiratorydata$incomedep2010, 
respiratorydata$incomedep2010), method="manhattan", diag=TRUE, upper=TRUE)) * 
W.mat / 2


