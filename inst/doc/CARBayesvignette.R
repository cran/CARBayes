### R code from vignette source 'CARBayesvignette.Rnw'

###################################################
### code chunk number 1: CARBayesvignette.Rnw:186-187
###################################################
library(CARBayes)


###################################################
### code chunk number 2: CARBayesvignette.Rnw:192-196
###################################################
library(shapefiles)
library(sp)
library(maptools)
library(spdep)


###################################################
### code chunk number 3: CARBayesvignette.Rnw:245-247
###################################################
data(spatialhousedata)
housedata <- spatialhousedata@data


###################################################
### code chunk number 4: CARBayesvignette.Rnw:297-306
###################################################
northarrow <- list("SpatialPolygonsRescale", layout.north.arrow(), offset = c(220000,647000), 
scale = 4000)
scalebar <- list("SpatialPolygonsRescale", layout.scale.bar(), offset = c(225000,647000), 
scale = 10000, fill=c("transparent","black"))
text1 <- list("sp.text", c(225000,649000), "0")
text2 <- list("sp.text", c(230000,649000), "5000 m")
spplot(spatialhousedata, c("price"), sp.layout=list(northarrow, scalebar, text1, text2), 
scales=list(draw = TRUE), at=seq(min(housedata$price)-1, max(housedata$price)+1, length.out=8), 
col.regions=hsv(0,seq(0.05,1,length.out=7),1), col="transparent")


###################################################
### code chunk number 5: CARBayesvignette.Rnw:321-323
###################################################
housedata$logprice <- log(housedata$price)
housedata$logdriveshop <- log(housedata$driveshop)


###################################################
### code chunk number 6: CARBayesvignette.Rnw:328-331
###################################################
library(splines)
form <- logprice~ns(crime,3)+rooms+sales+factor(type) + logdriveshop
model <- lm(formula=form, data=housedata)


###################################################
### code chunk number 7: CARBayesvignette.Rnw:339-343
###################################################
W.nb <- poly2nb(spatialhousedata, row.names = rownames(housedata))
W.list <- nb2listw(W.nb, style="B")
resid.model <- residuals(model)
moran.mc(x=resid.model, listw=W.list, nsim=1000)


###################################################
### code chunk number 8: CARBayesvignette.Rnw:356-360
###################################################
W.mat <- nb2mat(W.nb, style="B")
model.spatial <- S.CARleroux(formula=form, data=housedata, family="gaussian", 
W=W.mat, burnin=1000, n.sample=2000, verbose=FALSE)
print(model.spatial)


###################################################
### code chunk number 9: CARBayesvignette.Rnw:377-378
###################################################
summary(model.spatial)


###################################################
### code chunk number 10: CARBayesvignette.Rnw:388-389
###################################################
summarise.samples(model.spatial$samples$beta, quantiles=c(0.5, 0.025, 0.975))


###################################################
### code chunk number 11: CARBayesvignette.Rnw:394-396
###################################################
crime.effect <- summarise.lincomb(model=model.spatial, columns=c(2,3,4), 
quantiles=c(0.5, 0.025, 0.975), distribution=FALSE)


###################################################
### code chunk number 12: CARBayesvignette.Rnw:414-417
###################################################
plot(housedata$crime, crime.effect$quantiles[ ,1], pch=19, ylim=c(-0.55,0.05), xlab="Crime rate", ylab="Effect of crime")
points(housedata$crime, crime.effect$quantiles[ ,2], pch=19, col="red")
points(housedata$crime, crime.effect$quantiles[ ,3], pch=19, col="red")


###################################################
### code chunk number 13: CARBayesvignette.Rnw:439-442
###################################################
data(spatialrespdata)
respdata <- spatialrespdata@data
head(respdata)


###################################################
### code chunk number 14: CARBayesvignette.Rnw:448-452
###################################################
respdata$SIR2010 <- respdata$observed2010 / respdata$expected2010
spatialrespdata@data$SIR2010 <- respdata$SIR2010
W.nb <- poly2nb(spatialrespdata, row.names = rownames(respdata))
W.mat <- nb2mat(W.nb, style="B")


###################################################
### code chunk number 15: CARBayesvignette.Rnw:460-467
###################################################
northarrow <- list("SpatialPolygonsRescale", layout.north.arrow(), offset = c(220000,647000), scale = 4000)
scalebar <- list("SpatialPolygonsRescale", layout.scale.bar(), offset = c(225000,647000), scale = 10000, fill=c("transparent","black"))
text1 <- list("sp.text", c(225000,649000), "0")
text2 <- list("sp.text", c(230000,649000), "5000 m")
spplot(spatialrespdata, c("SIR2010"), sp.layout=list(northarrow, scalebar, text1, text2), 
       scales=list(draw = TRUE), at=seq(min(respdata$SIR2010)-0.05, max(respdata$SIR2010)+0.05, length.out=8), 
       col.regions=hsv(0,seq(0.05,1,length.out=7),1), col="transparent")


###################################################
### code chunk number 16: CARBayesvignette.Rnw:486-488
###################################################
Z.incomedep <- as.matrix(dist(cbind(respdata$incomedep2010, respdata$incomedep2010), 
method="manhattan", diag=TRUE, upper=TRUE)) * W.mat / 2


