### R code from vignette source 'CARBayes.Rnw'

###################################################
### code chunk number 1: CARBayes.Rnw:62-63
###################################################
options(prompt = "R> ")


###################################################
### code chunk number 2: CARBayes.Rnw:380-386
###################################################
library(CARBayesdata)
library(shapefiles)
library(sp)
data(lipdata)
data(lipdbf)
data(lipshp)


###################################################
### code chunk number 3: CARBayes.Rnw:391-394
###################################################
library(CARBayes)
lipdbf$dbf <- lipdbf$dbf[ ,c(2,1)]
data.combined <- combine.data.shapefile(data=lipdata, shp=lipshp, dbf=lipdbf)


###################################################
### code chunk number 4: CARBayes.Rnw:400-403
###################################################
library(spdep)
W.nb <- poly2nb(data.combined, row.names = rownames(lipdata))
W.mat <- nb2mat(W.nb, style="B")


###################################################
### code chunk number 5: CARBayes.Rnw:420-424
###################################################
library(CARBayesdata)
library(sp)
data(GGHB.IG)
data(pricedata)


###################################################
### code chunk number 6: CARBayes.Rnw:429-433
###################################################
missing.IG <- setdiff(rownames(GGHB.IG@data), pricedata$IG)
missing.IG.row <- which(missing.IG==rownames(GGHB.IG@data))
propertydata.spatial <- GGHB.IG[-missing.IG.row, ]
propertydata.spatial@data <- data.frame(propertydata.spatial@data, pricedata)


###################################################
### code chunk number 7: CARBayes.Rnw:461-471
###################################################
northarrow <- list("SpatialPolygonsRescale", layout.north.arrow(), 
    offset = c(220000,647000), scale = 4000)
scalebar <- list("SpatialPolygonsRescale", layout.scale.bar(), 
    offset = c(225000,647000), scale = 10000, fill=c("transparent","black"))
text1 <- list("sp.text", c(225000,649000), "0")
text2 <- list("sp.text", c(230000,649000), "5000 m")
breakpoints <- seq(min(pricedata$price)-1, max(pricedata$price)+1, length.out=8)
spplot(propertydata.spatial, c("price"), sp.layout=list(northarrow, scalebar, 
    text1, text2), scales=list(draw = TRUE), at=breakpoints,
    col.regions=terrain.colors(n=length(breakpoints)-1), col="transparent")


###################################################
### code chunk number 8: CARBayes.Rnw:488-490
###################################################
propertydata.spatial@data$logprice <- log(propertydata.spatial@data$price)
propertydata.spatial@data$logdriveshop <- log(propertydata.spatial@data$driveshop)


###################################################
### code chunk number 9: CARBayes.Rnw:495-498
###################################################
library(splines)
form <- logprice~ns(crime,3)+rooms+sales+factor(type) + logdriveshop
model <- lm(formula=form, data=propertydata.spatial@data)


###################################################
### code chunk number 10: CARBayes.Rnw:505-510
###################################################
library(spdep)
W.nb <- poly2nb(propertydata.spatial, row.names = rownames(propertydata.spatial@data))
W.list <- nb2listw(W.nb, style="B")
resid.model <- residuals(model)
moran.mc(x=resid.model, listw=W.list, nsim=1000)


###################################################
### code chunk number 11: CARBayes.Rnw:654-658
###################################################
library(CARBayesdata)
library(sp)
data(GGHB.IG)
data(respiratorydata)


###################################################
### code chunk number 12: CARBayes.Rnw:663-672
###################################################
missing.IG <- setdiff(rownames(GGHB.IG@data), respiratorydata$IG)
missing.IG.row <- rep(NA, length(missing.IG))
    for(i in 1:length(missing.IG))
    {
    missing.IG.row[i] <- which(missing.IG[i] == rownames(GGHB.IG@data))     
    }
respiratorydata.spatial <- GGHB.IG[-missing.IG.row, ]
respiratorydata.spatial@data <- data.frame(respiratorydata.spatial@data, 
    respiratorydata)


###################################################
### code chunk number 13: CARBayes.Rnw:677-678
###################################################
head(respiratorydata.spatial@data)


###################################################
### code chunk number 14: CARBayes.Rnw:685-696
###################################################
northarrow <- list("SpatialPolygonsRescale", layout.north.arrow(), offset = 
    c(220000,647000), scale = 4000)
scalebar <- list("SpatialPolygonsRescale", layout.scale.bar(), offset = 
    c(225000,647000), scale = 10000, fill=c("transparent","black"))
text1 <- list("sp.text", c(225000,649000), "0")
text2 <- list("sp.text", c(230000,649000), "5000 m")
breakpoints <- seq(min(respiratorydata.spatial@data$SMR)-0.05, 
    max(respiratorydata.spatial@data$SMR)+0.05, length.out=8)
spplot(respiratorydata.spatial, c("SMR"), sp.layout=list(northarrow, scalebar, 
    text1, text2), scales=list(draw = TRUE), at=breakpoints, 
    col.regions=terrain.colors(n=length(breakpoints)-1), col="transparent")


###################################################
### code chunk number 15: CARBayes.Rnw:711-714
###################################################
W.nb <- poly2nb(respiratorydata.spatial, row.names = 
    rownames(respiratorydata.spatial@data))
W <- nb2mat(W.nb, style="B")


###################################################
### code chunk number 16: CARBayes.Rnw:722-724
###################################################
income <- respiratorydata.spatial@data$incomedep
Z.incomedep <- as.matrix(dist(income, diag=TRUE, upper=TRUE)) 


