### R code from vignette source 'CARBayes.Rnw'

###################################################
### code chunk number 1: CARBayes.Rnw:63-64
###################################################
options(prompt = "R> ")


###################################################
### code chunk number 2: CARBayes.Rnw:423-427
###################################################
library(CARBayesdata)
library(sf)
data(pricedata)
data(GGHB.IZ)


###################################################
### code chunk number 3: CARBayes.Rnw:433-435
###################################################
head(pricedata)
head(GGHB.IZ)


###################################################
### code chunk number 4: CARBayes.Rnw:461-464
###################################################
library(dplyr)
pricedata <- pricedata %>% mutate(logprice = log(pricedata$price))
head(pricedata)


###################################################
### code chunk number 5: CARBayes.Rnw:470-472
###################################################
library(GGally)
ggpairs(data = pricedata, columns = c(8, 3:7))


###################################################
### code chunk number 6: CARBayes.Rnw:485-486
###################################################
pricedata.sf <- merge(x=GGHB.IZ, y=pricedata, by="IZ", all.x=FALSE)


###################################################
### code chunk number 7: CARBayes.Rnw:494-496
###################################################
pricedata.sf <- st_transform(x=pricedata.sf, 
                             crs='+proj=longlat +datum=WGS84 +no_defs')


###################################################
### code chunk number 8: CARBayes.Rnw:501-507
###################################################
library(mapview)
library(RColorBrewer)
map1 <- mapview(pricedata.sf, zcol = "price", col.regions=brewer.pal(9, "YlOrRd"), 
            alpha.regions=0.6, layer.name="Price", lwd=0.5, col="grey90", 
            homebutton=FALSE) 
removeMapJunk(map1, junk = c("zoomControl", "layersControl"))


###################################################
### code chunk number 9: CARBayes.Rnw:523-526
###################################################
form <- logprice~crime+rooms+sales+factor(type) + driveshop
model <- lm(formula=form, data=pricedata.sf)
summary(model)


###################################################
### code chunk number 10: CARBayes.Rnw:531-535
###################################################
library(spdep)
W.nb <- poly2nb(pricedata.sf, row.names = pricedata.sf$IZ)
W.list <- nb2listw(W.nb, style="B")
moran.mc(x=residuals(model), listw=W.list, nsim=1000)


###################################################
### code chunk number 11: CARBayes.Rnw:545-546
###################################################
W <- nb2mat(W.nb, style="B")


###################################################
### code chunk number 12: CARBayes.Rnw:703-707
###################################################
library(CARBayesdata)
library(sf)
data(respiratorydata)
data(GGHB.IZ)


###################################################
### code chunk number 13: CARBayes.Rnw:713-715
###################################################
head(pricedata)
head(GGHB.IZ)


###################################################
### code chunk number 14: CARBayes.Rnw:721-723
###################################################
respiratorydata.sf <- merge(x=GGHB.IZ, y=respiratorydata, by="IZ", all.x=FALSE)
head(respiratorydata.sf)


###################################################
### code chunk number 15: CARBayes.Rnw:728-730
###################################################
respiratorydata.sf <- st_transform(x=respiratorydata.sf, 
                                   crs='+proj=longlat +datum=WGS84 +no_defs')


###################################################
### code chunk number 16: CARBayes.Rnw:737-743
###################################################
library(mapview)
library(RColorBrewer)
map2 <- mapview(respiratorydata.sf, zcol = "SMR", col.regions=brewer.pal(9, "YlOrRd"), 
            alpha.regions=0.6, layer.name="SMR", lwd=0.5, col="grey90", 
            homebutton=FALSE) 
removeMapJunk(map2, junk = c("zoomControl", "layersControl"))


###################################################
### code chunk number 17: CARBayes.Rnw:758-760
###################################################
W.nb <- poly2nb(respiratorydata.sf, row.names = respiratorydata.sf$IZ)
W <- nb2mat(W.nb, style="B")


###################################################
### code chunk number 18: CARBayes.Rnw:775-777
###################################################
income <- respiratorydata.sf$incomedep
Z.incomedep <- as.matrix(dist(income, diag=TRUE, upper=TRUE)) 


