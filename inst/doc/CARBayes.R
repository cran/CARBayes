### R code from vignette source 'CARBayes.Rnw'

###################################################
### code chunk number 1: CARBayes.Rnw:63-64
###################################################
options(prompt = "R> ")


###################################################
### code chunk number 2: CARBayes.Rnw:419-423
###################################################
library(CARBayesdata)
library(sf)
data(pricedata)
data(GGHB.IZ)


###################################################
### code chunk number 3: CARBayes.Rnw:429-431
###################################################
head(pricedata)
head(GGHB.IZ)


###################################################
### code chunk number 4: CARBayes.Rnw:457-460
###################################################
library(dplyr)
pricedata <- pricedata %>% mutate(logprice = log(pricedata$price))
head(pricedata)


###################################################
### code chunk number 5: CARBayes.Rnw:466-468
###################################################
library(GGally)
ggpairs(data = pricedata, columns = c(8, 3:7))


###################################################
### code chunk number 6: CARBayes.Rnw:481-482
###################################################
pricedata.sf <- merge(x=GGHB.IZ, y=pricedata, by="IZ", all.x=FALSE)


###################################################
### code chunk number 7: CARBayes.Rnw:490-492
###################################################
pricedata.sf <- st_transform(x=pricedata.sf, 
                             crs='+proj=longlat +datum=WGS84 +no_defs')


###################################################
### code chunk number 8: CARBayes.Rnw:497-507
###################################################
library(leaflet)
colours <- colorNumeric(palette = "YlOrRd", domain = pricedata.sf$price)
map1 <- leaflet(data=pricedata.sf) %>% 
    addTiles() %>% 
    addPolygons(fillColor = ~colours(price), color="", weight=1, 
                fillOpacity = 0.7) %>%
    addLegend(pal = colours, values = pricedata.sf$price, opacity = 1, 
              title="Price") %>%
    addScaleBar(position="bottomleft")
map1


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
### code chunk number 12: CARBayes.Rnw:729-733
###################################################
library(CARBayesdata)
library(sf)
data(respiratorydata)
data(GGHB.IZ)


###################################################
### code chunk number 13: CARBayes.Rnw:739-741
###################################################
head(pricedata)
head(GGHB.IZ)


###################################################
### code chunk number 14: CARBayes.Rnw:747-749
###################################################
respiratorydata.sf <- merge(x=GGHB.IZ, y=respiratorydata, by="IZ", all.x=FALSE)
head(respiratorydata.sf)


###################################################
### code chunk number 15: CARBayes.Rnw:754-756
###################################################
respiratorydata.sf <- st_transform(x=respiratorydata.sf, 
                                   crs='+proj=longlat +datum=WGS84 +no_defs')


###################################################
### code chunk number 16: CARBayes.Rnw:763-773
###################################################
library(leaflet)
colours <- colorNumeric(palette = "YlOrRd", domain = respiratorydata.sf$SMR)
map2 <- leaflet(data=respiratorydata.sf) %>% 
  addTiles() %>% 
  addPolygons(fillColor = ~colours(SMR), color="", weight=1, 
              fillOpacity = 0.7) %>%
  addLegend(pal = colours, values = respiratorydata.sf$SMR, opacity = 1, 
            title="SMR") %>%
  addScaleBar(position="bottomleft")
map2


###################################################
### code chunk number 17: CARBayes.Rnw:788-790
###################################################
W.nb <- poly2nb(respiratorydata.sf, row.names = respiratorydata.sf$IZ)
W <- nb2mat(W.nb, style="B")


###################################################
### code chunk number 18: CARBayes.Rnw:805-807
###################################################
income <- respiratorydata.sf$incomedep
Z.incomedep <- as.matrix(dist(income, diag=TRUE, upper=TRUE)) 


