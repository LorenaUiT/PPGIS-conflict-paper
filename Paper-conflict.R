### Libraries ####
library(reshape2)
library(raster)
library(sp)
library(rgdal)
library(rgeos)
library(grid)
library(lattice) 
library(mvabund)  ###glm
library(spatstat)
library(maptools) ##to create window to use in spatstat

#options(stringsAsFactors=FALSE) # turn off automatic factor coersion
options(scipen=9999)            # turn off plotting axis lables in scientific notation
options(help_type="html")       # help as html
setwd("C:/Users/lma021/Desktop/Conflict paper")

#### Table preparation ####

### Load datasets
first<-read.delim("First-conflict.txt")
markers<-read.delim("Markers-conflict.txt")
second<-read.delim("Second-conflict.txt")
second<-subset(second,second$age>17) ##select respondents that are above 18
first<-first[!duplicated(first$LoginID),] ## eliminates duplicates of LoginID (each respondent should have a unique questionnaire)

## Select desired columns
first<-subset(first[,c(4,9,11,15,21)])
markers<-subset(markers[,c(4,10,11,15,16,17)])
second<-subset(second[,c(1,69:71,77,85)])

markers <- markers[!(markers$category2conflict == ""), ] ## delete rows that have no markers for the conflict analysis

## Merging tables:
firstsecond<-merge(first,second,by.x = "LoginID",by.y = "LoginID")
markerstype<-dcast(markers,LoginID~category2conflict)
firstsecondmarkerstype<-merge(firstsecond,markerstype,by.x = "LoginID",by.y = "LoginID")
markersfirstsecond<-merge(markers,firstsecond,by.x = "LoginID",by.y = "LoginID")

### No need to repeat this again ####
## create a table of markers and create a shp with them
write.csv(markersfirstsecond, "markersfirstsecond.csv")
## Be sure that the projection is set to WGS 84


#### Prepare point and buffer layers ####

#Create a shapefile out of the points 
pointES<-readOGR(dsn="markersfirstsecond.shp",layer="markersfirstsecond")
plot(pointES)
pointES@proj4string ## to check for the projection and the units (meters in this case)


## Read the PA layer and create buffer zones of 100m around each PA
PA<-readOGR(dsn="PA polygons.shp", layer="PA polygons")
PA<-SpatialPolygonsDataFrame(PA,data=PA@data)
plot(PA)
PA@proj4string ## to check for the projection and the units (meters in this case)
PAbuffer<-buffer(PA, width =1000,dissolve=FALSE)
PAbuffer<-SpatialPolygonsDataFrame(PAbuffer,data=PA@data)
plot(PAbuffer)
PAbuffer@proj4string
PAbufferdis<-unionSpatialPolygons(PAbuffer,PAbuffer@data$OBJTYPE) ### creates a unique polygon for the whole study area to be used with "ppp" later

## Clip ES points with PAbuffer
pointES$PA<-over(pointES,PAbuffer)$NAVN
plot(PAbufferdis, border="black")
plot(pointES,col=factor(pointES$PA),add=TRUE)
pointES <- subset(pointES, PA!=" ") ## delete rows that have no markers inside PAs
summary(pointES)
plot(pointES,col=factor(pointES$origin),add=TRUE)
plot(pointES,col=factor(pointES$category2c),add=TRUE)
plot(pointES,col=factor(pointES$type),add=TRUE)
pointES.values<-subset(pointES,pointES@data$type=="placevalues")
pointES.pref<-subset(pointES,!is.na(pointES@data$preference))

#### GLM ####
#### Check whether there are differences in place values and preferences between user groups
## glm using mvabund library
# values
pointES.values.table<-dcast(pointES.values@data,origin~category2c)
pointES.values.table<-subset(pointES.values.table[-4,]) ## Deleting the last row for origin=NA
values<-mvabund(pointES.values.table[,2:13])
origin<-pointES.values.table$origin
values.nb <- manyglm(values ~ origin, family="negative.binomial") 
anova(values.nb,p.uni="adjusted") ### Origin does not explain the variance in place values
plot(values.nb)
coef(values.nb)
residuals(values.nb)
values.mvformula<-mvformula(values~origin)
plot(values.mvformula)

# preferences
pointES.pref.table<-dcast(pointES@data,origin~category2c[type=="preferences"])
pointES.pref.table<-subset(pointES.pref.table[-4,]) ## Deleting the last row for origin=NA
preferences<-mvabund(pointES.pref.table[,2:8])
origin<-pointES.pref.table$origin
pref.nb <- manyglm(preferences ~ origin, family="negative.binomial") 
anova(pref.nb,p.uni="adjusted") ### Origin does not explain the variance in place values
plot(pref.nb)
preferences.mvformula<-mvformula(preferences~origin)
plot(preferences.mvformula)


####
##### PLACE VALUES #####
####
#### Create polygons with points that are nearby for each user group 
window<-as.owin(PAbufferdis) ## we need to create an observation window to delimit the study area
m<-pointES.values@data[,c(1,6,11)] # We put origin as marks, so we are able to identify areas of different user groups
pointES.values.ppp<-ppp(x=pointES.values@coords[,1],y=pointES.values@coords[,2], window = window,marks = m)
pointES.values.ppp$n/sum(sapply(slot(PAbufferdis, "polygons"), slot, "area")) ### Calculates average intensity of markers
plot(pointES.values.ppp)
plot(Kest(pointES.values.ppp)) ## How do I interpret this???

## Errors: there are duplicated points in pointES.values that appear if marks="m" is not included in the model
# solved vy including field_1, which are unique identifiers (I think)

plot(density(pointES.values.ppp,1000))
persp(density(pointES.values.ppp,1000))
plot(density(pointES.values.ppp[pointES.values.ppp$marks$origin=="Local"],1000))
plot(density(pointES.values.ppp[pointES.values.ppp$marks$origin=="Domestic"],1000))
plot(density(pointES.values.ppp[pointES.values.ppp$marks$origin=="International"],1000))

## Kolmogorov-Smirnov test to compare expected and observed distributions of x and y
# kstest is out of date, R recomends using cdf.test instead
cdf.test(pointES.values.ppp, function(x, y) {x}) ## Rejects H0, therefore expected and observed distributions of x are different
plot(cdf.test(pointES.values.ppp, function(x, y) {x}))
cdf.test(pointES.values.ppp, function(x, y) {y}) ## Rejects H0, therefore expected and observed distributions of y are different
plot(cdf.test(pointES.values.ppp, function(x, y) {y}))

## We cannot use MLE for poisson because our points are not completely independent (is this right?)
## MLE tests covariates?

### Nearest neighbour distances
ave(Gest(pointES.values.ppp[pointES.values.ppp$marks$origin=="Local"])$r) #Result: 1523.188
plot(Gest(pointES.values.ppp[pointES.values.ppp$marks$origin=="Local"])) #Cumulative distribution function of the nearest-neighbour distance
ave(Gest(pointES.values.ppp[pointES.values.ppp$marks$origin=="Domestic"])$r) #Result: 1326.815
plot(Gest(pointES.values.ppp[pointES.values.ppp$marks$origin=="Domestic"]))
ave(Gest(pointES.values.ppp[pointES.values.ppp$marks$origin=="International"])$r) #Result: 1295.643
plot(Gest(pointES.values.ppp[pointES.values.ppp$marks$origin=="International"]))
# Place values by international tourists are the ones with the shorterst avegare distance between neighbouring points
# This may mean that international tourists have a lower dispersion than Locals (highest average distance between points)
# Plots show that bG(r) > Gpois(r) suggest a clustered pattern.
# This can be concluded from the plots as they show a rapid increase on the accummulated distance, which means that points are very close together





#### STOPPED HERE####



#### MAPS ####
## Create a grid to be used for preference analysis
grid<-raster(extent(PAbuffer)) ## create a raster
proj4string(grid)<-proj4string(PAbuffer)## put the same resolution as the buffer layer
res(grid)<-2000 ## choose resolution
gridpolygon<-rasterToPolygons(grid) ##transform the raster into a polygon
gridbuffer<-intersect(PAbuffer,gridpolygon)
plot(gridbuffer)


#### Create a map for development index ####

function.dev.index <- function(currorigin, preference.positive, preference.negative, grid.size){
  pointES.pref.origin <- subset(pointES.pref, origin==currorigin)
  pos <- subset(pointES.pref.origin, preference==preference.positive)
  neg <- subset(pointES.pref.origin, preference==preference.negative)
  res(grid)<-grid.size
  posRast <- rasterize(pos, grid, fun="count")
  negRast <- rasterize(neg, grid, fun="count")
  outRast <- (posRast-negRast)*(posRast+negRast)
  writeRaster(outRast, filename=paste0(currorigin, "_" , preference.positive),overwrite=TRUE)
  return(outRast)
}

devindexRastLocal <- function.dev.index("Local", "development", "nodevelopment",2000)
devindexRastDom <- function.dev.index("Domestic", "development", "nodevelopment",2000)
devindexRastInter <- function.dev.index("International", "development", "nodevelopment",2000)



#### Create maps for development and no development separately####

# development
function.develop <- function(currorigin, preference.positive, grid.size){
  pointES.pref.origin <- subset(pointES.pref, origin==currorigin)
  pos <- subset(pointES.pref.origin, preference==preference.positive)
  res(grid)<-grid.size
  posRast <- rasterize(pos, grid, fun="count")
  outRast <- posRast
  writeRaster(outRast, filename=paste0(currorigin, "_" , preference.positive),overwrite=TRUE)
  return(outRast)
}

developRastLocal <- function.develop("Local", "development",2000)
developRastDom <- function.develop("Domestic", "development",2000)
developRastInter <- function.develop("International", "development",2000)

# no development
function.nodevelop <- function(currorigin, preference.negative, grid.size){
  pointES.pref.origin <- subset(pointES.pref, origin==currorigin)
  neg <- subset(pointES.pref.origin, preference==preference.negative)
  res(grid)<-grid.size
  negRast <- rasterize(neg, grid, fun="count")
  outRast <- negRast
  writeRaster(outRast, filename=paste0(currorigin, "_" , preference.negative),overwrite=TRUE)
  return(outRast)
}

nodevelopRastLocal <- function.nodevelop("Local", "nodevelopment",2000)
nodevelopRastDom <- function.nodevelop("Domestic", "nodevelopment",2000)
nodevelopRastInter <- function.nodevelop("International", "nodevelopment",2000)

# Plots
par(mfrow=c(2,3))
plot(developRastLocal$preference,main="Local development")
plot(PAbuffer,add=TRUE)
plot(developRastDom$preference,main="Domestic development")
plot(PAbuffer,add=TRUE)
plot(developRastInter$preference,main="International development")
plot(PAbuffer,add=TRUE)
plot(nodevelopRastLocal$preference,main="Local nodevelopment")
plot(PAbuffer,add=TRUE)
plot(nodevelopRastDom$preference,main="Domestic nodevelopment")
plot(PAbuffer,add=TRUE)
plot(nodevelopRastInter$preference,main="International nodevelopment")
plot(PAbuffer,add=TRUE)

### Create maps for contrasting development preferences between different user groups (prodev vs.against dev)
devRast.conflict.LD<-developRastLocal*((-1)*nodevelopRastDom)
devRast.conflict.LI<-developRastLocal*((-1)*nodevelopRastInter)
devRast.conflict.DL<-developRastDom*((-1)*nodevelopRastLocal)
devRast.conflict.DI<-developRastDom*((-1)*nodevelopRastInter)
devRast.conflict.IL<-developRastInter*((-1)*nodevelopRastLocal)
devRast.conflict.ID<-developRastInter*((-1)*nodevelopRastDom)

par(mfrow=c(3,3))
plot(devindexRastLocal$preference,main="Local dev index")
plot(PAbuffer,add=TRUE)
plot(devRast.conflict.LD$preference,main="Local dev vs. Dom nodev")
plot(PAbuffer,add=TRUE)
plot(devRast.conflict.LI$preference,main="Local dev vs. Inter nodev")
plot(PAbuffer,add=TRUE)
plot(devRast.conflict.DL$preference,main="Dom dev vs. Local nodev")
plot(PAbuffer,add=TRUE)
plot(devindexRastDom$preference,main="Dom dev index")
plot(PAbuffer,add=TRUE)
plot(devRast.conflict.DI$preference,main="Dom dev vs. Inter nodev")
plot(PAbuffer,add=TRUE)
plot(devRast.conflict.IL$preference,main="Inter dev vs. Local nodev")
plot(PAbuffer,add=TRUE)
plot(devRast.conflict.ID$preference,main="Inter dev vs. Dom nodev")
plot(PAbuffer,add=TRUE)
plot(devindexRastInter$preference,main="Inter dev index")
plot(PAbuffer,add=TRUE)

#### Histogram for preference frequencies ####

histogram(~category2c|origin,data=subset(pointES.pref@data, !pointES.pref@data$origin=="NA"),
          index.cond=list(c(3,1,2)),layout=c(1,3))


#### Convert grid rasters to polygons to import them into ArcGIS/Qgis ####
devindexPolyLocal<-rasterToPolygons(devindexRastLocal)
devindexPolyDom<-rasterToPolygons(devindexRastDom)
devindexPolyInter<-rasterToPolygons(devindexRastInter)
devPoly.conflict.LD<-rasterToPolygons(devRast.conflict.LD)
devPoly.conflict.LI<-rasterToPolygons(devRast.conflict.LI)
devPoly.conflict.DL<-rasterToPolygons(devRast.conflict.DL)
devPoly.conflict.DI<-rasterToPolygons(devRast.conflict.DI)
devPoly.conflict.IL<-rasterToPolygons(devRast.conflict.IL)
devPoly.conflict.ID<-rasterToPolygons(devRast.conflict.ID)


writeOGR(devindexPolyLocal,"Preference maps","devindexPolyLocal",driver="ESRI Shapefile")
writeOGR(devindexPolyDom,"Preference maps","devindexPolyDom",driver="ESRI Shapefile")
writeOGR(devindexPolyInter,"Preference maps","devindexPolyInter",driver="ESRI Shapefile")
writeOGR(devPoly.conflict.LD,"Preference maps","devPoly.conflict.LD",driver="ESRI Shapefile")
writeOGR(devPoly.conflict.LI,"Preference maps","devPoly.conflict.LI",driver="ESRI Shapefile")
writeOGR(devPoly.conflict.DL,"Preference maps","devPoly.conflict.DL",driver="ESRI Shapefile")
writeOGR(devPoly.conflict.DI,"Preference maps","devPoly.conflict.DI",driver="ESRI Shapefile")
writeOGR(devPoly.conflict.IL,"Preference maps","devPoly.conflict.IL",driver="ESRI Shapefile")
writeOGR(devPoly.conflict.ID,"Preference maps","devPoly.conflict.ID",driver="ESRI Shapefile")
writeOGR(pointES.pref,"Preference maps","pointES.pref", driver="ESRI Shapefile")



