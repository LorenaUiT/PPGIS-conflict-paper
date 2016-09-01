### Libraries ####
library(reshape2)
library(raster)
library(sp)
library(rgdal)
library(rgeos)
library(grid)
library(lattice) 
library(mvabund)  ###glm
library (MASS)  ###chi-sq
library(spatstat)
library(maptools) ##to create window to use in spatstat
library(graphics) ## pairs function
library(dbscan)   ## to calculate distance between clusters and find out how many clusters there are
library(sm) ## to compare point densities of different user groups

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
pointES.values@data$category2c<-factor(pointES.values@data$category2c)
pointES.pref<-subset(pointES,!is.na(pointES@data$preference))
pointES.pref@data$category2c<-factor(pointES.pref@data$category2c)

#### Chi-square ####
#### Check whether there are differences in place values and preferences between user groups
## 
# values
pointES.values.table<-table(pointES.values@data$origin,pointES.values@data$category2c)
pointES.values.table<-subset(pointES.values.table[-4,]) ## Deleting the last row for origin=NA
chisq.test(pointES.values.table) ### p<0.05 --> user groups differ in place values
chisq.test(pointES.values.table)$stdres

# preferences
pointES.pref.table<-table(pointES.pref@data$origin,pointES.pref@data$category2c)
pointES.pref.table<-subset(pointES.pref.table[-4,]) ## Deleting the last row for origin=NA
chisq.test(pointES.pref.table) ###there are >20% of <5 values, chi-sq is not optimal
chisq.test(pointES.pref.table)$stdres
kruskal.test(origin~category2c,data=pointES.pref@data) ### p<0.05 --> user groups differ in preferences


#### COMPARE DENSITIES OF USER GROUPS  #####
### To see whether they use the same areas or not.




##### PLACE VALUES #####
####
#### Create polygons with points that are nearby for each user group 
window<-as.owin(PAbufferdis) ## we need to create an observation window to delimit the study area
m<-pointES.values@data[,c(6,11)] # We put origin as marks, so we are able to identify areas of different user groups
m<-droplevels(m)
## Warning: there are duplicated points in pointES.values that appear if marks="m" is not included in the model(it's a common thing)

pointES.values.ppp<-ppp(x=pointES.values@coords[,1],y=pointES.values@coords[,2], window = window,marks = m)
is.multitype(pointES.values.ppp) ##Says true when there is a single mark
pointES.values.ppp<-subset(pointES.values.ppp[pointES.values.ppp$marks$origin!="NA"])
pointES.values.ppp$n/sum(sapply(slot(PAbufferdis, "polygons"), slot, "area")) ### Calculates average intensity of markers
plot(pointES.values.ppp)
plot(Kest(pointES.values.ppp)) ## How do I interpret this???

### STEP 1: Overal distribution of values
plot(density(pointES.values.ppp,1000))
persp(density(pointES.values.ppp,1000))
plot(density(pointES.values.ppp[pointES.values.ppp$marks$origin=="Local"],1000),main = "Local's place value distribution")
plot(density(pointES.values.ppp[pointES.values.ppp$marks$origin=="Domestic"],1000),main = "Domestic's place value distribution")
plot(density(pointES.values.ppp[pointES.values.ppp$marks$origin=="International"],1000),main = "International's place value distribution")

## Kolmogorov-Smirnov test to compare expected and observed distributions of x and y
# kstest is out of date, R recomends using cdf.test instead
cdf.test(pointES.values.ppp, function(x, y) {x}) ## Rejects H0, therefore expected and observed distributions of x are different
plot(cdf.test(pointES.values.ppp, function(x, y) {x}))
cdf.test(pointES.values.ppp, function(x, y) {y}) ## Rejects H0, therefore expected and observed distributions of y are different
plot(cdf.test(pointES.values.ppp, function(x, y) {y}))

###STEP 2: Segregation of values
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

## We cannot use MLE for poisson because our points are not completely independent (is this right?)
## MLE tests covariates?

#### COMPARE DENSITIES OF USER GROUPS  #####
### To see whether they use the same areas or not.
nn.sim = vector()
P.r = pointES.values.ppp
marks(P.r) = pointES.values.ppp$marks$origin
for(i in 1:999){
  marks(P.r) = sample(pointES.values.ppp$marks$origin,replace = FALSE)  # Reassign labels at random, point locations don't change
  nn.sim[i] = mean(nncross(subset(P.r[P.r$marks=="Local"]),subset(P.r[P.r$marks=="Domestic"]))$dist)
}

hist(nn.sim,breaks=30)
abline(v=mean(nncross(pointES.values.ppp[pointES.values.ppp$marks$origin=="Local"],pointES.values.ppp[pointES.values.ppp$marks$origin=="Domestic"])$dist),col="red")
mean(nncross(pointES.values.ppp[pointES.values.ppp$marks$origin=="Local"],pointES.values.ppp[pointES.values.ppp$marks$origin=="Domestic"])$dist)
# Compute empirical cumulative distribution
nn.sim.ecdf = ecdf(nn.sim)
plot(nn.sim.ecdf)
# See how the original stat compares to the simulated distribution
nn.sim.ecdf(mean(nncross(subset(P.r[P.r$marks=="Local"]),subset(P.r[P.r$marks=="Domestic"]))$dist)) 
summary(nn.sim.ecdf)







nn.sim1 = vector()
P.r = pointES.values.ppp
marks(P.r) = pointES.values.ppp$marks$origin
for(i in 1:999){
  marks(P.r) = sample(pointES.values.ppp$marks$origin,replace = FALSE)  # Reassign labels at random, point locations don't change
  nn.sim1[i] = mean(nncross(subset(P.r[P.r$marks=="Local"]),subset(P.r[P.r$marks=="International"]))$dist)
}

hist(nn.sim1,breaks=30)
abline(v=mean(nncross(pointES.values.ppp[pointES.values.ppp$marks$origin=="Local"],pointES.values.ppp[pointES.values.ppp$marks$origin=="International"])$dist),col="red")
mean(nncross(pointES.values.ppp[pointES.values.ppp$marks$origin=="Local"],pointES.values.ppp[pointES.values.ppp$marks$origin=="International"])$dist)
# Compute empirical cumulative distribution
nn.sim1.ecdf = ecdf(nn.sim1)
plot(nn.sim1.ecdf)
# See how the original stat compares to the simulated distribution
nn.sim1.ecdf(mean(nncross(subset(P.r[P.r$marks=="Local"]),subset(P.r[P.r$marks=="International"]))$dist)) 
summary(nn.sim1.ecdf)






nn.sim2 = vector()
P.r = pointES.values.ppp
marks(P.r) = pointES.values.ppp$marks$origin
for(i in 1:999){
  marks(P.r) = sample(pointES.values.ppp$marks$origin,replace = FALSE)  # Reassign labels at random, point locations don't change
  nn.sim2[i] = mean(nncross(subset(P.r[P.r$marks=="Domestic"]),subset(P.r[P.r$marks=="International"]))$dist)
}

hist(nn.sim2,breaks=30)
abline(v=mean(nncross(pointES.values.ppp[pointES.values.ppp$marks$origin=="Domestic"],pointES.values.ppp[pointES.values.ppp$marks$origin=="International"])$dist),col="red")
mean(nncross(pointES.values.ppp[pointES.values.ppp$marks$origin=="Domestic"],pointES.values.ppp[pointES.values.ppp$marks$origin=="International"])$dist)
# Compute empirical cumulative distribution
nn.sim2.ecdf = ecdf(nn.sim2)
plot(nn.sim2.ecdf)
# See how the original stat compares to the simulated distribution
nn.sim2.ecdf(mean(nncross(subset(P.r[P.r$marks=="Domestic"]),subset(P.r[P.r$marks=="International"]))$dist)) 
summary(nn.sim2.ecdf)



###STEP 3: Identify overlaping hotspots

## Calculate distance between clusters so we can use it in the dbscan function (eps)
## K=2*dimension-1
## to determine eps: http://www.sthda.com/english/wiki/dbscan-density-based-clustering-for-discovering-clusters-in-large-datasets-with-noise-unsupervised-machine-learning
kNNdist(pointES.values@coords,k=3)
kNNdistplot(pointES.values@coords,k=3) ### eps is around 1100

pointES.values.dbscan<-dbscan(pointES.values@coords,eps = 1100,minPts = 10)
plot(pointES.values, col=pointES.values.dbscan$cluster,main="Value clusters")
text(polygon.values,labels=getSpPPolygonsIDSlots(polygon.values),cex=4)
plot(PA,add=TRUE)



# add a column called cluster to the data frame of each user group values
pointES.values@data$cluster<-pointES.values.dbscan$cluster
pointES.values<-subset(pointES.values,pointES.values@data$origin!="NA")
pointES.values@data$origin<-droplevels(pointES.values@data$origin)



## Create polygons with each cluster identified in dbscan

polygon.values<-SpatialPolygons(list(),proj4string = pointES.values@proj4string)

for (  i in unique(pointES.values@data[,18])[unique(pointES.values@data[,18])!=0]) {
  polygon.values.hull<-gConvexHull(pointES.values[pointES.values@data$cluster==i,])
  polygon.values.hull@polygons[[1]]@ID=as.character(i)
  polygon.values<-spRbind(polygon.values,polygon.values.hull)

  }


plot(polygon.values,main="Value clusters",col=1:14)
text(polygon.values,labels=getSpPPolygonsIDSlots(polygon.values),cex=4)
plot(PA,add=TRUE)


## Compare categories of each user group in each polygon

for ( i in unique(pointES.values@data$cluster)){
  pointES.values.clusters<-table(pointES.values@data$origin[pointES.values@data$cluster==i],pointES.values@data$category2c[pointES.values@data$cluster==i])
  print(paste("Overlapping with local polygons in cluster no", i))
  if(!is.na(kruskal.test(origin~category2c,data=subset(pointES.values@data,pointES.values@data$cluster==i))$p.value)){
    if(kruskal.test(origin~category2c,data=subset(pointES.values@data,pointES.values@data$cluster==i))$p.value<0.05){
      print(pointES.values.clusters)
      print(chisq.test(pointES.values.clusters)$stdres)
      print(kruskal.test(origin~category2c,data=subset(pointES.values@data,pointES.values@data$cluster==i)))
    }
  }  else
    print(paste("Cluster",i,"not significantly different"))
}









##### DEVELOPMENT PREFERENCES #####
####
#### Create polygons with points that are nearby for each user group 
window<-as.owin(PAbufferdis) ## we need to create an observation window to delimit the study area
m1<-pointES.pref@data[,c(6,7,11)] # We put origin as marks, so we are able to identify areas of different user groups
m1<-droplevels(m1)
## Warning: there are duplicated points in pointES.pref that appear if marks="m" is not included in the model(it's a common thing)

pointES.pref.ppp<-ppp(x=pointES.pref@coords[,1],y=pointES.pref@coords[,2], window = window,marks = m1)
is.multitype(pointES.pref.ppp) ##Says true when there is a single mark
pointES.pref.ppp$n/sum(sapply(slot(PAbufferdis, "polygons"), slot, "area")) ### Calculates average intensity of markers
plot(pointES.pref.ppp)
plot(Kest(pointES.pref.ppp)) ## How do I interpret this???


### STEP1: Overal distribution of preferences
plot(density(pointES.pref.ppp,1000))
persp(density(pointES.pref.ppp,1000))
plot(density(pointES.pref.ppp[pointES.pref.ppp$marks$origin=="Local"],1000),main = "Local's development preference distribution")
plot(density(pointES.pref.ppp[pointES.pref.ppp$marks$origin=="Domestic"],1000),main = "Domestic's development preference distribution")
plot(density(pointES.pref.ppp[pointES.pref.ppp$marks$origin=="International"],1000),main = "International's development preference distribution")

## Kolmogorov-Smirnov test to compare expected and observed distributions of x and y
# kstest is out of date, R recomends using cdf.test instead
cdf.test(pointES.pref.ppp, function(x, y) {x}) ## Rejects H0, therefore expected and observed distributions of x are different
plot(cdf.test(pointES.pref.ppp, function(x, y) {x}))
cdf.test(pointES.pref.ppp, function(x, y) {y}) ## Rejects H0, therefore expected and observed distributions of y are different
plot(cdf.test(pointES.pref.ppp, function(x, y) {y}))


###STEP 2: Segregation of preferences
### Nearest neighbour distances
ave(Gest(pointES.pref.ppp[pointES.pref.ppp$marks$origin=="Local"])$r) #Result: 3645.003
plot(Gest(pointES.pref.ppp[pointES.pref.ppp$marks$origin=="Local"])) #Cumulative distribution function of the nearest-neighbour distance
ave(Gest(pointES.pref.ppp[pointES.pref.ppp$marks$origin=="Domestic"])$r) #Result: 3584.249
plot(Gest(pointES.pref.ppp[pointES.pref.ppp$marks$origin=="Domestic"]))
ave(Gest(pointES.pref.ppp[pointES.pref.ppp$marks$origin=="International"])$r) #Result: 3308.618
plot(Gest(pointES.pref.ppp[pointES.pref.ppp$marks$origin=="International"]))
# Development preferences by international tourists are the ones with the shorterst avegare distance between neighbouring points
# This may mean that international tourists have a lower dispersion than Locals (highest average distance between points)
# Plots show that bG(r) > Gpois(r) suggest a clustered pattern.
# This can be concluded from the plots as they show a rapid increase on the accummulated distance, which means that points are very close together


#### COMPARE DENSITIES OF USER GROUPS  #####
### To see whether they use the same areas or not.
nn.sim3 = vector()
P.rp = pointES.pref.ppp
marks(P.rp) = pointES.pref.ppp$marks$origin
for(i in 1:999){
  marks(P.rp) = sample(pointES.pref.ppp$marks$origin,replace = FALSE)  # Reassign labels at random, point locations don't change
  nn.sim3[i] = mean(nncross(subset(P.rp[P.rp$marks=="Local"]),subset(P.rp[P.rp$marks=="Domestic"]))$dist)
}

hist(nn.sim3,breaks=30)
abline(v=mean(nncross(pointES.pref.ppp[pointES.pref.ppp$marks$origin=="Local"],pointES.pref.ppp[pointES.pref.ppp$marks$origin=="Domestic"])$dist),col="red")
mean(nncross(pointES.pref.ppp[pointES.pref.ppp$marks$origin=="Local"],pointES.pref.ppp[pointES.pref.ppp$marks$origin=="Domestic"])$dist)
# Compute empirical cumulative distribution
nn.sim3.ecdf = ecdf(nn.sim)
plot(nn.sim3.ecdf)
# See how the original stat compares to the simulated distribution
nn.sim3.ecdf(mean(nncross(subset(P.rp[P.rp$marks=="Local"]),subset(P.rp[P.rp$marks=="Domestic"]))$dist)) 
summary(nn.sim3.ecdf)







nn.sim4 = vector()
P.rp = pointES.pref.ppp
marks(P.rp) = pointES.pref.ppp$marks$origin
for(i in 1:999){
  marks(P.rp) = sample(pointES.pref.ppp$marks$origin,replace = FALSE)  # Reassign labels at random, point locations don't change
  nn.sim4[i] = mean(nncross(subset(P.rp[P.rp$marks=="Local"]),subset(P.rp[P.rp$marks=="International"]))$dist)
}

hist(nn.sim4,breaks=30)
abline(v=mean(nncross(pointES.pref.ppp[pointES.pref.ppp$marks$origin=="Local"],pointES.pref.ppp[pointES.pref.ppp$marks$origin=="International"])$dist),col="red")
mean(nncross(pointES.pref.ppp[pointES.pref.ppp$marks$origin=="Local"],pointES.pref.ppp[pointES.pref.ppp$marks$origin=="International"])$dist)
# Compute empirical cumulative distribution
nn.sim4.ecdf = ecdf(nn.sim1)
plot(nn.sim4.ecdf)
# See how the original stat compares to the simulated distribution
nn.sim4.ecdf(mean(nncross(subset(P.rp[P.rp$marks=="Local"]),subset(P.rp[P.rp$marks=="International"]))$dist)) 
summary(nn.sim4.ecdf)






nn.sim5 = vector()
P.rp = pointES.pref.ppp
marks(P.rp) = pointES.pref.ppp$marks$origin
for(i in 1:999){
  marks(P.rp) = sample(pointES.pref.ppp$marks$origin,replace = FALSE)  # Reassign labels at random, point locations don't change
  nn.sim5[i] = mean(nncross(subset(P.rp[P.rp$marks=="Domestic"]),subset(P.rp[P.rp$marks=="International"]))$dist)
}

hist(nn.sim5,breaks=30)
abline(v=mean(nncross(pointES.pref.ppp[pointES.pref.ppp$marks$origin=="Domestic"],pointES.pref.ppp[pointES.pref.ppp$marks$origin=="International"])$dist),col="red")
mean(nncross(pointES.pref.ppp[pointES.pref.ppp$marks$origin=="Domestic"],pointES.pref.ppp[pointES.pref.ppp$marks$origin=="International"])$dist)
# Compute empirical cumulative distribution
nn.sim5.ecdf = ecdf(nn.sim2)
plot(nn.sim5.ecdf)
# See how the original stat compares to the simulated distribution
nn.sim5.ecdf(mean(nncross(subset(P.rp[P.rp$marks=="Domestic"]),subset(P.rp[P.rp$marks=="International"]))$dist)) 
summary(nn.sim5.ecdf)







###STEP 3: Identify overlaping hotspots
## For this I will separate points in three layers (user types)
pointES.pref.local<-subset(pointES.pref,pointES.pref@data$origin=="Local")
pointES.pref.dom<-subset(pointES.pref,pointES.pref@data$origin=="Domestic")
pointES.pref.inte<-subset(pointES.pref,pointES.pref@data$origin=="International")


plot(PA,main="Preferences of Locals")
plot(pointES.pref.local,add=TRUE,pch=2)
plot(PA,main="Preferences of Domestic")
plot(pointES.pref.dom,add=TRUE,pch=1)
plot(PA,main="Preferences of Internationals")
plot(pointES.pref.inte,add=TRUE,pch=3)


## Calculate distance between clusters so we can use it in the dbscan function (eps)
## K=2*dimension-1
kNNdist(pointES.pref@coords,k=3)
kNNdistplot(pointES.pref@coords,k=3)

pointES.pref.dbscan<-dbscan(pointES.pref@coords,eps = 2500,minPts = 10)
plot(PA,main="Preference clusters")
plot(pointES.pref, col=pointES.pref.dbscan$cluster,add=TRUE)


# add a column called cluster to the data frame of each user group preferences
pointES.pref@data$cluster<-pointES.pref.dbscan$cluster


## Create polygons with each cluster identified in dbscan

polygon.pref<-SpatialPolygons(list(),proj4string = pointES.pref@proj4string)

for (  i in unique(pointES.pref@data[,18])[unique(pointES.pref@data[,18])!=0]) {
  polygon.pref.hull<-gConvexHull(pointES.pref[pointES.pref@data$cluster==i,])
  polygon.pref.hull@polygons[[1]]@ID=as.character(i)
  polygon.pref<-spRbind(polygon.pref,polygon.pref.hull)
  
}


plot(PA,main="Preference clusters")
plot(polygon.pref,col=getSpPPolygonsIDSlots(polygon.pref),add=TRUE)
text(polygon.pref,labels=getSpPPolygonsIDSlots(polygon.pref),cex=4)

pointES.pref<-subset(pointES.pref,pointES.pref@data$origin!="NA")
pointES.pref@data$origin<-droplevels(pointES.pref@data$origin)



## Compare categories of each user group in each polygon

for ( i in unique(pointES.pref@data$cluster)){
  pointES.pref.clusters<-table(pointES.pref@data$origin[pointES.pref@data$cluster==i],pointES.pref@data$preference[pointES.pref@data$cluster==i])
  print(paste("Overlapping with cluster no", i))
  if(!is.na(kruskal.test(origin~preference,data=subset(pointES.pref@data,pointES.pref@data$cluster==i))$p.value)){
    if(kruskal.test(origin~preference,data=subset(pointES.pref@data,pointES.pref@data$cluster==i))$p.value<0.05){
      print(pointES.pref.clusters)
      print(chisq.test(pointES.pref.clusters)$stdres)
      print(kruskal.test(origin~preference,data=subset(pointES.pref@data,pointES.pref@data$cluster==i)))
    }
  }  else
    print(paste("Cluster",i,"not significantly different"))
}






