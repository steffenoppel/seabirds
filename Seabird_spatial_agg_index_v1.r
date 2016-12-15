############################################################################################################
####### DEVELOP SPATIAL AGGREGATION INDEX FOR SEABIRDS   ###################################################################
############################################################################################################
## this analysis is based on BirdLife International marine IBA processing scripts
## developed by steffen.oppel@rspb.org.uk in December 2016
## data provided by ana.carneiro@birdlife.org on 12 Dec 2016
## fixed some mIBA_script issues on 14 Dec 2016
## added function to remove duplicate time values on 14 Dec 2016
## initial phase to explore approaches


## based on Stats Club Discussion on 9 Dec we could proceed as follows:

1. Split trips, assess FTP and calculate home ranges for each individual
2. Use polycount output for spatial aggregation index

3. Loop over spatial scales from 100 m to 10000 km
4. For each spatial scale, draw 100 random samples of sample sizes 5,10,15,25,40
5. calculate polycount for each spatial scale
6. use values from polycount in dispindmorisita to calculate Morisita's index
7. plot Morisita's index (y-axis) over all spatial scales (x-axis) with different sample sizes (different points)
8. Fit non-linear function and determine the point at which the asymptote is reached (same function as for sample size assessment in bootstrap)

## Expectation is that the scale at which the asymptote is reached will separate species with scattered distributions from this with more aggregated distributions.







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD PACKAGES AND CUSTOM SCRIPTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(Hmisc)
require(maps)
require(mapdata)
require(adehabitat)
require(foreign)
require(maptools)
require(geosphere)
require(sp)
require(rgdal)
require(rgeos)
library(raster)
library(trip)
library(rworldmap)
library(plyr)
library(vegan)
data(countriesLow)


#source("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Statistics\\northarrow.r")
#source("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Statistics\\scalebar.r")
#source("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\IBA\\Analysis\\tripSplit.r")
#source("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\IBA\\Analysis\\tripSummary.r")
#source("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\IBA\\Analysis\\scaleARS.r")
#source("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\IBA\\Analysis\\BatchUD.r")
#source("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\IBA\\Analysis\\Bootstrap.r")
#source("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\IBA\\Analysis\\Bootstrap.r")
#source("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\IBA\\Analysis\\polyCount.r")
#source("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\IBA\\Analysis\\thresholdRaster.r")
#source("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\IBA\\Analysis\\varianceTest.r")
source("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\IBA\\Analysis\\mIBA_functions_par.r")





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD DATA FROM BIRDLIFE DATABASE AND MODIFY DATA TO MEET REQUIREMENTS FOR PROCESSING
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Load data from database
setwd("A:\\RSPB\\Marine\\SeabirdPrioritisation")
alldat<-read.table("Steffen_Spp_Prioritisation_Project_2016-12-12.csv", header=T, sep=",")
head(alldat)
dim(alldat)


### Convert Dates and Times and calculate the time difference (in s) between each location and the first location of that deployment
alldat$Time<-format(alldat$time_gmt,format="%H:%M:%S")
alldat$Date<-as.Date(alldat$date_gmt,format="%d/%m/%Y")
alldat$Loctime<-as.POSIXlt(paste(alldat$Date, alldat$Time), format = "%Y-%m-%d %H:%M:%S")
alldat$DateTime <- as.POSIXct(strptime(alldat$Loctime, "%Y-%m-%d %H:%M:%S"), "GMT")
alldat$TrackTime <- as.double(alldat$DateTime)
names(alldat)

alldat<-alldat[,c(1,3,4,6:10,12:14,19,18,24:25)]
names(alldat)[c(8,13:12)]<-c('ID','Latitude','Longitude')
head(alldat)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOT THE DATA TO CHECK THAT IT LOOKS OK
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

xlow<-min(alldat$Longitude)+0.2
xup<-max(alldat$Longitude)-0.2
yup<-max(alldat$Latitude)-0.5
ylow<-min(alldat$Latitude)+0.5

windows(600,400)
plot(Latitude~Longitude, data=alldat, pch=16, cex=0.3,col=dataset_id, asp=1, xlim=c(xlow,xup), ylim=c(ylow,yup), main="", frame=F, axes=F, xlab="", ylab="")
plot(countriesLow, col='darkgrey', add=T)






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ASSESS THE NUMBER OF DATA GROUPS (Species * Colony * life history stage) and create output table
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Create Overview table and assign data groups of equal species, site,age, and life history stage:
overview <- ddply(alldat, c("common_name","site_name","age","breed_stage"), summarise,n_individuals=length(unique(bird_id)), n_tracks=length(unique(ID)))
overview$DataGroup<-seq(1:dim(overview)[1])
alldat<-merge(alldat,overview[,c(1:4,7)],by=c("common_name","site_name","age","breed_stage"), all.x=T)


### CREATE SPATIAL SCALES FOR ANALYSIS
spatscales<-exp(seq(0,9.5,0.5))			## on log scale
spatscales<-c(0.5,1,1.5,2.5,5,7.5,10,12.5,15,17.5,20,25,30,40,50,75,100,150,250,500,1000,1500,2500,5000,10000)			## in km
spatscales<-spatscales/100			### roughly in decimal degrees, as required for the polyCount function


### Create Table that includes one line for each DataGroup at each spatial scale
### THIS WILL NEED A MANUAL ADJUSTMENT FOR THE TripSPlit function

out<-expand.grid(overview$DataGroup,spatscales)
names(out)<-c("DataGroup","Scale")
out$Morisita<-0





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# START LOOP OVER EACH DATA GROUP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



for (dg in 1:max(overview$DataGroup)){


#### SELECT DATA FOR ANALYSIS AND ORDER BY TIME AND REMOVE DUPLICATE TIME STAMPS ##########################
tracks<-alldat[alldat$DataGroup==dg,]
tracks<-tracks[order(tracks$ID, tracks$TrackTime),]
tracks$TrackTime<-adjust.duplicateTimes(tracks$TrackTime, tracks$ID)


windows(600,400)
plot(Latitude~Longitude, data=tracks, pch=16, cex=0.3,asp=1, col=tracks$ID)
plot(countriesLow, col='red', add=T)



#### CREATE COLONY LOCATION AND APPROPRIATE COORDINATE REFERENCE SYSTEM FOR PROJECTION ##########################
loc<-aggregate(lat_colony~ID, data=tracks, FUN=mean)		## Colony location is mean of all nest locations
loc$lon_colony<-aggregate(lon_colony~ID, data=tracks, FUN=mean)[,2]
names(loc)[2:3]<-c('Latitude','Longitude')


### PROJECT COORDINATES FOR SPATIAL ANALYSES
### SpatialPointsDataFrame cannot contain POSIXlt data type!

DataGroup.Wgs <- SpatialPoints(data.frame(tracks$Longitude, tracks$Latitude), proj4string=CRS("+proj=longlat + datum=wgs84"))
DgProj <- CRS(paste("+proj=laea +lon_0=", loc$Longitude, " +lat_0=", loc$Latitude, sep=""))
DataGroup.Projected <- spTransform(DataGroup.Wgs, CRS=DgProj)
input <- SpatialPointsDataFrame(DataGroup.Projected, data = tracks)
#localmap<-spTransform(countriesLow, CRS=DgProj) 



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SPLIT INTO FORAGING TRIPS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Trips <- NULL


for(i in 1:length(unique(tracks$ID)))
  {
  	Temp <- subset(input, ID == unique(tracks$ID)[i])
	Trip <- tripSplit(Track=Temp, Colony=loc[loc$ID==unique(tracks$ID)[i],2:3], InnerBuff=5, ReturnBuff=50, Duration = 5, plotit=T, nests=F)
  	if(i == 1) {Trips <- Trip} else
  	Trips <- spRbind(Trips,Trip)
  }

#str(Trips)
dim(Trips)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE SCALE FOR AREA RESTRICTED SEARCH
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DataGroup <- Trips
head(Trips@data)
ScaleOut <- scaleARS(DataGroup[DataGroup@data$trip_id!="-1",], Scales = c(seq(0.5, 50, 0.5)), Peak="Flexible")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DELINEATE CORE AREAS (KERNEL DENSITY ESTIMATOR - 50% Utilisation Distribution)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
UD<-50		## pick the % utilisation distribution (50%, 95% etc.)
Output <- batchUD(DataGroup[DataGroup@data$trip_id!="-1",], Scale = ScaleOut/2, UDLev = UD)
plot(localmap, col='darkolivegreen3', add=T)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# START LOOP OVER EACH SPATIAL SCALE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


for (sc in spatscales){


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# COUNTING THE POLYGONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ASI<-polyCount(Output, Res = sc)				### resolution in decimal degrees



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE DATA FRAME FOR MORISITA SPATIAL INDEX FROM polyCount output
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nind<-length(unique(DataGroup$ID[DataGroup@data$trip_id!="-1"]))
out$Morisita[out$DataGroup==dg & out$Scale==sc]<- dispindmorisita(ASI@data@values*nind)[1]

} ### end loop over spatial scales


} ### end loop over species 





sys.time()
system.time()













#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# USE TripGrid to calculate spatial aggregation index
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### will need interpolation?

library(trip)
#### SET THE BOUNDING BOX
long<-bbox(Trips)[1,]	### 20.360492 is westernmost nest, 44.9 is east of Yemen/Djibouti crossing
lat<-bbox(Trips)[2,]	### 11.5 is Djibouti cut-off for migration, 43.69655 is northern most EGVU nest
boundbox<-matrix(c(long,lat), ncol=2, byrow=F)


### CREATE A GRID AND COUNT NUMBER OF LOCATIONS IN EACH GRID CELL ###
grd<-makeGridTopology(boundbox, cellsize = c(10,10), adjust2longlat = F)			### CREATE THE GRID USING Trips 
extent(grd)

### CREATE A TRIPS OBJECT FOR THE GRID FROM Trips ###
all_trips<-trip(Trips, TORnames=c("DateTime","ID"))			### switch to "DateTime" when using the raw locations
trg <- tripGrid(all_trips, grid=grd,method="pixellate")						### this will provide the number of bird seconds spent in each grid cell
spplot(trg)			## plots the trips with a legend
proj4string(turbSPDF)<-proj4string(EVSP_all)


### CONVERT SPATIAL GRID TO SOMETHING WE CAN PLOT
spdf <- SpatialPixelsDataFrame(points=trg, data=trg@data)
HOTSPOTS<-data.frame(lat=spdf@coords[,2],long=spdf@coords[,1],time=spdf@data$z)
HOTSPOTS$time<-HOTSPOTS$time/(3600*24)								### this converts the seconds into bird days
summary(HOTSPOTS)

HOTSPOTS<-HOTSPOTS[HOTSPOTS$time>10,]

##### CONVERT TO SPATIAL POLYGONS FOR OVERLAY ###
ras<- raster(spdf)		# converts the SpatialPixelDataFrame into a raster
spoldf <- rasterToPolygons(ras, n=4) # converts the raster into quadratic polygons
proj4string(spoldf)<-proj4string(EVSP_all)


### PRODUCE NICE AND SHINY MAP WITH THE MOST IMPORTANT TEMPORARY CONGREGATION SITES ###

#MAP <- get_map(EGVUbox, source="google", zoom=4, color = "bw")		### retrieves a map from Google (requires live internet connection)
MAP <- get_map(location = c(lon = mean(long), lat = mean(lat)), source="google", zoom=4, color = "bw")		### retrieves a map from Google (requires live internet connection)

pdf("EGVU_MIGRATION_HOTSPOTS.pdf", width=12, height=11)
ggmap(MAP)+geom_tile(data=HOTSPOTS, aes(x=long,y=lat, fill = time)) +
	scale_fill_gradient(name = 'N bird days', low="white", high="red", na.value = 'transparent', guide = "colourbar", limits=c(10, 30))+
	theme(axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank())+
	theme(strip.text.y = element_text(size = 20, colour = "black"), strip.text.x = element_text(size = 15, colour = "black"))+
	geom_point(data=turbs, aes(x=Longitude, y=Latitude), pch=16, col='darkolivegreen', size=0.5)
dev.off()



