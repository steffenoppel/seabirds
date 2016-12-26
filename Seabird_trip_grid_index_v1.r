############################################################################################################
####### DEVELOP SPATIAL AGGREGATION INDEX FOR SEABIRDS   ###################################################################
############################################################################################################
## this analysis is based on Warwick-Evans et al. 2015 (tripGrid function)
## developed by steffen.oppel@rspb.org.uk in December 2016
## data provided by ana.carneiro@birdlife.org on 12 Dec 2016
## branched from Seabird_spatial_agg_index_v4.r on 26 Dec 2016


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD PACKAGES AND CUSTOM SCRIPTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(Hmisc)
require(maps)
require(mapdata)
require(adehabitatHR)
require(foreign)
require(maptools)
require(geosphere)
require(sp)
library(rgdal)
require(rgeos)
library(raster)
library(trip)
library(rworldmap)
library(plyr)
library(vegan)
library(adehabitatLT)
data(countriesLow)
library(inflection)
library(spatstat)

source("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\IBA\\Analysis\\mIBA_functions_upd2016.r")
source("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\seabird_index.r")




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

alldat<-alldat[,c(1,2,4,6:10,12:14,19,18,24:25)]
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
overview <- ddply(alldat, c("scientific_name","site_name","age","breed_stage"), summarise,n_individuals=length(unique(bird_id)), n_tracks=length(unique(ID)))
overview$DataGroup<-seq(1:dim(overview)[1])
alldat<-merge(alldat,overview[,c(1:4,7)],by=c("scientific_name","site_name","age","breed_stage"), all.x=T)


### CREATE SPATIAL SCALES FOR ANALYSIS
### order not ascending so that parallel loop will not calculate the most demanding subsets simultaneously
spatscales<-exp(seq(0,9.5,0.5))			## on log scale
spatscales<-c(2.5,500,1000,1500,5,100,150,250,7.5,10,12.5,15,17.5,20,25,30,40,50,75)			## in km
#spatscales<-c(500,30,1000,1500,40,75,100,250,5,10,15,20,50)			## in km


### Create Table that includes one line for each DataGroup at each spatial scale
### THIS WILL NEED A MANUAL ADJUSTMENT FOR THE TripSPlit function

outI<-expand.grid(overview$DataGroup,spatscales)			### sets up the lines for which we need to calculate Morisita's I
names(outI)<-c("DataGroup","Scale")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# COMBINE OVERVIEW TABLE WITH FAMILY INFO AND SET TRIP SPLIT PARAMETERS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("A://RSPB/Marine/GlobalSeabirdDistribution")
setwd("C:\\STEFFEN\\RSPB\\Marine\\GlobalSeabirdDistribution")
species_list<-read.table("SeabirdRangeMaps_SpeciesList.csv", header=T, sep=',')
species_list$GISname<-paste(substr(species_list$Family,1,3), species_list$SpcRecID, sep='')
head(species_list)
head(overview)

overview$Family<-species_list$Family[match(overview$scientific_name, species_list$Scientific_name)]
overview$Seabird_type<-species_list$Seabird_type[match(overview$scientific_name, species_list$Scientific_name)]
overview$device<-alldat$device[match(overview$scientific_name, alldat$scientific_name)]

overview$InnerBuff<-ifelse(overview$device=="PTT",10,ifelse(overview$Seabird_type=="Pelagic seabird", 5,0.5))
overview$ReturnBuff<-ifelse(overview$device=="PTT",50,ifelse(overview$Seabird_type=="Pelagic seabird", 50,5))
overview$Duration<-ifelse(overview$device=="PTT",10,ifelse(overview$Seabird_type=="Pelagic seabird", 5,1))






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# START LOOP OVER EACH DATA GROUP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

OUTPUT<-data.frame()							### set up blank frame to collate Morisita index values

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





#### IF SPARSE DATA (PTT)THEN INTERPOLATE TO EVERY 1 HR 
if(tracks$device[1]=="PTT"){
traj <- as.ltraj(xy=data.frame(tracks$Longitude, tracks$Latitude), date=as.POSIXct(tracks$TrackTime, origin="1970/01/01", tz="GMT"), id=tracks$ID, typeII = TRUE)

## Rediscretization every 3600 seconds
tr <- adehabitatLT::redisltraj(traj, 3600, type="time")

## Convert output into a data frame
tracks.intpol<-data.frame()
for (l in 1:length(unique(tracks$ID))){
out<-tr[[l]]
out$MID<-as.character(attributes(tr[[l]])[4])				#### extracts the MigID from the attribute 'id'
tracks.intpol<-rbind(tracks.intpol,out)				#### combines all data
}

### re-insert year and season

tracks.intpol$age<-tracks$age[match(tracks.intpol$MID,tracks$ID)]
tracks.intpol$bird_id<-tracks$bird_id[match(tracks.intpol$MID,tracks$ID)]
tracks.intpol$sex<-tracks$sex[match(tracks.intpol$MID,tracks$ID)]
tracks.intpol$TrackTime <- as.double(tracks.intpol$date)

## recreate data frame 'tracks' that is compatible with original data
tracks<-data.frame(scientific_name=tracks$scientific_name[1],
		site_name=tracks$site_name[1],
		age=tracks.intpol$age,
		breed_stage=tracks$breed_stage[1],
		dataset_id=tracks$dataset_id[1],
		device="PTT",
		bird_id=tracks.intpol$bird_id,
		ID=tracks.intpol$MID,
		sex=tracks.intpol$sex,
		Longitude=tracks.intpol$x,
		Latitude=tracks.intpol$y,
		DateTime=tracks.intpol$date,
		TrackTime=tracks.intpol$TrackTime,
		DataGroup=dg)
		

tracks<-tracks[order(tracks$ID, tracks$TrackTime),]
tracks$TrackTime<-adjust.duplicateTimes(tracks$TrackTime, tracks$ID)

}			## close IF loop for PTT tracks



### PROJECT COORDINATES FOR SPATIAL ANALYSES - this now needs WGS84 in CAPITAL LETTERS (not wgs84 anymore!!)
### SpatialPointsDataFrame cannot contain POSIXlt data type!

DataGroup.Wgs <- SpatialPoints(data.frame(tracks$Longitude, tracks$Latitude), proj4string=CRS("+proj=longlat +datum=WGS84"))
input <- SpatialPointsDataFrame(SpatialPoints(data.frame(tracks$Longitude, tracks$Latitude), proj4string=CRS("+proj=longlat +datum=WGS84")), data = tracks, match.ID=F)
DgProj <- CRS(paste("+proj=laea +lon_0=", loc$Longitude, " +lat_0=", loc$Latitude, sep=""))
DataGroup.Projected <- spTransform(input, CRS=DgProj)
input <- DataGroup.Projected
localmap<-spTransform(countriesLow, CRS=DgProj) 




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SPLIT INTO FORAGING TRIPS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Trips <- NULL
for(i in 1:length(unique(tracks$ID)))
  {
  	Temp <- subset(input, ID == unique(tracks$ID)[i])
	Trip <- tripSplit(Track=Temp, Colony=loc[loc$ID==unique(tracks$ID)[i],2:3], InnerBuff=overview$InnerBuff[overview$DataGroup==dg], ReturnBuff=overview$ReturnBuff[overview$DataGroup==dg], Duration = overview$Duration[overview$DataGroup==dg], plotit=T, nests=F)
  	if(i == 1) {Trips <- Trip} else
  	Trips <- spRbind(Trips,Trip)
  }




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SETTING UP A GRID AT THE APPROPRIATE RESOLUTION AND CALCULATING SPATIAL DISPERSION INDEX
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SeaPositions<-Trips[Trips@data$trip_id!="-1",]
SeaPositions<-SeaPositions[order(SeaPositions$ID, SeaPositions$TrackTime),]
nind<-length(unique(Trips@data$ID))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# START LOOP OVER EACH SPATIAL SCALE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

outMorisita<-data.frame()
for (sc in spatscales){


### CREATE A GRID BASED ON THE RESPECTIVE TripSplit output ###
grd<-makeGridTopology(SeaPositions, cellsize = c(sc*1000,sc*1000), adjust2longlat = F)			### This function assumes that cellsize is in km, but because projection is in m we need to change km into m for compatibility   

### CREATE A TRIPS OBJECT FOR THE GRID FROM Trips ###
all_trips<-trip(SeaPositions, TORnames=c("DateTime","trip_id"))			### switch to "DateTime" when using the raw locations
trg <- tripGrid(all_trips, grid=grd,method="pixellate")		### this will provide the number of bird seconds spent in each grid cell

#### CALCULATE STATISTICS FOR SPATIAL AGGREGATION ####
report<-data.frame(DataGroup=dg, Scale=sc,n_cells=length(trg@data$z),maxtime=max(trg@data$z/3600))
index<- as.numeric(dispindmorisita(trg@data$z))		### calculate index on bird seconds - if converted to hours then index gets negative because values <1 occur
report$Morisita<- index[1]
report$mclu<- index[2]
report$muni<- index[3]
report$imst<- index[4]
report$pchisq<- index[5]
outMorisita<-rbind(outMorisita,report)

} ### end loop over spatial scales

OUTPUT<-rbind(OUTPUT,outMorisita)


} ### end loop over species 



########################################







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# INSPECTING OUTPUT: EXPLORING VARIOUS SLOPES ACROSS SPECIES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### for overall results
OUTPUT$Species<-overview$scientific_name[match(OUTPUT$DataGroup,overview$DataGroup)]
#OUTPUT$DataGroup<-as.factor(OUTPUT$DataGroup)
#levels(OUTPUT$DataGroup) <- paste(c("RAZO","RAZO","FRIG","FRIG","RBTB","RBTB","SHAG","SHAG","SOAL","MUPE","MABO","MABO","MABO","MABO","CGUIL","CGUIL"),OUTPUT$DataGroup,sep="")
OUTPUT<-OUTPUT[order(OUTPUT$Scale,OUTPUT$DataGroup),] 
OUTPUT$label<-paste(c("RAZO","RAZO","FRIG","FRIG","RBTB","RBTB","SHAG","SHAG","SOAL","MUPE","MABO","MABO","MABO","MABO","CGUIL","CGUIL"),OUTPUT$DataGroup,sep="")



ggplot(OUTPUT, aes(x=log(Scale), y=Morisita), size=2)+
geom_smooth(fill="lightblue", size=1.5)+facet_wrap("label", ncol=4, scales = "free")+
  geom_point(colour="black", size=2.5) +
  xlab("log(spatial scale (km))") +
  ylab("Morisita aggregation index") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank())



ggplot(OUTPUT, aes(x=log(Scale), y=imst), size=2)+
geom_smooth(fill="lightblue", size=1.5)+facet_wrap("label", ncol=4)+
  geom_point(colour="black", size=2.5) +
  xlab("log(spatial scale (km))") +
  ylab("Standardised morisita index") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank())


####


