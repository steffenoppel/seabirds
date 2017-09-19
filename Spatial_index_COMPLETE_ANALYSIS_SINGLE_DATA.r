############################################################################################################
####### SPATIAL AGGREGATION INDEX ANALYSIS FOR SEABIRDS   ##################################################
############################################################################################################
## this analysis is based on BirdLife International marine IBA processing scripts
## developed by steffen.oppel@rspb.org.uk in December 2016
## data provided by ana.carneiro@birdlife.org on 12 Dec 2016 and in 2017

## RE-ANALYSIS FOR ALL SPECIES
## see 'Seabird_spatial_agg_index_v8 - 9' for previous analysis descriptions and notes
## some datasets in the main analysis were found to include erroneous data
## 22 Aug 2017 - complete cleaning of all data with speed filter in 'Seabird_spatial_agg_data_preparation_v3.R'
## 23 Aug 2017 - build entire analysis in one loop, includes raw data analysis, and simulation across sample sizes

## 24 Aug 2017 - split from COMPLETE REANALYSIS script to provide standalone script for single straggler datasets - including data preparation

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD PACKAGES AND CUSTOM SCRIPTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(Hmisc)
library(data.table)
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
library(vegan)
library(adehabitatLT)
library(inflection)
library(spatstat)
library(move)
library(R.utils)
require(foreach)
require(doParallel)
require(parallel)
library(lubridate)
library(tidyverse)
library(dplyr)
library(ggplot2)

#install.packages("devtools")
#require(devtools)
#install_github("TakahiroShimada/SDLfilter")
library(SDLfilter)


source("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\IBA\\Analysis\\mIBA_functions_par.r")
source("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Analysis\\seabird_index3.r")

#source("C:\\STEFFEN\\RSPB\\Marine\\IBA\\Analysis\\mIBA_functions_par.r")
#source("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation\\seabird_index3.r")




#################################################################### DATA PREPARATION ##############################################################################




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD DATA FROM BIRDLIFE DATABASE AND MODIFY DATA TO MEET REQUIREMENTS FOR PROCESSING
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Load data from database
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation")
#setwd("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation")
alldat<-fread("Data\\Falklands_raw_data.csv", header=T, sep=",")

alldat <- alldat %>%
	mutate(DateTime=ymd_hms(paste(date_gmt, time_gmt))) %>%
	mutate(TrackTime=as.double(DateTime)) %>%
	dplyr::select(c(1,2,4:10,12:14,19,18,20),DateTime,TrackTime)				### position for last two columns doesn't work because of inconsistency among data

names(alldat)
names(alldat)[c(8,9,14:13)]<-c('bird_id','ID','Latitude','Longitude')
head(alldat)
dim(alldat)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# REMOVE DATA FROM OUTSIDE THE BREEDING SEASON 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
removal<-c("migration","non-breeding","pre-egg","pre-moult")
alldat<-alldat %>%
  filter(!(breed_stage %in% removal)) %>%
  mutate(breed_stage=as.character(breed_stage)) %>%
  mutate(breed_stage=ifelse(breed_stage=="chick-rearing","brood-guard",breed_stage)) %>%
  mutate(breed_stage=ifelse(breed_stage=="unknown","breeding",breed_stage)) %>%
  mutate(breed_stage=ifelse(breed_stage=="creche","post-guard",breed_stage))

dim(alldat)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# COMBINE OVERVIEW TABLE WITH FAMILY INFO AND SET TRIP SPLIT PARAMETERS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Data")
#setwd("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation\\Data")
overview<-fread("Seabird_priority_overview_v4.csv")
DG<-max(overview$DataGroup)+1
head(overview)


setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation")
#setwd("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation")

species_list<-read.table("Seabird_SpeciesList.csv", header=T, sep=',')
species_list$GISname<-paste(substr(species_list$Family,1,3), species_list$SpcRecID, sep='')
head(species_list)


### Create Overview table and assign data groups of equal species, site,age, and life history stage:
overview2 <- alldat %>%
	group_by(scientific_name,site_name,colony_name,age,breed_stage,device) %>%
	summarise(n_individuals=length(unique(bird_id)), n_tracks=length(unique(ID))) %>%
	mutate(Family=species_list$Family[match(scientific_name, species_list$Scientific_name)]) %>%
	mutate(Seabird_type=species_list$Seabird_type[match(scientific_name, species_list$Scientific_name)]) %>%
	mutate(InnerBuff=ifelse(device=="PTT",10,ifelse(Seabird_type=="Pelagic seabird", 5,0.5))) %>%
	mutate(ReturnBuff=ifelse(device=="PTT",50,ifelse(Seabird_type=="Pelagic seabird", 50,5))) %>%
	mutate(Duration=ifelse(device=="PTT",10,ifelse(Seabird_type=="Pelagic seabird", 5,1)))

overview2$DataGroup<-seq(DG,(DG+dim(overview2)[1]-1),1)
overview2<-as.data.frame(overview2)

names(overview2)

overview<-rbind(overview[,c(1:14)],overview2[,c(1:5,7,8,14,9,10,6,11:13)])

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Data")
#setwd("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation\\Data")

fwrite(overview,"Seabird_priority_overview_v4.csv")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# COMBINE DATA WITH OVERVIEW TABLE AND REMOVE LOCATIONS 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
names(overview)


alldat<-merge(alldat,overview[,c(1:5,8)],by=c("scientific_name","site_name","colony_name","age","breed_stage"), all.x=T)
dim(alldat)
head(alldat)

ALL<- alldat %>%
	filter(!is.na(Longitude)) %>%
	filter(!is.na(Latitude)) %>%
	filter(!is.na(DataGroup))	

dim(ALL)
names(ALL)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOT THE DATA TO CHECK THAT IT LOOKS OK
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(rworldmap)
data(countriesLow)
xlow<-min(ALL$Longitude)+0.2
xup<-max(ALL$Longitude)-0.2
yup<-max(ALL$Latitude)-0.5
ylow<-min(ALL$Latitude)+0.5

windows(600,400)
plot(Latitude~Longitude, data=ALL, pch=16, cex=0.3,col=dataset_id, asp=1, xlim=c(xlow,xup), ylim=c(ylow,yup), main="", frame=F, axes=F, xlab="", ylab="")
plot(countriesLow, add=T)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXPORT INPUT DATA FOR EACH NEW DATAGROUP INTO SEPARATE CSV FILE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm('alldat')
gc(verbose = TRUE)


dgs<-overview2$DataGroup[overview2$n_tracks>4]
pttqi<-data.frame(argos_lc=c("3","2","1","0","A","B","Z"), qi=c(10,8,6,3,2,1,0))


#setup parallel backend to use 4 processors
cl<-makeCluster(4)
registerDoParallel(cl)

foreach(dg=dgs, .packages=c("tidyverse","SDLfilter","adehabitatHR","geosphere","rgdal","sp","adehabitatLT","trip","data.table")) %dopar% {


#dg=201

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SELECT DATA FOR ANALYSIS AND MANUALLY CLEAN THE DATASET BY FILTERING OUT SPECIFIED LOCATIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Data\\DataGroups")


### SIMPLE SPEED FILTER
## based on: Shimada T, Jones R, Limpus C, Hamann M (2012) Improving data retention and home range estimates by data-driven screening. Mar Ecol Prog Ser 457:171-180 http://dx.doi.org/10.3354/meps09747

tracks<-ALL %>%
	filter(DataGroup==dg) %>%
	filter(argos_quality %in% c("3","2","1","0","A","B",NA,"")) %>%
	arrange(ID, TrackTime)

sdata <- tracks
names(sdata)[c(10,11,14:13,15)]<-c('bird_id','id','lat','lon','qi')
sdata$qi<-ifelse(is.na(sdata$qi),3,sdata$qi)
sdata$qi<-ifelse(sdata$device=="PTT",pttqi$qi[match(sdata$qi,pttqi$argos_lc)],sdata$qi)
tripout<-ddfilter.speed(sdata, vmax = ifelse(overview$Family[overview$DataGroup==dg]=="Spheniscidae",15,120), method = 2)
names(tripout)[c(10,11,14:13,15)]<-names(ALL)[c(10,11,14:13,15)]
dim(tripout)


### WRITE GOOGLE EARTH FILE FOR EASIER INSPECTION OF OUTLIERS
#input <- SpatialPointsDataFrame(SpatialPoints(data.frame(tripout$Longitude, tripout$Latitude), proj4string=CRS("+proj=longlat +datum=WGS84")), data = tripout, match.ID=F)
#writeOGR(input, sprintf("SeabirdTracks_DataGroup%s.kml",dg), layer=sprintf("SeabirdTracks_DataGroup%s",dg), driver="KML", overwrite_layer=T)



### SEE WHETHER IT HAS WORKED

#windows(600,400)
#plot(Latitude~Longitude, data=tripout, pch=1, cex=0.8, asp=1, main="", frame=F, axes=F, xlab="", ylab="")
#points(Latitude~Longitude, data=tracks, pch=4, cex=0.3,col='red')
#plot(countriesLow, col='darkgrey', add=T)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FURTHER PROCESSING OF DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tracks<-tripout[!(is.infinite(tripout$sSpeed)),]						### this removes rows with an 'infinite' speed due to time diff = 0
tracks$TrackTime<-adjust.duplicateTimes(tracks$TrackTime, tracks$ID)


### troubleshoot the dreaded non-unique error
#duplicates <- tracks %>%
#	mutate(count=1) %>%
#	group_by(ID, TrackTime) %>%
#	summarise(nlocs=sum(count)) %>%
#	filter(nlocs>1)
#duplicates<-merge(duplicates, tracks, by=c("ID","TrackTime"),all.x=T)
 


#### REMOVE MYSTERIOUSLY CREATED NA LINES
dim(tracks)
tracks<-tracks[!is.na(tracks$Longitude),]
tracks<-tracks[!is.na(tracks$Latitude),]
dim(tracks)


#### IF SPARSE DATA (PTT)THEN INTERPOLATE TO EVERY 1 HR 
if(tracks$device[1]=="PTT"){
if(NA %in% tracks$argos_quality){}else{tracks<-tracks[!(tracks$argos_quality=="Z"),]}
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
		colony_name=tracks$colony_name[1],
		age=tracks.intpol$age,
		breed_stage=tracks$breed_stage[1],
		dataset_id=tracks$dataset_id[1],
		lat_colony=tracks$lat_colony[1],
		lon_colony=tracks$lon_colony[1],
		device="PTT",
		bird_id=tracks.intpol$bird_id,
		ID=tracks.intpol$MID,
		sex=tracks.intpol$sex,
		Longitude=tracks.intpol$x,
		Latitude=tracks.intpol$y,
		argos_quality="2",
		DateTime=tracks.intpol$date,
		TrackTime=tracks.intpol$TrackTime,
		DataGroup=dg)
		

tracks<-tracks[order(tracks$ID, tracks$TrackTime),]
tracks$TrackTime<-adjust.duplicateTimes(tracks$TrackTime, tracks$ID)

}			## close IF loop for PTT tracks


############ WRITE TRACKS TO CSV TABLE THAT CAN BE READ IN FOR FUTURE ANALYSIS ###########
fwrite(tracks,sprintf("CLEAN_SeabirdTracks_DataGroup%s.csv",dg))



}			## close loop over all data groups


## stop the cluster
stopCluster(cl)
print(Sys.time())




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# OVERVIEW TABLE FOR ALL DATAGROUPS AND SPECIES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

overview<- overview %>%
  arrange(DataGroup) %>%
  select(DataGroup, Family, scientific_name,site_name,colony_name,breed_stage,n_individuals,n_tracks,InnerBuff,ReturnBuff,Duration)

head(overview)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SET SPATIAL AND ITERATION PARAMETERS FOR BOOTSTRAPPING ANALYSIS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## number of iterations and sample sizes that are randomly drawn
### CREATE SERIES from 3 to max n trips THAT IS NOT COMPUTATIONALLY EXCESSIVE

ITERA=25
nloops<-c(3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,23,26,29,32,35,40,45,50,55,60,65,70,75,85,100,125,150,175,200,250,300,310)## max number is 301 from previous analysis


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD PREVIOUSLY SAVED RESULTS AND PREPARE OUTPUT DATA FRAMES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs")

rm('ALL','tripout','bootstrap','cl','DG','dg','tracks','species_list','sdata')
gc()


### CHECK WHAT DATA ALREADY EXIST
#donefiles<-list.files(path="S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs", pattern=".*SeabirdOverlapIndex_DataGroup.*\\.csv")

					
#dgs<-dgs[-1]
#dgs<-dgs[order(dgs, decreasing=T)]


TRIP_SUMMARY<-data.frame()
SpatInd_SUMMARY<-data.frame()






#################################################################### ANALYSIS ##############################################################################






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD OVERVIEW TABLE FOR ALL DATAGROUPS AND SPECIES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Data")
#setwd("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation\\Data")

overview<-fread("Seabird_priority_overview_v4.csv", header=T, sep=',')


overview<- overview %>%
  mutate(ReturnBuff=as.numeric(ReturnBuff)) %>%
  mutate(Duration=as.numeric(Duration)) %>%
  arrange(DataGroup) %>%
  select(DataGroup, Family, scientific_name,site_name,colony_name,breed_stage,n_individuals,n_tracks,InnerBuff,ReturnBuff,Duration)

head(overview)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# START OF SPECIES-SPECIFIC INDEX CALCULATION FOR EACH DATA GROUP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (dg in dgs){        ### START OVERALL LOOP OVER EACH DATA GROUP (deliberately kept serial to facilitate troubleshooting)

#### SELECT DATA FOR ANALYSIS ##########################
tracksname<-sprintf("CLEAN_SeabirdTracks_DataGroup%s.csv",dg)
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Data\\DataGroups")
tracks<-read.table(tracksname, header=T, sep=",")
tracks$DateTime <- ymd_hms(tracks$DateTime, tz="GMT")
  

#### CREATE COLONY LOCATION AND APPROPRIATE COORDINATE REFERENCE SYSTEM FOR PROJECTION ##########################
loc<-aggregate(lat_colony~ID, data=tracks, FUN=mean)		## Colony location is mean of all nest locations
loc$lon_colony<-aggregate(lon_colony~ID, data=tracks, FUN=mean)[,2]
names(loc)[2:3]<-c('Latitude','Longitude')


### PROJECT COORDINATES FOR SPATIAL ANALYSES - this now needs WGS84 in CAPITAL LETTERS (not wgs84 anymore!!)
### SpatialPointsDataFrame cannot contain POSIXlt data type!

DataGroup.Wgs <- SpatialPoints(data.frame(tracks$Longitude, tracks$Latitude), proj4string=CRS("+proj=longlat +datum=WGS84"))
input <- SpatialPointsDataFrame(SpatialPoints(data.frame(tracks$Longitude, tracks$Latitude), proj4string=CRS("+proj=longlat +datum=WGS84")), data = tracks, match.ID=F)
DgProj <- CRS(paste("+proj=laea +lon_0=", loc$Longitude[1], " +lat_0=", loc$Latitude[1], sep=""))
DataGroup.Projected <- spTransform(input, CRS=DgProj)
input <- DataGroup.Projected



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SPLIT INTO FORAGING TRIPS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Trips <- input[input@data$Latitude>90,]
for(i in 1:length(unique(tracks$ID)))
  {
  	Temp <- subset(input, ID == unique(tracks$ID)[i])
	  Trip <- tripSplit(Track=Temp, Colony=loc[loc$ID==unique(tracks$ID)[i],2:3], InnerBuff=overview$InnerBuff[overview$DataGroup==dg], ReturnBuff=overview$ReturnBuff[overview$DataGroup==dg], Duration = overview$Duration[overview$DataGroup==dg], plotit=F, nests=F)
  	if(dim(Trips)[1] == 0) {Trips <- Trip[Trip@data$trip_id!="-1",]} else
  	Trips <- rbind(Trips,Trip[Trip@data$trip_id!="-1",])
  }


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE TRIP SUMMARIES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tas<-tripSummary(Trips, Colony=loc)
tas$DataGroup<-dg
TRIP_SUMMARY<-rbind(TRIP_SUMMARY,tas)
fwrite(TRIP_SUMMARY,"S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Seabird_trip_summaries_part4.csv")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PROCEED ONLY WHEN THERE ARE AT LEAST 5 TRIPS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DataGroup <- Trips
if(length(unique(DataGroup@data$trip_id))>4){


### WRITE GOOGLE EARTH FILE FOR EASIER INSPECTION OF OUTLIERS
#outp<-spTransform(Trips, CRS=CRS("+proj=longlat +datum=WGS84"))
#input <- SpatialPointsDataFrame(SpatialPoints(data.frame(outp@data$Longitude, outp@data$Latitude), proj4string=CRS("+proj=longlat +datum=WGS84")), data = outp@data[,1:10], match.ID=F)
#writeOGR(input, sprintf("SeabirdTracks_DataGroup%s.kml",dg), layer=sprintf("SeabirdTracks_DataGroup%s",dg), driver="KML", overwrite_layer=T)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE SCALE FOR AREA RESTRICTED SEARCH
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ScaleOut <- scaleARS(DataGroup, Scales = c(seq(0, 250, 0.5)), Peak="Flexible")


######################################## ORIGINAL DATA SECTION ######################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DELINEATE CORE AREAS (KERNEL DENSITY ESTIMATOR - 50% Utilisation Distribution)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### creating blank output file
SampleLoopOutput <- data.frame(DataGroup = dg, Nind = length(unique(Trips@data$ID)),Ntrips = length(unique(Trips@data$trip_id)))

UD<-50		## pick the % utilisation distribution (50%, 95% etc.)
Output <- batchUDOL(DataGroup=Trips, Scale = ScaleOut/2, UDLev = UD)
SampleLoopOutput$BA<-Output$OverlapIndex


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE 95% MCP AS AREA REQUIRED FOR MANAGEMENT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mcp<-as.numeric(mcp.area(SpatialPoints(coordinates(Trips)), percent = 95,unin = "m",unout = "km2", plotit = T))
SampleLoopOutput$MMA<-mcp


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE AREA OF MARINE IBAs based on 10,15, and 20% threshold AND EXPORT SPATIAL INDEX OUTPUT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

AREAs<-spatInd(Output$UDpolygons, Res = ScaleOut/50)				### resolution in decimal degrees
SampleLoopOutput <- cbind(SampleLoopOutput,AREAs)
ADD_OUT<-merge(SampleLoopOutput,overview,by="DataGroup",all.x=T)
SpatInd_SUMMARY<-rbind(SpatInd_SUMMARY,ADD_OUT)
fwrite(SpatInd_SUMMARY,"S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\SpatIndex_OrigData_all_part4.csv")






######################################## SIMULATED DATA SECTION ######################################################




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### SPECIFY THE SAMPLE SIZES AND CLEAN THE WORKSPACE ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  DataGroup$X <- DataGroup@coords[,1]
  DataGroup$Y <- DataGroup@coords[,2]
  BoundBox <- bbox(DataGroup)
  UIDs <- unique(DataGroup$ID)
  Ntrips <- length(UIDs)
  if(Ntrips<4){Ntrips <- length(unique(DataGroup@data$trip_id))			### fail safe for small data sets with few individuals
			names(DataGroup)[c(10,20)]<-c("Ind","ID")}

### CREATE SERIES THAT IS NOT COMPUTATIONALLY EXCESSIVE
	Nloop<-nloops[nloops<Ntrips]

### CREATE DATA FRAME THAT DETERMINES WHICH LOOP IS RUN HOW OFTEN
  DoubleLoop <- data.frame(SampleSize = rep(Nloop,each=ITERA), Iteration=rep(seq(1:ITERA),length(Nloop)))
  LoopNr <- seq(1:dim(DoubleLoop)[1])
print(sprintf("starting simulation with Datagroup %s",dg))
print(Sys.time())


# remove large objects that are not needed for processing
rm("ADD_OUT","AREAs","bootstrap",'tracks','DataGroup.Wgs','Trips','DataGroup.Projected','input','varianceTest','Nloop','tas','Output','SampleLoopOutput','mcp')
gc(verbose = TRUE)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### ALL CALCULATIONS BELOW HERE ARE REPEATED OVER A RANGE OF SAMPLE SIZES ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#setup parallel backend to use 8 processors
cl<-makeCluster(8)
registerDoParallel(cl)

Result <- foreach(LoopN=LoopNr, .packages=c('vegan',"sp","adehabitatHR","geosphere","rgdal","raster"),.combine = rbind) %dopar% {		#.combine = rbind, 

#Result<-data.frame()
#for(LoopN in LoopNr){

### setting iteration number and number of sampled individuals
N<-DoubleLoop$SampleSize[LoopN]
i<-DoubleLoop$Iteration[LoopN]

### creating blank output file
SampleLoopOutput <- data.frame(SampleSize = N, InclusionMean = 0,Iteration=i)

### selecting the random number of samples
RanNum <- sample(UIDs, N, replace=F)
Selected <- DataGroup[DataGroup$ID %in% RanNum,]
NotSelected <- DataGroup[!DataGroup$ID %in% RanNum,]



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DELINEATE CORE AREAS (KERNEL DENSITY ESTIMATOR - 50% Utilisation Distribution)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
UD<-50		## pick the % utilisation distribution (50%, 95% etc.)
Output <- batchUDOL(DataGroup=Selected, Scale = ScaleOut/2, UDLev = UD)
SampleLoopOutput$BA<-Output$OverlapIndex
SampleLoopOutput$DataGroup<-dg


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE 95% MCP AS AREA REQUIRED FOR MANAGEMENT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mcp<-as.numeric(mcp.area(SpatialPoints(coordinates(Selected)), percent = 95,unin = "m",unout = "km2", plotit = FALSE))
SampleLoopOutput$MMA<-mcp


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ASSESS INCLUSION VALUE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#plot(Output$UDpolygons)
#plot(NotSelected, pch=16, col=2, add=T)
Overlain <- over(NotSelected, Output$UDpolygons)   							## overlay unselected points with 50%UD polygons
SampleLoopOutput$InclusionMean <- length(which(!is.na(Overlain$ID)))/nrow(NotSelected)  	## proportion of unselected locations that are within 50%UD


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE AREA OF MARINE IBAs based on 10,15, and 20% threshold
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

AREAs<-spatInd(Output$UDpolygons, Res = ScaleOut/50)				### resolution in decimal degrees
SampleLoopOutput <- cbind(SampleLoopOutput,AREAs) 
rownames(SampleLoopOutput)<-LoopN
#Result<-rbind(Result,SampleLoopOutput)

rm('Selected','Overlain','AREAs','mcp','Output','NotSelected')
gc(verbose = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# WRITE OUTPUTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
return(SampleLoopOutput)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# END LOOP OVER VARIOUS SAMPLE SIZES AND ITERATIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

} ### end foreach loop over sample sizes

## stop the cluster
stopCluster(cl)
#print(Sys.time())

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs")
OLname<-sprintf("SimulatedIndexOutput_DataGroup%s.csv",dg)
fwrite(Result,OLname, row.names=F, sep=",")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# END SAFETY LOOP FOR SMALL NUMBER OF TRIPS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

} ### end if statement for <5 trips




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# END OF SPECIES-SPECIFIC INDEX CALCULATION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

} ### end of loop across all data groups

########################################

dg
Sys.time()






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXPORT OUTPUT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation")
save.image("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\AggIndex_output_v10_part3.RData")






