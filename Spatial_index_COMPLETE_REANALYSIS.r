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

source("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\IBA\\Analysis\\mIBA_functions_par.r")
source("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Analysis\\seabird_index3.r")

#source("C:\\STEFFEN\\RSPB\\Marine\\IBA\\Analysis\\mIBA_functions_par.r")
#source("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation\\seabird_index3.r")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD OVERVIEW TABLE FOR ALL DATAGROUPS AND SPECIES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Data")
#setwd("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation\\Data")

overview<-fread("Seabird_priority_overview_v2.csv", header=T, sep=',')


overview<- overview %>%
  mutate(ReturnBuff=as.numeric(ReturnBuff)) %>%
  mutate(Duration=as.numeric(Duration)) %>%
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


### CHECK WHAT DATA ALREADY EXIST
#donefiles<-list.files(path="S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs", pattern=".*SeabirdOverlapIndex_DataGroup.*\\.csv")

dgs<-overview$DataGroup							
#dgs<-dgs[-1]
#dgs<-dgs[order(dgs, decreasing=T)]


TRIP_SUMMARY<-data.frame()
SpatInd_SUMMARY<-data.frame()




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
fwrite(TRIP_SUMMARY,"S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Seabird_trip_summaries_all_backend.csv")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PROCEED ONLY WHEN THERE ARE AT LEAST 5 TRIPS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DataGroup <- Trips
if(length(unique(DataGroup@data$trip_id))>4){


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
mcp<-as.numeric(mcp.area(SpatialPoints(coordinates(Trips)), percent = 95,unin = "m",unout = "km2", plotit = FALSE))
SampleLoopOutput$MMA<-mcp


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE AREA OF MARINE IBAs based on 10,15, and 20% threshold AND EXPORT SPATIAL INDEX OUTPUT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

AREAs<-spatInd(Output$UDpolygons, Res = ScaleOut/50)				### resolution in decimal degrees
SampleLoopOutput <- cbind(SampleLoopOutput,AREAs)
ADD_OUT<-merge(SampleLoopOutput,overview,by="DataGroup",all.x=T)
SpatInd_SUMMARY<-rbind(SpatInd_SUMMARY,ADD_OUT)
fwrite(SpatInd_SUMMARY,"S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\SpatIndex_OrigData_all_backend.csv")






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
save.image("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\AggIndex_output_v10_backend.RData")






