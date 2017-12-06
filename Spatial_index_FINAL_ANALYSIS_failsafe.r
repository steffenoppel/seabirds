############################################################################################################
####### SPATIAL AGGREGATION INDEX ANALYSIS FOR SEABIRDS   ##################################################
############################################################################################################
## this analysis is based on BirdLife International marine IBA processing scripts
## developed by steffen.oppel@rspb.org.uk in December 2016
## data provided by ana.carneiro@birdlife.org on 12 Dec 2016 and in 2017

## FINAL ANALYSIS FOR ALL SPECIES
## see 'Seabird_spatial_agg_index_v8 - 9' for previous analysis descriptions and notes
## some datasets in the main analysis were found to include erroneous data
## 22 Aug 2017 - complete cleaning of all data with speed filter in 'Seabird_spatial_agg_data_preparation_v3.R'
## 23 Aug 2017 - build entire analysis in one loop, includes raw data analysis, and simulation across sample sizes

## updated on 11 Oct 2017 after feedback from all collaborators
## added single trip analysis suggested by Juan Masello

## troubleshooting on 19 Oct 2017 to re-run for a few data groups
## fixed fatal error on 20 Oct 2017 - need to re-run DG up to 220 for SpatInd_SUMMARYSingle as this has not been saved


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



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD CUSTOM FUNCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\IBA\\Analysis\\mIBA_functions_par.r")
source("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Analysis\\seabird_index6.r")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD OVERVIEW TABLE FOR ALL DATAGROUPS AND SPECIES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Data")

overview<-fread("Seabird_priority_overview_v9.csv", header=T, sep=',')


overview<- overview %>%
  mutate(ReturnBuff=as.numeric(ReturnBuff)) %>%
  mutate(Duration=as.numeric(Duration)) %>%
  filter(n_individuals>4) %>%
  arrange(n_individuals) %>%
  dplyr::select(DataGroup, Family, scientific_name,site_name,colony_name,LATITUDE,LONGITUDE,breed_stage,n_individuals,n_tracks,device,InnerBuff,ReturnBuff,Duration)

head(overview)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD PREVIOUSLY SAVED RESULTS AND PREPARE OUTPUT DATA FRAMES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Outputv15")
dgs<-overview$DataGroup							


### CREATE OUTPUT DATA FRAMES ###
TRIP_SUMMARY<-data.frame()
SpatInd_SUMMARY<-data.frame()
SpatInd_SUMMARYSingle<-data.frame()
AreaIncrease_SUMMARY<-data.frame()
AreaIncrease_SUMMARYSingle<-data.frame()


### IF PROCESS WAS INTERRUPTED, CHECK WHAT DATA ALREADY EXIST AND AVOID RE-RUNNING THE LOOP OVER THESE
TRIP_SUMMARY<-fread("Seabird_trip_summaries_all.csv")
donedgs<-as.numeric(unique(TRIP_SUMMARY$DataGroup))
dgs<-dgs[!(dgs %in% donedgs)]


### READ IN TABLES THAT ALREADY EXIST

TRIP_SUMMARY<-as.data.frame(TRIP_SUMMARY[!(TRIP_SUMMARY$DataGroup %in% dgs),])
SpatInd_SUMMARY<-fread("SpatIndex_OrigData_all.csv")
SpatInd_SUMMARY<-as.data.frame(SpatInd_SUMMARY[!(SpatInd_SUMMARY$DataGroup %in% dgs),])
SpatInd_SUMMARYSingle<-fread("SpatIndex_OrigDataSingle_all.csv")
SpatInd_SUMMARYSingle<-as.data.frame(SpatInd_SUMMARYSingle[!(SpatInd_SUMMARYSingle$DataGroup %in% dgs),])		
AreaIncrease_SUMMARY<-fread("AreaIncrease_OrigData_all.csv")		
AreaIncrease_SUMMARY<-as.data.frame(AreaIncrease_SUMMARY[!(AreaIncrease_SUMMARY$DataGroup %in% dgs),])
AreaIncrease_SUMMARYSingle<-fread("AreaIncrease_OrigDataSingle_all.csv")		
AreaIncrease_SUMMARYSingle<-as.data.frame(AreaIncrease_SUMMARYSingle[!(AreaIncrease_SUMMARYSingle$DataGroup %in% dgs),])






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SET SPATIAL AND ITERATION PARAMETERS FOR BOOTSTRAPPING ANALYSIS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## number of iterations and sample sizes that are randomly drawn
### CREATE SERIES from 3 to max n trips THAT IS NOT COMPUTATIONALLY EXCESSIVE

ITERA=30
nloops<-c(3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,23,26,29,32,35,40,45,50,55,60,65,70,75,85,100,125,150,175,200,250,300,310)## max number is 301 from previous analysis




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# START OF SPECIES-SPECIFIC INDEX CALCULATION FOR EACH DATA GROUP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (dg in dgs){        ### START OVERALL LOOP OVER EACH DATA GROUP (deliberately kept serial to facilitate troubleshooting)

tryCatch({			### try to capture any errors and let the loop continue even if some error occurs

#### SELECT DATA FOR ANALYSIS ##########################
tracksname<-sprintf("CLEAN_SeabirdTracks_DataGroup%s.csv",dg)
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Data\\DataGroups")
tracks<-read.table(tracksname, header=T, sep=",")
tracks$DateTime <- ymd_hms(tracks$DateTime, tz="GMT")
try(tracks$TrackTime<-adjust.duplicateTimes(tracks$TrackTime, tracks$ID),silent=T)


#### REPLACE NON-NUMERIC IDs ##########################

if(is.numeric(tracks$bird_id)==F){tracks$bird_id<-as.numeric(as.factor(tracks$bird_id))}
if(is.numeric(tracks$ID)==F){tracks$ID<-as.numeric(as.factor(tracks$ID))}


  
#### CREATE COLONY LOCATION AND APPROPRIATE COORDINATE REFERENCE SYSTEM FOR PROJECTION ##########################
## FOR GPS DATA - take the first location for each Individual

if(overview$device[overview$DataGroup==dg]=="GPS"){
loc <- tracks %>%
	group_by(ID) %>%
	summarise(Latitude=mean(Latitude[1:3], na.rm=T), Longitude=mean(Longitude[1:3], na.rm=T))
loc<-as.data.frame(loc)

}else{

loc <- overview %>%
	filter(DataGroup==dg) %>%
	mutate(Latitude=LATITUDE, Longitude=LONGITUDE)%>%
	dplyr::select(Latitude,Longitude)
loc<-as.data.frame(loc)
}


### PROJECT COORDINATES FOR SPATIAL ANALYSES - this now needs WGS84 in CAPITAL LETTERS (not wgs84 anymore!!)
### SpatialPointsDataFrame cannot contain POSIXlt data type!

DataGroup.Wgs <- SpatialPoints(data.frame(tracks$Longitude, tracks$Latitude), proj4string=CRS("+proj=longlat +datum=WGS84"))
input <- SpatialPointsDataFrame(SpatialPoints(data.frame(tracks$Longitude, tracks$Latitude), proj4string=CRS("+proj=longlat +datum=WGS84")), data = tracks, match.ID=F)
DgProj <- CRS(paste("+proj=laea +lon_0=", loc$Longitude[1], " +lat_0=", loc$Latitude[1], sep=""))
DataGroup.Projected <- spTransform(input, CRS=DgProj)

### FOR TROUBLESHOOTING
#writeOGR(input, sprintf("SeabirdTracks_DataGroup%s.kml",dg), layer=sprintf("SeabirdTracks_DataGroup%s",dg), driver="KML", overwrite_layer=T)
#overview[overview$DataGroup==dg,]
input <- DataGroup.Projected



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SPLIT INTO FORAGING TRIPS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Trips <- input[input@data$Latitude>90,]
for(i in 1:length(unique(tracks$ID)))
  {
  	Temp <- subset(input, ID == unique(tracks$ID)[i])
	  Trip <- tripSplit(Track=Temp, Colony=loc,
					InnerBuff=overview$InnerBuff[overview$DataGroup==dg],
					ReturnBuff=overview$ReturnBuff[overview$DataGroup==dg],
					Duration = overview$Duration[overview$DataGroup==dg],
					plotit=F, nests=ifelse(overview$device[overview$DataGroup==dg]=="GPS",T,F))
  	if(dim(Trips)[1] == 0) {Trips <- Trip[Trip@data$trip_id!="-1",]} else
  	Trips <- rbind(Trips,Trip[Trip@data$trip_id!="-1",])
  }



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE TRIP SUMMARIES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tas<-tripSummary(Trips, Colony=loc)
tas$DataGroup<-dg
TRIP_SUMMARY<-rbind(TRIP_SUMMARY,tas)
fwrite(TRIP_SUMMARY,"S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Outputv15\\Seabird_trip_summaries_all.csv")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PROCEED ONLY WHEN THERE ARE AT LEAST 5 TRIPS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DataGroup <- Trips
if(length(unique(DataGroup@data$trip_id))>4){



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE SCALE FOR AREA RESTRICTED SEARCH
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### this has a huge impact on the size of IBAs as it influences the size of the kernels
### make this dependent on the max dist from colony for any given dataset

ULscale<-ifelse(max(tas$max_dist)>100,ifelse(max(tas$max_dist)>1000,500,250),max(tas$max_dist))
scaleres<-ifelse(max(tas$max_dist)>50,ifelse(max(tas$max_dist)>1000,1.5,0.5),ifelse(max(tas$max_dist)>25,0.2,0.1))
scalevec = c(seq(0,ULscale,scaleres))
ScaleOut <- scaleARS(DataGroup, Scales = scalevec, Peak="Flexible")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DELINEATE CORE AREAS (KERNEL DENSITY ESTIMATOR - 50% Utilisation Distribution)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### creating blank output file
SampleLoopOutput <- data.frame(DataGroup = dg, Nind = length(unique(Trips@data$ID)),Ntrips = length(unique(Trips@data$trip_id)))
UD<-50		## pick the % utilisation distribution (50%, 95% etc.)
Output <- batchUDOL(DataGroup=Trips, Scale = ScaleOut/2, UDLev = UD)

SampleLoopOutput$BA<-Output$OverlapIndex
ADDareaincrease<-Output$AreaIncrease %>% mutate(DataGroup=dg)
AreaIncrease_SUMMARY<-rbind(AreaIncrease_SUMMARY,ADDareaincrease)
fwrite(AreaIncrease_SUMMARY,"S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Outputv15\\AreaIncrease_OrigData_all.csv")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DELINEATE CORE AREAS (BUT ONLY FOR ONE TRIP PER INDIVIDUAL)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## this function is equivalent to the function above, but uses only 1 trip per individual
## suggested by Juan Masello to avoid pseudoreplication

## select only the first trip_id for each individual
firsttrips<-as.numeric(paste(unique(tas$ID),1,sep=""))
SingleTrips<-Trips[Trips@data$trip_id %in% firsttrips,]

### creating blank output file
SampleLoopOutputSingle <- data.frame(DataGroup = dg, Nind = length(unique(SingleTrips@data$ID)),PropData=round((dim(SingleTrips)[1]/dim(Trips)[1])*100,2))
OutputSingle <- batchUDOL(DataGroup=SingleTrips, Scale = ScaleOut/2, UDLev = UD)

SampleLoopOutputSingle$BA<-OutputSingle$OverlapIndex
ADDareaincreaseSingle<-OutputSingle$AreaIncrease %>% mutate(DataGroup=dg)
AreaIncrease_SUMMARYSingle<-rbind(AreaIncrease_SUMMARYSingle,ADDareaincreaseSingle)
fwrite(AreaIncrease_SUMMARYSingle,"S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Outputv15\\AreaIncrease_OrigDataSingle_all.csv")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE 95% MCP AS AREA REQUIRED FOR MANAGEMENT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mcp<-as.numeric(mcp.area(SpatialPoints(coordinates(Trips)), percent = 95,unin = "m",unout = "km2", plotit = FALSE))
SampleLoopOutput$MMA<-mcp

mcp<-as.numeric(mcp.area(SpatialPoints(coordinates(SingleTrips)), percent = 95,unin = "m",unout = "km2", plotit = FALSE))
SampleLoopOutputSingle$MMA<-mcp


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE AREA OF MARINE IBAs based on 10,15, and 20% threshold AND EXPORT SPATIAL INDEX OUTPUT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

AREAs<-spatInd(Output$UDpolygons, Res = ScaleOut/50)				### resolution in decimal degrees
SampleLoopOutput <- cbind(SampleLoopOutput,AREAs)
ADD_OUT<-merge(SampleLoopOutput,overview,by="DataGroup",all.x=T)
SpatInd_SUMMARY<-rbind(SpatInd_SUMMARY,ADD_OUT)
fwrite(SpatInd_SUMMARY,"S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Outputv15\\SpatIndex_OrigData_all.csv")

AREAs<-spatInd(OutputSingle$UDpolygons, Res = ScaleOut/50)				### resolution in decimal degrees
SampleLoopOutputSingle <- cbind(SampleLoopOutputSingle,AREAs)
ADD_OUTSingle<-merge(SampleLoopOutputSingle,overview,by="DataGroup",all.x=T)
names(SpatInd_SUMMARYSingle)<-names(ADD_OUTSingle)
SpatInd_SUMMARYSingle<-rbind(SpatInd_SUMMARYSingle,ADD_OUTSingle)
fwrite(SpatInd_SUMMARYSingle,"S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Outputv15\\SpatIndex_OrigDataSingle_all.csv")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOT MAP OF MARINE IBAs THAT EXCEED 20% THRESHOLD
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Outputv15")
pdf(sprintf("Potential_IBA_DG%i_%s.pdf",dg,overview$scientific_name[overview$DataGroup==dg]), width=8, height=8)
IBAcount<-polyCount(Polys=Output$UDpolygons, Res = ScaleOut/50)
title(sprintf("%s | %s | %s",overview$scientific_name[overview$DataGroup==dg],overview$breed_stage[overview$DataGroup==dg],overview$site_name[overview$DataGroup==dg]))
try(IBA<-thresholdRaster(IBAcount, 20), silent=T)
maps::map("worldHires", add=T, fill=T, col="darkolivegreen3")
title(sprintf("%s | %s | %s",overview$scientific_name[overview$DataGroup==dg],overview$breed_stage[overview$DataGroup==dg],overview$site_name[overview$DataGroup==dg]))
dev.off()






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
rm("Temp","Trip",'i','datetime2011',"IBA","IBAcount","ADDareaincrease","ADDareaincreaseSingle","firsttrips","ADD_OUT","ADD_OUTSingle","AREAs","bootstrap",'tracks','DataGroup.Wgs','Trips','SingleTrips','DataGroup.Projected','input','varianceTest','Nloop','tas','Output','OutputSingle','SampleLoopOutput','SampleLoopOutputSingle','mcp')
gc(verbose = TRUE)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### ALL CALCULATIONS BELOW HERE ARE REPEATED OVER A RANGE OF SAMPLE SIZES ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#setup parallel backend to use 8 processors
cl<-makeCluster(ifelse(Ntrips>100,4,8))
registerDoParallel(cl)



######### USE AN RBIND LOOP FOR SMALLER DATASETS #########
if(length(TRIP_SUMMARY$trip[TRIP_SUMMARY$DataGroup==dg])<100){

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
write.table(Result,sprintf("SeabirdOverlapIndex_DataGroup%s.csv",dg), row.names=F, sep=",")


}else{						######### USE A WRITE INDIVIDUAL OUTPUT LOOP FOR LARGE SAMPLES #########

Result <- foreach(LoopN=LoopNr, .packages=c('vegan',"sp","adehabitatHR","geosphere","rgdal","raster")) %dopar% {		#.combine = rbind, 

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
Overlain <- over(NotSelected, Output$UDpolygons)   							## overlay unselected points with 50%UD polygons
SampleLoopOutput$InclusionMean <- length(which(!is.na(Overlain$ID)))/nrow(NotSelected)  	## proportion of unselected locations that are within 50%UD


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE AREA OF MARINE IBAs based on 10,15, and 20% threshold
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

AREAs<-spatInd(Output$UDpolygons, Res = ScaleOut/50)				### resolution in decimal degrees
SampleLoopOutput <- cbind(SampleLoopOutput,AREAs) 
rownames(SampleLoopOutput)<-LoopN

rm('Selected','Overlain','AREAs','mcp','Output','NotSelected')
gc(verbose = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# WRITE OUTPUTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Outputv15")
OLname<-sprintf("SeabirdOverlapIndex_DataGroup%s_Part%i.csv",dg,LoopN)
write.table(SampleLoopOutput,OLname, row.names=F, sep=",")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# END LOOP OVER VARIOUS SAMPLE SIZES AND ITERATIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

} ### end foreach loop over sample sizes


} ### end else loop for large samples





## stop the cluster
stopCluster(cl)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# END SAFETY LOOP FOR SMALL NUMBER OF TRIPS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

} ### end if statement for <5 trips




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# END OF SPECIES-SPECIFIC INDEX CALCULATION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


 }, error=function(e){
		cat(sprintf("DataGroup %i failed on %s with %s. \n", dg,Sys.time(),conditionMessage(e), sep=",",collapse=","))})



} ### end of loop across all data groups







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SECOND INSTANCE OF LOOP - JUST IN CASE SOMETHING WENT WRONG!! 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

donedgs<-as.numeric(unique(TRIP_SUMMARY$DataGroup))
dgs<-dgs[!(dgs %in% donedgs)]

if(length(dgs>3)){	### START SECOND INSTANCE OF LOOP ONLY IF SUFFICIENT DATA GROUPS ARE LEFT


for (dg in dgs){        ### START OVERALL LOOP OVER EACH DATA GROUP (deliberately kept serial to facilitate troubleshooting)

tryCatch({			### try to capture any errors and let the loop continue even if some error occurs

#### SELECT DATA FOR ANALYSIS ##########################
tracksname<-sprintf("CLEAN_SeabirdTracks_DataGroup%s.csv",dg)
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Data\\DataGroups")
tracks<-read.table(tracksname, header=T, sep=",")
tracks$DateTime <- ymd_hms(tracks$DateTime, tz="GMT")
try(tracks$TrackTime<-adjust.duplicateTimes(tracks$TrackTime, tracks$ID),silent=T)


#### REPLACE NON-NUMERIC IDs ##########################

if(is.numeric(tracks$bird_id)==F){tracks$bird_id<-as.numeric(as.factor(tracks$bird_id))}
if(is.numeric(tracks$ID)==F){tracks$ID<-as.numeric(as.factor(tracks$ID))}


  
#### CREATE COLONY LOCATION AND APPROPRIATE COORDINATE REFERENCE SYSTEM FOR PROJECTION ##########################
## FOR GPS DATA - take the first location for each Individual

if(overview$device[overview$DataGroup==dg]=="GPS"){
loc <- tracks %>%
	group_by(ID) %>%
	summarise(Latitude=mean(Latitude[1:3], na.rm=T), Longitude=mean(Longitude[1:3], na.rm=T))
loc<-as.data.frame(loc)

}else{

loc <- overview %>%
	filter(DataGroup==dg) %>%
	mutate(Latitude=LATITUDE, Longitude=LONGITUDE)%>%
	dplyr::select(Latitude,Longitude)
loc<-as.data.frame(loc)
}


### PROJECT COORDINATES FOR SPATIAL ANALYSES - this now needs WGS84 in CAPITAL LETTERS (not wgs84 anymore!!)
### SpatialPointsDataFrame cannot contain POSIXlt data type!

DataGroup.Wgs <- SpatialPoints(data.frame(tracks$Longitude, tracks$Latitude), proj4string=CRS("+proj=longlat +datum=WGS84"))
input <- SpatialPointsDataFrame(SpatialPoints(data.frame(tracks$Longitude, tracks$Latitude), proj4string=CRS("+proj=longlat +datum=WGS84")), data = tracks, match.ID=F)
DgProj <- CRS(paste("+proj=laea +lon_0=", loc$Longitude[1], " +lat_0=", loc$Latitude[1], sep=""))
DataGroup.Projected <- spTransform(input, CRS=DgProj)

### FOR TROUBLESHOOTING
#writeOGR(input, sprintf("SeabirdTracks_DataGroup%s.kml",dg), layer=sprintf("SeabirdTracks_DataGroup%s",dg), driver="KML", overwrite_layer=T)
#overview[overview$DataGroup==dg,]
input <- DataGroup.Projected



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SPLIT INTO FORAGING TRIPS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Trips <- input[input@data$Latitude>90,]
for(i in 1:length(unique(tracks$ID)))
  {
  	Temp <- subset(input, ID == unique(tracks$ID)[i])
	  Trip <- tripSplit(Track=Temp, Colony=loc,
					InnerBuff=overview$InnerBuff[overview$DataGroup==dg],
					ReturnBuff=overview$ReturnBuff[overview$DataGroup==dg],
					Duration = overview$Duration[overview$DataGroup==dg],
					plotit=F, nests=ifelse(overview$device[overview$DataGroup==dg]=="GPS",T,F))
  	if(dim(Trips)[1] == 0) {Trips <- Trip[Trip@data$trip_id!="-1",]} else
  	Trips <- rbind(Trips,Trip[Trip@data$trip_id!="-1",])
  }



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE TRIP SUMMARIES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tas<-tripSummary(Trips, Colony=loc)
tas$DataGroup<-dg
TRIP_SUMMARY<-rbind(TRIP_SUMMARY,tas)
fwrite(TRIP_SUMMARY,"S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Outputv15\\Seabird_trip_summaries_all.csv")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PROCEED ONLY WHEN THERE ARE AT LEAST 5 TRIPS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DataGroup <- Trips
if(length(unique(DataGroup@data$trip_id))>4){



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE SCALE FOR AREA RESTRICTED SEARCH
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### this has a huge impact on the size of IBAs as it influences the size of the kernels
### make this dependent on the max dist from colony for any given dataset

ULscale<-ifelse(max(tas$max_dist)>100,ifelse(max(tas$max_dist)>1000,500,250),max(tas$max_dist))
scaleres<-ifelse(max(tas$max_dist)>50,ifelse(max(tas$max_dist)>1000,1.5,0.5),ifelse(max(tas$max_dist)>25,0.2,0.1))
scalevec = c(seq(0,ULscale,scaleres))
ScaleOut <- scaleARS(DataGroup, Scales = scalevec, Peak="Flexible")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DELINEATE CORE AREAS (KERNEL DENSITY ESTIMATOR - 50% Utilisation Distribution)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### creating blank output file
SampleLoopOutput <- data.frame(DataGroup = dg, Nind = length(unique(Trips@data$ID)),Ntrips = length(unique(Trips@data$trip_id)))
UD<-50		## pick the % utilisation distribution (50%, 95% etc.)
Output <- batchUDOL(DataGroup=Trips, Scale = ScaleOut/2, UDLev = UD)

SampleLoopOutput$BA<-Output$OverlapIndex
ADDareaincrease<-Output$AreaIncrease %>% mutate(DataGroup=dg)
AreaIncrease_SUMMARY<-rbind(AreaIncrease_SUMMARY,ADDareaincrease)
fwrite(AreaIncrease_SUMMARY,"S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Outputv15\\AreaIncrease_OrigData_all.csv")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DELINEATE CORE AREAS (BUT ONLY FOR ONE TRIP PER INDIVIDUAL)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## this function is equivalent to the function above, but uses only 1 trip per individual
## suggested by Juan Masello to avoid pseudoreplication

## select only the first trip_id for each individual
firsttrips<-as.numeric(paste(unique(tas$ID),1,sep=""))
SingleTrips<-Trips[Trips@data$trip_id %in% firsttrips,]

### creating blank output file
SampleLoopOutputSingle <- data.frame(DataGroup = dg, Nind = length(unique(SingleTrips@data$ID)),PropData=round((dim(SingleTrips)[1]/dim(Trips)[1])*100,2))
OutputSingle <- batchUDOL(DataGroup=SingleTrips, Scale = ScaleOut/2, UDLev = UD)

SampleLoopOutputSingle$BA<-OutputSingle$OverlapIndex
ADDareaincreaseSingle<-OutputSingle$AreaIncrease %>% mutate(DataGroup=dg)
AreaIncrease_SUMMARYSingle<-rbind(AreaIncrease_SUMMARYSingle,ADDareaincreaseSingle)
fwrite(AreaIncrease_SUMMARYSingle,"S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Outputv15\\AreaIncrease_OrigDataSingle_all.csv")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE 95% MCP AS AREA REQUIRED FOR MANAGEMENT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mcp<-as.numeric(mcp.area(SpatialPoints(coordinates(Trips)), percent = 95,unin = "m",unout = "km2", plotit = FALSE))
SampleLoopOutput$MMA<-mcp

mcp<-as.numeric(mcp.area(SpatialPoints(coordinates(SingleTrips)), percent = 95,unin = "m",unout = "km2", plotit = FALSE))
SampleLoopOutputSingle$MMA<-mcp


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE AREA OF MARINE IBAs based on 10,15, and 20% threshold AND EXPORT SPATIAL INDEX OUTPUT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

AREAs<-spatInd(Output$UDpolygons, Res = ScaleOut/50)				### resolution in decimal degrees
SampleLoopOutput <- cbind(SampleLoopOutput,AREAs)
ADD_OUT<-merge(SampleLoopOutput,overview,by="DataGroup",all.x=T)
SpatInd_SUMMARY<-rbind(SpatInd_SUMMARY,ADD_OUT)
fwrite(SpatInd_SUMMARY,"S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Outputv15\\SpatIndex_OrigData_all.csv")

AREAs<-spatInd(OutputSingle$UDpolygons, Res = ScaleOut/50)				### resolution in decimal degrees
SampleLoopOutputSingle <- cbind(SampleLoopOutputSingle,AREAs)
ADD_OUTSingle<-merge(SampleLoopOutputSingle,overview,by="DataGroup",all.x=T)
names(SpatInd_SUMMARYSingle)<-names(ADD_OUTSingle)
SpatInd_SUMMARYSingle<-rbind(SpatInd_SUMMARYSingle,ADD_OUTSingle)
fwrite(SpatInd_SUMMARYSingle,"S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Outputv15\\SpatIndex_OrigDataSingle_all.csv")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOT MAP OF MARINE IBAs THAT EXCEED 20% THRESHOLD
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Outputv15")
pdf(sprintf("Potential_IBA_DG%i_%s.pdf",dg,overview$scientific_name[overview$DataGroup==dg]), width=8, height=8)
IBAcount<-polyCount(Polys=Output$UDpolygons, Res = ScaleOut/50)
title(sprintf("%s | %s | %s",overview$scientific_name[overview$DataGroup==dg],overview$breed_stage[overview$DataGroup==dg],overview$site_name[overview$DataGroup==dg]))
try(IBA<-thresholdRaster(IBAcount, 20), silent=T)
maps::map("worldHires", add=T, fill=T, col="darkolivegreen3")
title(sprintf("%s | %s | %s",overview$scientific_name[overview$DataGroup==dg],overview$breed_stage[overview$DataGroup==dg],overview$site_name[overview$DataGroup==dg]))
dev.off()






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
rm("Temp","Trip",'i','datetime2011',"IBA","IBAcount","ADDareaincrease","ADDareaincreaseSingle","firsttrips","ADD_OUT","ADD_OUTSingle","AREAs","bootstrap",'tracks','DataGroup.Wgs','Trips','SingleTrips','DataGroup.Projected','input','varianceTest','Nloop','tas','Output','OutputSingle','SampleLoopOutput','SampleLoopOutputSingle','mcp')
gc(verbose = TRUE)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### ALL CALCULATIONS BELOW HERE ARE REPEATED OVER A RANGE OF SAMPLE SIZES ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#setup parallel backend to use 8 processors
cl<-makeCluster(ifelse(Ntrips>100,4,8))
registerDoParallel(cl)



######### USE AN RBIND LOOP FOR SMALLER DATASETS #########
if(length(TRIP_SUMMARY$trip[TRIP_SUMMARY$DataGroup==dg])<100){

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
write.table(Result,sprintf("SeabirdOverlapIndex_DataGroup%s.csv",dg), row.names=F, sep=",")


}else{						######### USE A WRITE INDIVIDUAL OUTPUT LOOP FOR LARGE SAMPLES #########

Result <- foreach(LoopN=LoopNr, .packages=c('vegan',"sp","adehabitatHR","geosphere","rgdal","raster")) %dopar% {		#.combine = rbind, 

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
Overlain <- over(NotSelected, Output$UDpolygons)   							## overlay unselected points with 50%UD polygons
SampleLoopOutput$InclusionMean <- length(which(!is.na(Overlain$ID)))/nrow(NotSelected)  	## proportion of unselected locations that are within 50%UD


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE AREA OF MARINE IBAs based on 10,15, and 20% threshold
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

AREAs<-spatInd(Output$UDpolygons, Res = ScaleOut/50)				### resolution in decimal degrees
SampleLoopOutput <- cbind(SampleLoopOutput,AREAs) 
rownames(SampleLoopOutput)<-LoopN

rm('Selected','Overlain','AREAs','mcp','Output','NotSelected')
gc(verbose = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# WRITE OUTPUTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Outputv15")
OLname<-sprintf("SeabirdOverlapIndex_DataGroup%s_Part%i.csv",dg,LoopN)
write.table(SampleLoopOutput,OLname, row.names=F, sep=",")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# END LOOP OVER VARIOUS SAMPLE SIZES AND ITERATIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

} ### end foreach loop over sample sizes


} ### end else loop for large samples





## stop the cluster
stopCluster(cl)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# END SAFETY LOOP FOR SMALL NUMBER OF TRIPS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

} ### end if statement for <5 trips




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# END OF SPECIES-SPECIFIC INDEX CALCULATION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


 }, error=function(e){
		cat(sprintf("DataGroup %i failed on %s with %s. \n", dg,Sys.time(),conditionMessage(e), sep=",",collapse=","))})



} ### end of loop across all data groups

} ### end IF loop for second try






######################################## PROVIDE INTELLIGIBLE ASSESSMENT OF COMPLETION #######################

completed<-unique(SpatInd_SUMMARY$DataGroup)
stilltodo<-dgs[!(dgs %in% donedgs)]

if(length(stilltodo)>0){
cat(sprintf("THE ANALYSIS HAS NOT COMPLETED!
The process failed with DataGroup %i on %s
The following DataGroups still need to be analysed: %s
Please re-start the process. \n", dg,Sys.time(),paste(c(stilltodo), sep=",",collapse=",")))
}else{
save.image("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Outputv15\\FINAL_ANALYSIS_raw_output.RData")
cat(sprintf("THE ANALYSIS HAS COMPLETED on %s.
Everything is saved. Please shut down the computer. \n", Sys.time()))
}


