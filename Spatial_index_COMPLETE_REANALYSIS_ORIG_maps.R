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

## 26 Aug 2017: modified to re-instate single file writing of results because memory limit reached in DG 36 and 181


## SPLIT INTO NEW SCRIPT FOR ORIGINAL DATA ONLY 13 Sept 2017
## included as new output maps of IBAs (for Maria)
## used a different colony definition for GPS data based on first locations
## used different scaleARS sequence for different families





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
source("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Analysis\\seabird_index5.r")
spatInd

#source("C:\\STEFFEN\\RSPB\\Marine\\IBA\\Analysis\\mIBA_functions_par.r")
#source("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation\\seabird_index3.r")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD OVERVIEW TABLE FOR ALL DATAGROUPS AND SPECIES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Data")
#setwd("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation\\Data")

overview<-fread("Seabird_priority_overview_v5.csv", header=T, sep=',')


overview<- overview %>%
  mutate(ReturnBuff=as.numeric(ReturnBuff)) %>%
  mutate(Duration=as.numeric(Duration)) %>%
  arrange(DataGroup) %>%
  select(DataGroup, Family, scientific_name,site_name,colony_name,LATITUDE,LONGITUDE,breed_stage,n_individuals,n_tracks,device,InnerBuff,ReturnBuff,Duration)

head(overview)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD PREVIOUSLY SAVED RESULTS AND PREPARE OUTPUT DATA FRAMES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs")


### CHECK WHAT DATA ALREADY EXIST
#donefiles<-list.files(path="S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs", pattern=".*SeabirdOverlapIndex_DataGroup.*\\.csv")

dgs<-overview$DataGroup							
#dgs<-dgs[!(dgs<205)]
#dgs<-dgs[order(dgs, decreasing=T)]


TRIP_SUMMARY<-data.frame()
SpatInd_SUMMARY<-data.frame()
AreaIncrease_SUMMARY<-data.frame()



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
	select(Latitude,Longitude)
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
					plotit=T, nests=ifelse(overview$device[overview$DataGroup==dg]=="GPS",T,F))
  	if(dim(Trips)[1] == 0) {Trips <- Trip[Trip@data$trip_id!="-1",]} else
  	Trips <- rbind(Trips,Trip[Trip@data$trip_id!="-1",])
  }


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE TRIP SUMMARIES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tas<-tripSummary(Trips, Colony=loc)
tas$DataGroup<-dg
TRIP_SUMMARY<-rbind(TRIP_SUMMARY,tas)
fwrite(TRIP_SUMMARY,"S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Seabird_trip_summaries_all.csv")



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
fwrite(AreaIncrease_SUMMARY,"S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\AreaIncrease_OrigData_all.csv")


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
fwrite(SpatInd_SUMMARY,"S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\SpatIndex_OrigData_all.csv")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOT MAP OF MARINE IBAs THAT EXCEED 20% THRESHOLD
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs")
pdf(sprintf("Potential_IBA_DG%i_%s.pdf",dg,overview$scientific_name[overview$DataGroup==dg]), width=8, height=8)
IBAcount<-polyCount(Polys=Output$UDpolygons, Res = ScaleOut/50)
title(sprintf("%s | %s | %s",overview$scientific_name[overview$DataGroup==dg],overview$breed_stage[overview$DataGroup==dg],overview$site_name[overview$DataGroup==dg]))
try(IBA<-thresholdRaster(IBAcount, 20), silent=T)
maps::map("worldHires", add=T, fill=T, col="darkolivegreen3")
title(sprintf("%s | %s | %s",overview$scientific_name[overview$DataGroup==dg],overview$breed_stage[overview$DataGroup==dg],overview$site_name[overview$DataGroup==dg]))
dev.off()



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






