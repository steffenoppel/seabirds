############################################################################################################
####### DEVELOP SPATIAL AGGREGATION INDEX FOR SEABIRDS   ###################################################################
############################################################################################################
## this analysis is based on BirdLife International marine IBA processing scripts
## developed by steffen.oppel@rspb.org.uk in December 2016
## data provided by ana.carneiro@birdlife.org on 12 Dec 2016
## fixed some mIBA_script issues on 14 Dec 2016
## added function to remove duplicate time values on 14 Dec 2016
## initial phase to explore approaches

## v2 on 16 Dec 2016 includes major update/fix of marine IBA scripts for polyCount and batchUD - see marineIBA_function_fixes.r for attempts
## 20 Dec 2016: fixed projection problem by changing wgs84 to WGS84; FIXED polyCount function and NA problem in Morisita input values

## v3 on 22 Dec 2016 includes revision of polyCount function to return all cells and not just occupied cells
## removed parallel computation because it resulted in downstream processing issues and time gain was marginal or negative

## v4 on 22 Dec 2016 improves batchUD function to use adehabitatHR and various overlap indices
## updated v4 on 25 Dec to include final ranking and switch overlap index from mean to median

## updated on 6 Jan 2016 to fix overlap index error (based on same4all=F)

## v5 on 17 Feb 2017 included EMD index from move package

## UPDATE 3 MARCH 2017: CHANGE BOOTSTRAP FUNCTION TO INCLUDE SPATIAL INDICES - REMOVED EMD AND SPATIAL REPLICATION OF MORISITA

## v7 on 7 MARCH 2017: after discussion with Richard Phillips decided on area and BA as meaningful indices. Need to include tracking period and focus only on adult breeding season tracks



### UPDATE 5 JUNE 2017: INCLUDED ALL DATA AND FIXED LOOP




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
library(move)
library(R.utils)
require(foreach)
require(doParallel)
require(parallel)

source("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\IBA\\Analysis\\mIBA_functions_par.r")
source("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Analysis\\seabird_index3.r")

source("C:\\STEFFEN\\RSPB\\Marine\\IBA\\Analysis\\mIBA_functions_par.r")
source("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation\\Analysis\\seabird_index3.r")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SET ITERATION PARAMETERS FOR BOOTSTRAPPING
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ITERA=25



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD DATA FROM BIRDLIFE DATABASE AND MODIFY DATA TO MEET REQUIREMENTS FOR PROCESSING
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Load data from database
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation")
#setwd("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation")
alldat<-read.table("Data\\Steffen_Spp_Prioritisation_Project_2016-12-12.csv", header=T, sep=",")
head(alldat)
alldat2<-read.table("Data\\RawData_2017-04-12.csv", header=T, sep=",")
head(alldat2)



alldat<-rbind(alldat,alldat2)
head(alldat)
dim(alldat)


### Convert Dates and Times and calculate the time difference (in s) between each location and the first location of that deployment
alldat$Time<-format(alldat$time_gmt,format="%H:%M:%S")
alldat$Date<-as.Date(alldat$date_gmt,format="%d/%m/%Y")
alldat$Loctime<-as.POSIXlt(paste(alldat$Date, alldat$Time), format = "%Y-%m-%d %H:%M:%S")
alldat$DateTime <- as.POSIXct(strptime(alldat$Loctime, "%Y-%m-%d %H:%M:%S"), "GMT")
alldat$TrackTime <- as.double(alldat$DateTime)


### FAME DATA SET HAS DIFFERENT DATE FORMAT AND NEEDS TO BE TREATED SEPERATELY
alldat3<-read.table("Data\\FAME_datasets_2017-06-01.csv", header=T, sep=",")
head(alldat3)
alldat4<-read.table("Data\\Jonathan_Louise_datasets_2017-06-01.csv", header=T, sep=",")
head(alldat4)
alldat5<-read.table("Data\\data_JeroenCreuwels_591592.csv", header=T, sep=",")
head(alldat5)
alldat3<-rbind(alldat3,alldat4, alldat5[,1:20])
alldat3$Time<-format(alldat3$time_gmt,format="%H:%M:%S")
alldat3$Date<-as.Date(alldat3$date_gmt,format="%Y-%m-%d")
alldat3$Loctime<-as.POSIXlt(paste(alldat3$Date, alldat3$Time), format = "%Y-%m-%d %H:%M:%S")
alldat3$DateTime <- as.POSIXct(strptime(alldat3$Loctime, "%Y-%m-%d %H:%M:%S"), "GMT")
alldat3$TrackTime <- as.double(alldat3$DateTime)
names(alldat)
names(alldat3)
alldat<-rbind(alldat,alldat3)



alldat<-alldat[,c(1,2,4,6:10,12:14,19,18,24:25)]
names(alldat)[c(8,13:12)]<-c('ID','Latitude','Longitude')
head(alldat)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# REMOVE DATA FROM OUTSIDE THE BREEDING SEASON 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
removal<-c("migration","non-breeding")
alldat<-alldat[!(alldat$breed_stage %in% removal),]
dim(alldat)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# REMOVE DATA FROM OUTSIDE THE BREEDING SEASON 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
removal<-c("migration","non-breeding")
alldat<-alldat[!(alldat$breed_stage %in% removal),]
dim(alldat)



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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# REMOVE DATA SETS THAT ARE TOO SMALL (<5 individuals)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
head(overview)
dim(overview)
overview<-overview[overview$n_individuals>4,]
dim(overview)


alldat<-merge(alldat,overview[,c(1:4,7)],by=c("scientific_name","site_name","age","breed_stage"), all.x=T)

alldat<-alldat[!is.na(alldat$Longitude),]
alldat<-alldat[!is.na(alldat$Latitude),]
dim(alldat)


### CREATE SPATIAL SCALES FOR ANALYSIS
### order not ascending so that parallel loop will not calculate the most demanding subsets simultaneously
spatscales<-exp(seq(0,9.5,0.5))			## on log scale
spatscales<-c(2.5,500,1000,1500,5,100,150,250,7.5,10,12.5,15,17.5,20,25,30,40,50,75)			## in km
#spatscales<-c(500,30,1000,1500,40,75,100,250,5,10,15,20,50)			## in km
spatscales<-spatscales/100			### roughly in decimal degrees, as required for the polyCount function


### Create Table that includes one line for each DataGroup at each spatial scale
### THIS WILL NEED A MANUAL ADJUSTMENT FOR THE TripSPlit function

outI<-expand.grid(overview$DataGroup,spatscales)			### sets up the lines for which we need to calculate Morisita's I
names(outI)<-c("DataGroup","Scale")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# COMBINE OVERVIEW TABLE WITH FAMILY INFO AND SET TRIP SPLIT PARAMETERS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation")
#setwd("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation")

species_list<-read.table("Seabird_SpeciesList.csv", header=T, sep=',')
species_list$GISname<-paste(substr(species_list$Family,1,3), species_list$SpcRecID, sep='')
head(species_list)
head(overview)

overview$Family<-species_list$Family[match(overview$scientific_name, species_list$Scientific_name)]
overview$Seabird_type<-species_list$Seabird_type[match(overview$scientific_name, species_list$Scientific_name)]
overview$device<-alldat$device[match(overview$scientific_name, alldat$scientific_name)]

overview$InnerBuff<-ifelse(overview$device=="PTT",10,ifelse(overview$Seabird_type=="Pelagic seabird", 5,0.5))
overview$ReturnBuff<-ifelse(overview$device=="PTT",50,ifelse(overview$Seabird_type=="Pelagic seabird", 50,5))
overview$Duration<-ifelse(overview$device=="PTT",10,ifelse(overview$Seabird_type=="Pelagic seabird", 5,1))

fix(overview)







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# START OF SPECIES-SPECIFIC INDEX CALCULATION FOR EACH DATA GROUP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation")

OUTPUT<-data.frame()							### set up blank frame to collate TIME AT SEA DATA
#SUMMARY<-data.frame()							### set up blank frame to collate inflection point summaries
OL_INDEX<-data.frame()							### set up blank frame to collect overlap indices for each data set

dim(alldat)
alldat<-alldat[!is.na(alldat$DataGroup),]					### this removes the datasets of <5 individuals
dim(alldat)


### for re-start after some have completed ###
#OL_INDEX<-read.table("SeabirdPrioritisation_Areas_v8.csv", header=T, sep=",")
#OUTPUT<-read.table("SeabirdPrioritisation_TrackDurations.csv", header=T, sep=",")
#head(OUTPUT)
#OUTPUT$start <- as.POSIXct(strptime(OUTPUT$start, "%Y-%m-%d %H:%M:%S"), "GMT")
#OUTPUT$return <- as.POSIXct(strptime(OUTPUT$return, "%Y-%m-%d %H:%M:%S"), "GMT")

### specify the datasets that are worth analysing
dgs<-unique(alldat$DataGroup)
done<-unique(OL_INDEX$DataGroup)
dgs<-dgs[!(dgs %in% done)]



for (dg in dgs){


#### SELECT DATA FOR ANALYSIS AND ORDER BY TIME AND REMOVE DUPLICATE TIME STAMPS ##########################
tracks<-alldat[alldat$DataGroup==dg,]
tracks<-tracks[order(tracks$ID, tracks$TrackTime),]
tracks$TrackTime<-adjust.duplicateTimes(tracks$TrackTime, tracks$ID)


#### REMOVE MYSTERIOUSLY CREATED NA LINES
dim(tracks)
tracks<-tracks[!is.na(tracks$Longitude),]
tracks<-tracks[!is.na(tracks$Latitude),]
dim(tracks)

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
DgProj <- CRS(paste("+proj=laea +lon_0=", loc$Longitude[1], " +lat_0=", loc$Latitude[1], sep=""))
DataGroup.Projected <- spTransform(input, CRS=DgProj)
input <- DataGroup.Projected
localmap<-spTransform(countriesLow, CRS=DgProj) 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SPLIT INTO FORAGING TRIPS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Trips <- input[input@data$Latitude>90,]
for(i in 1:length(unique(tracks$ID)))
  {
  	Temp <- subset(input, ID == unique(tracks$ID)[i])
	Trip <- tripSplit(Track=Temp, Colony=loc[loc$ID==unique(tracks$ID)[i],2:3], InnerBuff=overview$InnerBuff[overview$DataGroup==dg], ReturnBuff=overview$ReturnBuff[overview$DataGroup==dg], Duration = overview$Duration[overview$DataGroup==dg], plotit=T, nests=F)
  	if(dim(Trips)[1] == 0) {Trips <- Trip[Trip@data$trip_id!="-1",]} else
  	Trips <- rbind(Trips,Trip[Trip@data$trip_id!="-1",])
  }


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE TIME AT SEA FOR EACH INDIVIDUAL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tas<-ddply(Trips@data, c('ID','trip_id'), summarise, time_at_sea=mean(as.numeric(difftime(max(DateTime),min(DateTime),units='days'))),start=min(DateTime), return=max(DateTime))
tas$DataGRoup<-dg
OUTPUT<-rbind(OUTPUT,tas)
write.table(OUTPUT,"SeabirdPrioritisation_TrackDurations.csv", row.names=F, sep=",")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE SCALE FOR AREA RESTRICTED SEARCH
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DataGroup <- Trips
ScaleOut <- scaleARS(DataGroup, Scales = c(seq(0, 250, 0.5)), Peak="Flexible")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### ALL CALCULATIONS BELOW HERE ARE REPEATED OVER A RANGE OF SAMPLE SIZES ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  DataGroup$X <- DataGroup@coords[,1]
  DataGroup$Y <- DataGroup@coords[,2]
  BoundBox <- bbox(DataGroup)
  UIDs <- unique(DataGroup$ID)
  Ntrips <- length(UIDs)
  if(Ntrips<4){Ntrips <- length(unique(DataGroup$data$trip_id))			### fail safe for small data sets with few individuals
			names(DataGroup)[c(10,20)]<-c("Ind","ID")}
  Nloop<- seq(3,(Ntrips-1),ifelse(Ntrips>100,10,1))
  DoubleLoop <- data.frame(SampleSize = rep(Nloop,each=ITERA), Iteration=rep(seq(1:ITERA),length(Nloop)))
  LoopNr <- seq(1:dim(DoubleLoop)[1])	
  UDLev <- 50
print(sprintf("starting with Datagroup %s",dg))
print(Sys.time())
#setup parallel backend to use 4 processors
cl<-makeCluster(detectCores())
registerDoParallel(cl)

Result <- foreach(LoopN=LoopNr, .combine = rbind, .packages=c("sp","adehabitatHR","geosphere","rgdal","raster")) %dopar% {

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
return(SampleLoopOutput)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# END LOOP OVER VARIOUS SAMPLE SIZES AND ITERATIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

} ### end foreach loop over sample sizes

## stop the cluster
stopCluster(cl)
print(Sys.time())

#### COMBINE THE RESULTS FOR A DATA GROUP ####
OL_INDEX<-rbind(OL_INDEX,Result)
write.table(OL_INDEX,"SeabirdPrioritisation_Areas_v8.csv", row.names=F, sep=",")

} ### end loop over species 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# END OF SPECIES-SPECIFIC INDEX CALCULATION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



########################################

dg
Sys.time()







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOTTING THE OUTPUT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

head(OL_INDEX)
species_list<-list('1'="RAZO",'2'="RAZO",'3'="FRIG",'4'="FRIG",'5'="RBTB",'6'="RBTB",'7'="SHAG",'8'="SHAG",'9'="SOAL",'10'="MUPE",'11'="MABO",'12'="MABO",'13'="MABO",'14'="MABO",'15'="CGUIL",'16'="CGUIL")
species_labeller <- function(variable,value){return(species_list[value])}


pdf("Seabird_prioritisation_InclusionValue.pdf", height=7, width=9)
ggplot(OL_INDEX, aes(x=SampleSize, y=InclusionMean), size=1)+
geom_smooth(fill="lightblue", size=1.5)+facet_wrap("DataGroup", ncol=4, scales = "fixed", labeller=species_labeller)+
  geom_point(colour="black", size=1.5) +
  xlab("N tracked individuals") +
  ylab("prop. of untracked locations in 50%UD") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())
dev.off()



pdf("Seabird_prioritisation_MMA.pdf", height=7, width=9)
ggplot(OL_INDEX, aes(x=SampleSize, y=log(MMA)), size=1)+
geom_smooth(fill="lightblue", size=1.5)+facet_wrap("DataGroup", ncol=4, scales = "fixed", labeller=species_labeller)+
  geom_point(colour="black", size=1.5) +
  xlab("N tracked individuals") +
  ylab("log(Size) of marine management area (sq km)") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())
dev.off()


pdf("Seabird_prioritisation_IBA10.pdf", height=7, width=9)
ggplot(OL_INDEX, aes(x=SampleSize, y=log(IBA10)), size=1)+
geom_smooth(fill="lightblue", size=1.5)+facet_wrap("DataGroup", ncol=4, scales = "fixed", labeller=species_labeller)+
  geom_point(colour="black", size=1.5) +
  xlab("N tracked individuals") +
  ylab("log(Size) of marine management area (sq km)") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())
dev.off()


pdf("Seabird_prioritisation_MCP10.pdf", height=7, width=9)
ggplot(OL_INDEX, aes(x=SampleSize, y=log(MCP10)), size=1)+
geom_smooth(fill="lightblue", size=1.5)+facet_wrap("DataGroup", ncol=4, scales = "fixed", labeller=species_labeller)+
  geom_point(colour="black", size=1.5) +
  xlab("N tracked individuals") +
  ylab("log(Size) of marine management area (sq km)") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())
dev.off()


pdf("Seabird_prioritisation_BA.pdf", height=7, width=9)
ggplot(OL_INDEX, aes(x=SampleSize, y=BA), size=1)+
geom_smooth(fill="lightblue", size=1.5)+facet_wrap("DataGroup", ncol=4, scales = "fixed", labeller=species_labeller)+
  geom_point(colour="black", size=1.5) +
  xlab("N tracked individuals") +
  ylab("Bhattacharya's affinity index") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())
dev.off()







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXPORT OUTPUT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation")
save.image("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\AggIndex_output_v8.RData")









