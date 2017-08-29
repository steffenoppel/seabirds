############################################################################################################
####### DATA PREPARATION FOR SPATIAL AGGREGATION INDEX FOR SEABIRDS   ######################################
############################################################################################################
## this analysis is based on BirdLife International marine IBA processing scripts
## developed by steffen.oppel@rspb.org.uk in December 2016
## data provided by ana.carneiro@birdlife.org on 12 Dec 2016
## fixed some mIBA_script issues on 14 Dec 2016
## added function to remove duplicate time values on 14 Dec 2016
## initial phase to explore approaches

### UPDATE 18 JULY: included more data and fixed some memory issues because the loop collapsed
## for (thing in ls()) { message(thing); print(object.size(get(thing)), units='auto') }
## identified and removed large objects before loop processing

### SPLIT FROM ANALYTICAL SCRIPT ON 19 JULY 2017 - analysis would stall, so data need to be saved on harddrive and removed from workspace


### UPDATE 4 and 7 AUGUST 2017: DataGroups 72 and 97 (old PTT data from albatrosses in 1990s) required manual elimination of some individuals with very few locations

### UPDATE 17 AUGUST 2017: discovered that site_name and colony_name are very different, so that FAME data are split up into 'Scotland' etc., which leads to low overlap and huge areas.

### UPDATE 22 AUGUST 2017: multiple data sets have poor quality locations, hence included a speedfilter to remove the worst locations




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD PACKAGES AND CUSTOM SCRIPTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(Hmisc)
require(foreign)
require(maptools)
require(geosphere)
require(sp)
library(rgdal)
library(trip)
library(rworldmap)
library(vegan)
library(adehabitatLT)
data(countriesLow)
library(data.table)
library(lubridate)
library(tidyverse)
require(foreach)
require(doParallel)
require(parallel)
#install.packages("devtools")
#require(devtools)
#install_github("TakahiroShimada/SDLfilter")
library(SDLfilter)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD DATA FROM BIRDLIFE DATABASE AND MODIFY DATA TO MEET REQUIREMENTS FOR PROCESSING
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Load data from database
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation")
#setwd("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation")
alldat<-fread("Data\\Steffen_Spp_Prioritisation_Project_2016-12-12.csv", header=T, sep=",")
head(alldat)
alldat2<-fread("Data\\RawData_2017-04-12.csv", header=T, sep=",")
head(alldat2)
alldat<-rbind(alldat,alldat2)
head(alldat)
dim(alldat)

alldat <- alldat %>%
	mutate(DateTime=dmy_hms(paste(date_gmt, time_gmt))) %>%
	mutate(TrackTime=as.double(DateTime))


### FAME DATA SET HAS DIFFERENT DATE FORMAT AND NEEDS TO BE TREATED SEPERATELY
alldat3<-fread("Data\\FAME_datasets_2017-06-01.csv", header=T, sep=",")
head(alldat3)
alldat4<-fread("Data\\Jonathan_Louise_datasets_2017-06-01.csv", header=T, sep=",")
head(alldat4)
alldat5<-fread("Data\\data_JeroenCreuwels_591592.csv", header=T, sep=",")
head(alldat5)
alldat3<-rbind(alldat3,alldat4, alldat5[,1:20])

alldat3 <- alldat3 %>%
	mutate(DateTime=ymd_hms(paste(date_gmt, time_gmt))) %>%
	mutate(TrackTime=as.double(DateTime))

alldat<- alldat %>%
	bind_rows(alldat3) %>%
	dplyr::select(c(1,2,4:10,12:14,19,18,20:22))

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
# REMOVE DATA SETS THAT ARE TOO SMALL (<5 individuals)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Data")
#setwd("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation\\Data")
overview<-fread("Seabird_priority_overview_v2.csv")
head(overview)

alldat<-merge(alldat,overview[,c(1:5,8)],by=c("scientific_name","site_name","colony_name","age","breed_stage"), all.x=T)
dim(alldat)

ALL<- alldat %>%
	filter(!is.na(Longitude)) %>%
	filter(!is.na(Latitude)) %>%
	filter(!is.na(DataGroup))	

dim(ALL)
names(ALL)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOT THE DATA TO CHECK THAT IT LOOKS OK
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

xlow<-min(ALL$Longitude)+0.2
xup<-max(ALL$Longitude)-0.2
yup<-max(ALL$Latitude)-0.5
ylow<-min(ALL$Latitude)+0.5

windows(600,400)
plot(Latitude~Longitude, data=ALL, pch=16, cex=0.3,col=dataset_id, asp=1, xlim=c(xlow,xup), ylim=c(ylow,yup), main="", frame=F, axes=F, xlab="", ylab="")
plot(countriesLow, col='darkgrey', add=T)






##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SPECIFY THE DATA GROUPS THAT NEED CLEANING AND REANALYSIS AFTER CLEANING
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

REANALYSIS <- c(4,8,27,29,30,32,33,34,36,38,39,41,45,48,51,72,73,91,111,118,121,122,141,156,158,159,167,170,172)






##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXPORT INPUT DATA FOR EACH NEW DATAGROUP INTO SEPARATE CSV FILE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm('alldat','alldat2','alldat3','alldat4','alldat5','countriesLow')
gc(verbose = TRUE)


dgs<-overview$DataGroup
pttqi<-data.frame(argos_lc=c("3","2","1","0","A","B","Z"), qi=c(10,8,6,3,2,1,0))


#setup parallel backend to use 4 processors
cl<-makeCluster(8)
registerDoParallel(cl)

foreach(dg=TODO, .packages=c("tidyverse","SDLfilter","adehabitatHR","geosphere","rgdal","sp","adehabitatLT","trip","data.table")) %dopar% {


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






############## FIND THE FAILURE ###########

donefiles<-list.files(path="S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Data\\DataGroups", pattern=".*CLEAN_SeabirdTracks_DataGroup.*\\.csv")

m <- gregexpr('[0-9]+',donefiles)
TODO<-dgs[!(dgs %in% as.numeric(regmatches(donefiles,m)))]


