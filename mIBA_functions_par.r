## tripSplit    #####################################################################################################################

## Phil Taylor & Mark Miller, 2011

## this script splits central place foraging animal movement data
## into individual trips away from the colony based on distance and time
## away from a defined colony. A distance buffer is set, under which data is
## assumed to be either roosting or device error and is ignored.

## Track must be either a DataFrame or SpatialPointsDataFrame with Latitude, Longitude,
## ID and TrackTime fields
## Colony must be a DataFrame with Latitudes and Longitudes
## InnerBuff is a number indicating the distance in km that must be travelled for the
## movement to be considered a trip
## ReturnBuff is a number indicating the proximity in km that is required for a trip
## to be considered as returning.
## Duration is the length of time, in hours, that the birds must be atlarge for for the
## movement to be considered a trip.
## if plotit = TRUE a map will be drawn.
## if MidPoint = TRUE the calculations will be done projected on the data's centroid.


tripSplit <- function(Track, Colony, InnerBuff = 15, ReturnBuff = 45, Duration = 12, plotit = TRUE, MidPoint = FALSE, nests=FALSE)
  {

  if(!"Latitude" %in% names(Track)) stop("Latitude field does not exist")
  if(!"Longitude" %in% names(Track)) stop("Longitude field does not exist")
  if(!"ID" %in% names(Track)) stop("ID field does not exist")
  if(!"TrackTime" %in% names(Track)) stop ("TrackTime field does not exist")

  if(!"Latitude" %in% names(Colony)) stop("Colony missing Latitude field")
  if(!"Longitude" %in% names(Colony)) stop("Colony missing Longitude field")

  if(!(is.double(InnerBuff) & is.double(ReturnBuff))) stop ("InnerBuff and ReturnBuff should be numbers")

  require(sp)
  require(maps)
  require(mapdata)
  require(maptools)
  require(rgdal)

  if(class(Track) != "SpatialPointsDataFrame")
    {
    Track.Wgs <- SpatialPoints(data.frame(Track$Longitude, Track$Latitude), proj4string=CRS("+proj=longlat + datum=WGS84"))
    Track.Projected <- spTransform(Track.Wgs, CRS=CRS(paste("+proj=laea +lon_0=", Colony$Longitude[1], " +lat_0=", Colony$Latitude[1], sep="")))
    Track <- SpatialPointsDataFrame(Track.Projected, data = Track)
    }

    ### added section by Steffen Oppel to facilitate nest-specific distance calculations######
if(nests == TRUE)
    {  if(!"ID" %in% names(Colony)) stop("Colony missing ID field")
    nest<- Colony[match(unique(Track$ID), Colony$ID),]
    Colony.Wgs <- SpatialPoints(data.frame(nest$Longitude, nest$Latitude), proj4string=CRS("+proj=longlat + datum=WGS84"))
    usedCRS<-CRS(proj4string(Track))
    Colony.Projected <- spTransform(Colony.Wgs, CRS=usedCRS)
} else{

if(MidPoint == FALSE)
    {
    Colony.Wgs <- SpatialPoints(data.frame(Colony$Longitude, Colony$Latitude), proj4string=CRS("+proj=longlat + datum=WGS84"))
    Colony.Projected <- spTransform(Colony.Wgs, CRS=CRS(paste("+proj=laea +lon_0=", Colony$Longitude, " +lat_0=", Colony$Latitude, sep="")))
    } else
    {
   mid_point<-data.frame(centroid(cbind(DataGroup$Longitude, DataGroup$Latitude)))
    Colony.Wgs <- SpatialPoints(data.frame(Colony$Longitude, Colony$Latitude), proj4string=CRS("+proj=longlat + datum=WGS84"))
    Colony.Projected <- spTransform(Colony.Wgs, CRS=CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep="")))
    Track <- spTransform(Track, CRS=CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep="")))### Transform DataGroup too.
        }
    } 		## ends the else loop for nests=FALSE


  Track$X <- Track@coords[,1]
  Track$Y <- Track@coords[,2]

  Track$Returns <- ""
  Track$trip_id <- 0
  Track$ColDist <- spDists(Track, Colony.Projected)
  Trip.Sequence <- 0
  Time.Diff <- 0
  Max.Dist <- 0
  ReturnBuff <- ReturnBuff * 1000
  InnerBuff <- InnerBuff * 1000

  if(plotit == TRUE)
    {
    plot(Track, pch=1, cex=0.5)
    legend("topleft", paste(Track$ID[1]))
    points(Colony.Projected, pch=18, cex=1.5, col=2)
    }

  i <- 0
  while(i < nrow(Track))
    {
    i <- i + 1
  if(Track$ColDist[i] < InnerBuff) {Track$trip_id[i] <- -1} else
    {
    k <- i
    if(i == nrow(Track)) {Track$trip_id[i] <- -1; break}      ### need to look at how these breaks affect the DataGroup loop
    Dist <- Track$ColDist[i]
    while(Dist >= InnerBuff)
      {
      if(k == nrow(Track) & Dist < ReturnBuff) {break} else
       {
      if(k == nrow(Track))
       {
       print(paste("track ", Track$ID[1], Trip.Sequence + 1, " does not return to the colony", sep=""))
       Track$Returns[i:k] <- "N" ; break
       }
       }
      k <- k + 1
      points(Track[k,], col=2, pch=16, cex=0.5)
      Dist <- Track$ColDist[k]
      }
    Time.Diff <- (Track$TrackTime[k] - Track$TrackTime[i]) / 3600
    Max.Dist <- max(Track$ColDist[i:k])
    if(Time.Diff < Duration |  Max.Dist < InnerBuff)
      {
      Track$trip_id[i:k] <- -1;
      i <- k;
      print(paste("trip ", Track$ID[1], Trip.Sequence + 1, " is too small a trip"))
      next
      }
    Trip.Sequence <- Trip.Sequence + 1
    Track$trip_id[i:k] <- paste(Track$ID[1], Trip.Sequence, sep="")
    i <- k
    print(paste(Track$ID[1], Trip.Sequence, sep=""))
    }
    }
    points(Track, pch=16, cex=0.75, col=as.factor(Track$trip_id))
    return(Track)
  }




## tripSummary   #####################################################################################################################

## STEFFEN OPPEL, 2011

## this script provides a simple summary for the foraging trips of central place foraging animals
## direction can be provided from individual nests if desired (nests=TRUE), default is for colony (mean lat and long across all nests if specified)
## output is a table that provides trip length, distance, and direction for each trip


## Trips must be a SpatialPointsDataFrame generated by the tripSplit function
## Colony must be a DataFrame with Latitudes and Longitudes
#### IF nests=TRUE, Colony must be a DataFrame with ID (the same ID as in Trips), Latitudes and Longitudes
### modified 27 December 2016 to make function more robust to different data frame structure
### updated 9 January 2017 to allow non-numeric trip ID


require(geosphere)
                                         
tripSummary <- function(Trips, Colony=Colony, nests=FALSE)
  {

  if(!"Latitude" %in% names(Colony)) stop("Colony missing Latitude field")
  if(!"Longitude" %in% names(Colony)) stop("Colony missing Longitude field")



### SUMMARISE MAX DIST FROM COLONY AND TRIP TRAVELLING TIME FOR EACH TRIP

trip_distances<-data.frame(trip=unique(Trips@data$trip_id), max_dist=0, duration=0, total_dist=0)			### removed as.numeric as this only works with numeric ID
trip_distances<-trip_distances[!(trip_distances$trip==-1),]			### removed as.numeric as this only works with numeric ID
trip_distances$ID<-Trips@data$ID[match(trip_distances$trip, Trips@data$trip_id)]

for (i in unique(trip_distances$trip)){			### removed as.numeric as this only works with numeric ID
x<-Trips@data[Trips@data$trip_id==i,]
maxdist<-cbind(x$Longitude[x$ColDist==max(x$ColDist)],x$Latitude[x$ColDist==max(x$ColDist)])	
if(dim(maxdist)[1]>1){maxdist<-maxdist[1,]}
trip_distances[trip_distances$trip==i,2]<-max(Trips@data$ColDist[Trips@data$trip_id==i])/1000
trip_distances[trip_distances$trip==i,3]<-(max(Trips@data$TrackTime[Trips@data$trip_id==i])-min(Trips@data$TrackTime[Trips@data$trip_id==i]))/3600


## Calculate distances from one point to the next and total trip distance
x$Dist[1]<-x$ColDist[1]/1000				### distance to first point is assumed a straight line from the nest/colony
for (p in 2:dim(x)[1]){
p1<-c(x$Longitude[p-1],x$Latitude[p-1])
p2<-c(x$Longitude[p],x$Latitude[p])
#x$Dist[p]<-pointDistance(p1,p2, lonlat=T, allpairs=FALSE)/1000			### no longer works in geosphere
x$Dist[p]<-distMeeus(p1,p2)/1000						### great circle distance according to Meeus, converted to km

}
trip_distances[trip_distances$trip==i,4]<-sum(x$Dist)+(x$ColDist[p]/1000)				## total trip distance is the sum of all steps plus the dist from the nest of the last location


trip_distances$departure[trip_distances$trip==i]<-format(min(x$DateTime),format="%Y-%d-%m %H:%M:%S") 	## departure time of trip
trip_distances$return[trip_distances$trip==i]<-format(max(x$DateTime),format="%Y-%d-%m %H:%M:%S")		## return time of trip
trip_distances$n_locs[trip_distances$trip==i]<-dim(x)[1]		## number of locations per trip


if(nests == TRUE) {
origin<- Colony[match(unique(x$ID), Colony$ID),]}

origin<-data.frame(mean(Colony$Longitude), mean(Colony$Latitude)) # CHANGED BY MARIA 8DEC14: was "mean(origin$Latitude)"
trip_distances$bearing[trip_distances$trip==i]<-bearing(origin,maxdist)			## great circle route bearing of trip
trip_distances$bearingRhumb[trip_distances$trip==i]<-bearingRhumb(origin,maxdist) 	## constant compass bearing of trip

}

trip_distances<-trip_distances[,c(5,1,6,7,3,2,4,8,9,10)]
return(trip_distances)
}
  

## scaleARS     #####################################################################################################################

## Phil Taylor & Mark Miller, 2012

## scaleARS undertakes First Passage Time (FPT) analysis on each trip in DataGroup (defined
## by the ID field) and identifies the scale at which each trip is interacting with the
## environment, based on the maximum variance in FPT value (Fauchard & Taveraa; Pinuad &
## Weimerskirch). The function relies on Adehabitat package for FPT calculation and returns a
## object of type numeric indicating the average scale across the datagroup.

## DataGroup must be either a DataFrame or SpatialPointsDataFrame with Latitude,
## Longitude and ID as fields.
## Scales should be a vector with the scales (in kilometres) to be tested. Scales
## should be set according to the movements shown in the data (from 1:maximum distance travelled)
## and should not be too regular as this will begin to identify variances in sample rate rather
## than behaviour.
## Peak determines how the scale will be identified for each trip, this must be a character and from
## these options "Flexible", "First", "Max", "User". "Flexible" will make an automated decision,
## "First" will take the first peak in variance, "Max" will take the Maximum variance and "User"
## will allow the user to select the value on the graph.
    
## tested on 27 Dec 2016 with adehabitatLT, no issues encountered

scaleARS <- function(DataGroup, Scales = c(seq(1, 25, 1), seq(30, 50, 5), 75, seq(100, 250, 50)), Peak = "Flexible")
  {


  require(geosphere)
  require(sp)
  require(rgdal)
  require(rgeos)
  require(adehabitatLT)     ### updated to avoid loading deprecated adehabitat - tested and ok on 27 Dec 2016

  if(!"Latitude" %in% names(DataGroup)) stop("Latitude field does not exist")
  if(!"Longitude" %in% names(DataGroup)) stop("Longitude field does not exist")
  if(!"ID" %in% names(DataGroup)) stop("ID field does not exist")
  if(!"TrackTime" %in% names(DataGroup)) stop("TrackTime field does not exist")

  if(class(DataGroup)!= "SpatialPointsDataFrame")     ## convert to SpatialPointsDataFrame and project
    {
    mid_point<-data.frame(centroid(cbind(DataGroup$Longitude, DataGroup$Latitude)))
    DataGroup.Wgs <- SpatialPoints(data.frame(DataGroup$Longitude, DataGroup$Latitude), proj4string=CRS("+proj=longlat + datum=wgs84"))
    DgProj <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep=""))
    DataGroup.Projected <- spTransform(DataGroup.Wgs, CRS=DgProj)
    DataGroup <- SpatialPointsDataFrame(DataGroup.Projected, data = DataGroup)
    }else{DgProj<-DataGroup@proj4string}

  DataGroup$X <- DataGroup@coords[,1]
  DataGroup$Y <- DataGroup@coords[,2]


  DataGrouplt <- as.ltraj(data.frame(DataGroup$X, DataGroup$Y), date=as.POSIXct(DataGroup$TrackTime, origin="1970/01/01", tz="GMT"), id=DataGroup$ID, typeII = TRUE)

  Scales <- Scales * 1000

  fpt.out <- fpt(DataGrouplt, radii = Scales, units = "seconds")
  fpt.scales <- varlogfpt(fpt.out, graph = FALSE)
  Temp <- as.double(fpt.scales[1,])
  plot(Scales, Temp, type="l", ylim=c(0, max(fpt.scales, na.rm=T)))

  ars.scales <- NULL
  UIDs <- unique(DataGroup$ID)
  for(i in 1:length(UIDs))
    {
    if(length(Scales) == length(which(is.na(fpt.scales[i,])))) {print(paste("Warning: ID", UIDs[i], "is smaller than smallest scale and will be ignored")); next}
    Temp <- as.double(fpt.scales[i,])
    #lines(Scales,Temp)
    plot(Scales, Temp, type="l")

    q <- which(!is.na(Temp))
    p <- 2
    while(!is.na(Temp[q[p]]) & Temp[q[p]] < Temp[q[p-1]] & q[p] != length(Temp)) {p <- p + 1}
    while(!is.na(Temp[q[p]]) & Temp[q[p]] > Temp[q[p-1]]) {p <- p + 1}

    rfpt <- Scales[q[p-1]]
    if(suppressWarnings(min(which(is.na(Temp))) == p)) {print(paste("ID", UIDs[i], "has no peak")); next}
    FirstPeak <- Scales[q[p-1]]
    MaxPeak <- Scales[which(Temp == max(Temp[q[p-1]:length(Temp)], na.rm=T))]
    if(Peak == "Flexible")
    {
    if(FirstPeak < MaxPeak[1])
      {
      MaxPeak <- MaxPeak[MaxPeak >= FirstPeak]
      ifelse(MaxPeak[1] < FirstPeak + (max(Scales)/3), ars.sc <- MaxPeak[1], ars.sc <- FirstPeak)
      }  else  {ars.sc <- FirstPeak}
    }
    if(Peak == "Max") {ars.sc <- MaxPeak}
    if(Peak == "First")  {ars.sc <- FirstPeak}
    if(Peak == "User")
    {
    print("Select Peak on Graph")
    N <- identify(Scales, Temp, n=1)
    ars.sc <- Scales[N]
    }
    abline(v=ars.sc, col="red", lty=2)
    ars.scales <- c(ars.scales, ars.sc)
    #print(ars.sc)
    #readline("proceed?")
    }

  AprScale <- mean(ars.scales)
  AprScale <- round(AprScale/1000,3)
  plot((Scales/1000), Temp, type="l", ylim=c(0, max(fpt.scales, na.rm=T)), xlab="Scales (km)", ylab="")
  for(i in 1:length(UIDs))
    {
    Temp <- as.double(fpt.scales[i,])
    lines((Scales/1000),Temp)
    }
  abline(v=ars.scales/1000, col="red", lty=2)
  abline(v=AprScale, col="darkred", lty=1, lwd=3)
  #print(ars.scales)
  #print(AprScale)
  text(max(Scales/1000)/2, 1, paste(AprScale, "km"), col="darkred", cex=3)
  return(AprScale)
  }


## batchUD  #####################################################################################################################

## Phil Taylor & Mark Miller, 2011

## batchUD calculates the Utilisation Distribution for groups of data in
## a larger dataset. The function relies on Adehabitat package for the
## calculations but returns the UD as SpatialPolygons so they can be exported
## as shapefiles and be more versatile in the script.

## DataGroup must be either a DataFrame or SpatialPointsDataFrame with Latitude,
## Longitude and ID as fields. UD will be calculated for each unique ID value and
## each row should be a location. For each UD to be comparable, the data should be
## regularly sampled or interpolated.
## Scale should be the smoothing factor to be used in the Kernel Density Estimation
## and should be provided in Km.
## UDLev should be the quantile to be used for the Utilisation Distribution.
    
## UPDATED BY STEFFEN OPPEL ON 16 Dec 2016 to fix orphaned holes in output geometry
## UPDATED BY STEFFEN OPPEL ON 23 Dec 2016 to switch to adehabitatHR
## UPDATED BY ANA CARNEIRO AND STEFFEN OPPEL 31 JAN 2017 to same4all=T


batchUD <- function(DataGroup, Scale = 50, UDLev = 50)
    {
    require(sp)
    require(maptools)
    require(rgdal)
    require(adehabitatHR)   ### adehabitat is deprecated, switched to adehabitatHR on 23 Dec 2016
    require(geosphere)

    if(!"Latitude" %in% names(DataGroup)) stop("Latitude field does not exist")
    if(!"Longitude" %in% names(DataGroup)) stop("Longitude field does not exist")
    if(!"ID" %in% names(DataGroup)) stop("ID field does not exist")

    if(class(DataGroup)!= "SpatialPointsDataFrame")     ## convert to SpatialPointsDataFrame and project
    {
    mid_point<-data.frame(centroid(cbind(DataGroup$Longitude, DataGroup$Latitude)))
    DataGroup.Wgs <- SpatialPoints(data.frame(DataGroup$Longitude, DataGroup$Latitude), proj4string=CRS("+proj=longlat + datum=wgs84"))
    DgProj <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep=""))
    DataGroup.Projected <- spTransform(DataGroup.Wgs, CRS=DgProj)
    DataGroup <- SpatialPointsDataFrame(DataGroup.Projected, data = DataGroup)
    }else{DgProj<-DataGroup@proj4string}

    DataGroup$X <- DataGroup@coords[,1]
    DataGroup$Y <- DataGroup@coords[,2]
      
      
      ###### ENSURE THAT ONLY TRACKS WITH 5 OR MORE LOCATIONS ARE RETAINED (added 27 Feb 2017) ####
      UIDs <- names(which(table(DataGroup$ID)>5))             ### previous implementation restored
      #DataGroup@data$count<-1
      #nlocs<-aggregate(count~ID,DataGroup@data,sum)
      #retain<-as.character(nlocs$ID[nlocs$count>5])           ### kernelUD will fail for any ID with <5 locations
      DataGroup<-DataGroup[(DataGroup@data$ID %in% UIDs),]
      DataGroup@data$ID<-droplevels(DataGroup@data$ID)        ### encountered weird error when unused levels were retained (27 Feb 2017)

    UIDs <- unique(DataGroup$ID)
    #note<-0          removed loop to check whether >5 data points exist per trip - considered unnecessary
    #KDE.Sp <- NULL
    TripCoords<-SpatialPointsDataFrame(DataGroup, data=data.frame(ID=DataGroup@data$ID,TrackTime=DataGroup@data$TrackTime))		
    TripCoords@data$TrackTime<-NULL
    Ext <- (min(coordinates(TripCoords)[,1]) + 3 * diff(range(coordinates(TripCoords)[,1])))
    if(Ext < (Scale * 1000 * 2)) {BExt <- ceiling((Scale * 1000 * 3)/(diff(range(coordinates(TripCoords)[,1]))))} else {BExt <- 5} #changed from 3 to 5 on 23 Dec 2016 to avoid 'too small extent' error

    ### NEED TO EXPLORE: CREATE CUSTOM GRID TO feed into kernelUD (instead of same4all=T)  
      
KDE.Surface <- adehabitatHR::kernelUD(TripCoords, h=(Scale * 1000), grid=1000, extent=BExt, same4all=TRUE)		## newer version needs SpatialPoints object and id no longer required in adehabitatHR, also removed 'extent' as it caused problems
KDE.Sp <- adehabitatHR::getverticeshr(KDE.Surface, percent = UDLev,unin = "m", unout = "km2")	## syntax differs from older version; THIS FUNCTION CAN FAIL WHEN same4all=FALSE

    
    KDE.Sp@proj4string <- DgProj
    KDE.Wgs <- spTransform(KDE.Sp, CRS=CRS("+proj=longlat +ellps=WGS84"))
    Tbl <- data.frame(Name_0 = rep(1, length(UIDs)), Name_1 = 1:length(UIDs), ID = UIDs)
    row.names(Tbl) <- UIDs
    KDE.Spdf <- SpatialPolygonsDataFrame(KDE.Sp, data=Tbl)

    plot(KDE.Spdf, border=factor(UIDs))
      
      ##### OVERLAYS IN THE polyCount FUNCTION WILL NOT WORK IF THE POLYGONS CONTAIN HOLES OR ARE ORPHANED
      ## simple fix to remove holes from polygon object
      va90a <- spChFIDs(KDE.Spdf, paste(KDE.Spdf$Name_0, KDE.Spdf$Name_1, KDE.Spdf$ID, sep = ""))
      va90a <- va90a[, -(1:4)]
      va90_pl <- slot(va90a, "polygons")
      va90_pla <- lapply(va90_pl, checkPolygonsHoles)
      p4sva <- CRS(proj4string(va90a))
      vaSP <- SpatialPolygons(va90_pla, proj4string = p4sva)
      va90b <- SpatialPolygonsDataFrame(vaSP, data = as(va90a, "data.frame"))   ### this returns an empty data frame
      va90b@data<-KDE.Spdf@data                                                   ### this adds the original data back into the data frame - may not work if entire polygons are removed
      
      
    return(va90b)     ## changed from KDE.Spdf to replace with cleaned version
    }


## varianceTest    ###########################################################################################

## Phil Taylor & Mark Miller, 2012

## this script tests the spatial variance between polygons, investigating
## whether the variance between polygons in a group is significantly different
## to the variance between polygons from different groups.
## groups are determined using the DID field.

## Polys must be a series of polygons of class SpatialPolygonsDataFrame.
## ID must be a vector of equal length the Polys, identifying the group each Polys belongs to.
## Iteration must be an integer indicating the number of times to iterate the random sample.


varianceTest <- function(Polys, DID, Iteration=100)
  {



  require(rgeos)
  require(sp)
  #set.seed(1)
  if(length(Polys) != length(DID)) stop("DID must be the same length as Polys")
  if(class(Iteration) != "numeric") stop("Iteration must be an integer")
  TripTable <- data.frame(Index = 1:length(Polys), DID = DID)
  TripTable$Multiple <- TripTable[,2] %in% names(which(table(TripTable[,2]) > 1))
  if(length(Polys) ==  length(unique(TripTable[,2]))) {stop("Each Poly has a unique DID, Variance is Zero, Polys must be categorised")}
  if(class(Polys) != "SpatialPolygonsDataFrame") stop("Polys must be of class SpatialPolygonsDataFrame")
  Output <- NULL
  for(i in unique(TripTable[TripTable$Multiple == TRUE,2]))
    {
    PValue <- numeric()   #### changed
        for(k in 1:Iteration)
    {
    Selected <- TripTable[TripTable$DID == i,]
    Sample <- Polys[which(TripTable$DID == i),]

    if(nrow(Selected) > length(unique(TripTable[TripTable$DID != i,2])))
      {
      TRand <- sample(nrow(Selected), length(unique(TripTable[TripTable[,2] != i,2])))
      Selected <- Selected[TRand,]
      Sample <- Sample[TRand,]
      }
    Sample <- gSimplify(Sample, tol=100)
    Reps <- 1:nrow(Selected)
    TDist <- NULL
    for(l in Reps)
      {
      Pairs <- Reps[-(which(Reps == l))]
      if(length(which(Pairs > l)) < 1) {break}
      Pairs <- Pairs[which(Pairs > l)]
      for(m in Pairs)
        {
        Temp <- as.double(gDistance(Sample[l,], Sample[m,], hausdorff = T))
        TDist <- rbind(TDist, Temp)
        }
      }
    RanBirds <- sample(unique(TripTable[TripTable[,2] != i,2]), nrow(Selected))
    ITrips <- NULL
    for(l in RanBirds)
      {
      RTrips <- TripTable[TripTable[,2] == l,]
      STrips <- sample(as.character(RTrips[,1]),1)
      ITrips <- c(ITrips, STrips)
      }
    ISample <- Polys[as.double(ITrips),]
    ISample <- gSimplify(ISample, tol=100)
    Reps <- 1:nrow(Selected)
    IDist <- NULL              ### changed location
    for(l in Reps)
      {
      Pairs <- Reps[-(which(Reps == l))]
      if(length(which(Pairs > l)) < 1) {break}
      Pairs <- Pairs[which(Pairs > l)]
      for(m in Pairs)
        {
        Temp <- as.double(gDistance(ISample[l,], ISample[m,], hausdorff = T))
        IDist <- rbind(IDist, Temp)
        }
      }
    MW1 <- suppressWarnings(wilcox.test(TDist, IDist, alternative = "two.sided", paired =F))
    Result <- suppressWarnings(data.frame(DependentValue = TDist[,1], IndependentValue = IDist[,1]))
    boxplot(Result, main= paste("p.value =", round(MW1$p.value, 6)))
    PValue <- c(PValue, round(MW1$p.value, 6))  ### changed

    }

  Temp <- data.frame(DID = i, PVal = mean(PValue))   ### changed
  Output <- rbind(Output, Temp)
  }
  ## The Null Hypothesis is that the within individual DIDs are part of the same
  ## distribution as the between individual.
  ## Rejection level is set at 0.25
  PValueF <- round(mean(Output$PVal, na.rm=T),6)
  par(mfrow=c(1,1), mai=c(0.5,0.5,0.5,0.5))
  hist(Output$PVal, main= paste("p.value =", PValueF), border="grey", breaks=20)
  if(PValueF < 0.25)
    {
    legend("center", PValueF, "Do not reject the Null Hypothesis,
    site fidelity is significant", cex=2, bty="n" )
    } else {legend("center", PValueF, "Reject the Null Hypothesis,
    site fidelity is not significant", cex=2, bty="n" )
    }
  TripTable <- merge(TripTable, Output, all.x = T)
  return(PValueF)
  }





## bootstrap     ######################################################################################################

## Phil Taylor & Mark Miller, 2012

## this script Iteratively subsamples a dataset of tracking data, investigating the
## affect of sample size. At each iteration the data is split, one half is used as the
## 'training' data and the 50%UD is calculated from this. The second half is used as
## 'testing' and the proportion of it, captured within the 50%UD is calculated.
## A perfect dataset would tend towards 0.5.
## By fitting a trend line to this relationship we can establish the sample size at which
## any new data would simply add to the existing knowledge. This script indicates how close to
## this value the inputted data are.

## DataGroup must be a dataframe or SpatialPointsDataFrame with Latitude, Longitude and ID as fields.
## Scale determines the smoothing factor used in the kernel analysis.
## Iteration determines the number of times each sample size is iterated.
    
## REVISED BY Steffen Oppel in 2015 to facilitate parallel processing
## updated to adehabitatHR by Steffen Oppel on 27 Dec 2016
## changed to same4all=TRUE on 4 Feb 2017
    
## REVISED in 2017 to avoid error in nls function of singular gradient
## added mean output for inclusion value even if nls fails

bootstrap <- function(DataGroup, Scale=100, Iteration=50)
  {

  require(sp)
  require(geosphere)
  require(rgdal)
  require(adehabitatHR)   #### NEED TO FIX
  require(foreach)
  require(doParallel)
  require(parallel)

  if(!"Latitude" %in% names(DataGroup)) stop("Latitude field does not exist")
  if(!"Longitude" %in% names(DataGroup)) stop("Longitude field does not exist")
  if(!"ID" %in% names(DataGroup)) stop("ID field does not exist")

  if(class(DataGroup)!= "SpatialPointsDataFrame")     ## convert to SpatialPointsDataFrame and project
    {
    mid_point<-data.frame(centroid(cbind(DataGroup$Longitude, DataGroup$Latitude)))
    DataGroup.Wgs <- SpatialPoints(data.frame(DataGroup$Longitude, DataGroup$Latitude), proj4string=CRS("+proj=longlat + datum=wgs84"))
    DgProj <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep=""))
    DataGroup.Projected <- spTransform(DataGroup.Wgs, CRS=DgProj)
    DataGroup <- SpatialPointsDataFrame(DataGroup.Projected, data = DataGroup)
    }else{DgProj<-DataGroup@proj4string}

   DataGroup$X <- DataGroup@coords[,1]
  DataGroup$Y <- DataGroup@coords[,2]
  BoundBox <- bbox(DataGroup)
  UIDs <- unique(DataGroup$ID)
  Ntrips <- length(UIDs)
  Nloop<- seq(1,(Ntrips-1),ifelse(Ntrips>100,10,1))
  DoubleLoop <- data.frame(SampleSize = rep(Nloop,each=Iteration), Iteration=rep(seq(1:Iteration),length(Nloop)))
  LoopNr <- seq(1:dim(DoubleLoop)[1])	
  UDLev <- 50

#setup parallel backend to use 4 processors
cl<-makeCluster(detectCores())
registerDoParallel(cl)
Result<-data.frame()

Result <- foreach(LoopN=LoopNr, .combine = rbind, .packages=c("sp","adehabitatHR","geosphere","rgdal")) %dopar% {

    N<-DoubleLoop$SampleSize[LoopN]
    i<-DoubleLoop$Iteration[LoopN]
    Coverage <- NULL
    Inclusion <- NULL
    History <- NULL

    Output <- data.frame(SampleSize = N, InclusionMean = 0,Iteration=i)

     RanNum <- sample(UIDs, N, replace=F)
      SelectedCoords <- coordinates(DataGroup[DataGroup$ID %in% RanNum,])
      NotSelected <- DataGroup[!DataGroup$ID %in% RanNum,]
      Temp <- data.frame(SelectedCoords[,1], SelectedCoords[,2])
      Ext <- (min(Temp[,1]) + 3 * diff(range(Temp[,1])))
      if(Ext < (Scale * 1000 * 2)) {BExt <- ceiling((Scale * 1000 * 3)/(diff(range(Temp[,1]))))} else {BExt <- 5}
      Temp <- SpatialPoints(Temp,proj4string=DgProj)      ### added because adehabitatHR requires SpatialPoints object

      KDE.Surface <- adehabitatHR::kernelUD(Temp, h=(Scale * 1000), grid=1000, extent=BExt, same4all=TRUE)		## newer version needs SpatialPoints object and id no longer required in adehabitatHR
    try(      ## inserted to avoid function failing if one iteration in bootstrap has crazy extent
      KDE.UD <- adehabitatHR::getverticeshr(KDE.Surface, percent = UDLev,unin = "m", unout = "km2"), silent = TRUE)			## syntax differs from older version 
      #KDE.Spl <- kver2spol(KDE.UD)     ## deprecated, newer function would be khr2estUDm, but not required
      #KDE.Spl@proj4string <- DgProj    ## no longer necessary after update to adehabitatHR
      Overlain <- over(NotSelected, KDE.UD)   ## changed from KDE.Spl
      Output$InclusionMean <- length(which(!is.na(Overlain$area)))/nrow(NotSelected)  ## updated because overlay will not yield predictable number
    return(Output)
    }

## stop the cluster
stopCluster(cl)

  par(mfrow=c(1,1), mai=c(1,1,1,1))
  #Result <- Output[1:nrow(Output)-1,]
  write.table(Result,"bootout_temp.csv", row.names=F, sep=",")
  try(M1 <- nls((Result$InclusionMean ~ (a*Result$SampleSize)/(1+b*Result$SampleSize)), data=Result, start=list(a=1,b=0.1)), silent = TRUE)
 if ('M1' %in% ls()){       ### run this only if nls was successful
  PredData <- data.frame(SampleSize = unique(Result$SampleSize))
  Result$pred<-predict(M1)
  P2 <- aggregate(pred~SampleSize, Result, FUN=mean)
  P2$sd <- aggregate(InclusionMean~SampleSize, Result, FUN=sd)[,2]
  plot(InclusionMean~SampleSize, data=Result, pch=16, cex=0.2, col="darkgray", ylim=c(0,1), xlim=c(0,nrow(PredData)), ylab="Inclusion", xlab="SampleSize")
  yTemp <- c((P2[,2] + P2[,3]), rev(P2[,2] - P2[,3]))
  xTemp <- c(1:nrow(P2), nrow(P2):1)
  polygon(x=xTemp, y=yTemp, col="gray93", border=F)
  points(InclusionMean~SampleSize, data=Result, pch=16, cex=0.2, col="darkgray")
  lines(P2, lty=1,lwd=2)
  Asymptote <- (summary(M1)$coefficients[1]/summary(M1)$coefficients[2])
  RepresentativeValue <- max(P2$pred)/Asymptote*100
  print(RepresentativeValue)
  text(x=0, y=1,paste(round(RepresentativeValue,2), "%", sep=""), cex=2, col="gray45", adj=0)
}else{RepresentativeValue <- mean(Result$InclusionMean[Result$SampleSize==(Ntrips-1)])}   ### if nls is unsuccessful then use mean output for largest sample size
  Result$RepresentativeValue <- RepresentativeValue
  return(Result)
  }





## polyCount  ######################################################################################################

## Phil Taylor & Mark Miller, 2012

## polyCount calculates the number of overlapping Polygons. The calculation is
## done by overlapping the polygons within a grid and counting the number
## falling within each gridcell. The Res parameter sets the size of the grid
## (in decimal degrees) to be used. Smaller grids will allow for more detailed
## counts but will slow computing time. The function returns a raster object with
## values showing the proportion of Polys (i.e. nPolys/total nPolys) overlapping
## each cell. The raster extent is set at the bounding limits of the Polys or, when
## data crosses the dateline, set to the northern and southern most limits of the
## Polys but longitudinally crossing the circumference of the world.

## Polys must be a SpatialPolygonsDataFrame of the polygons to be counted.
## Res must be a numeric object indicating the resolution in decimal degrees.

## Steffen Oppel revision on 14 Dec 2016 - removed + (Res * 100) from NCol in L. 664
## Steffen Oppel revision on 16 Dec 2016: fixed the point overlay problem by using a proper grid instead
## Steffen Oppel revision on 27 Dec 206: ensured that output was projected in WGS84

## version 1.2    05-04-2012

polyCount <- function(Polys, Res = 0.1)
  {

  require(raster)
  require(maps)

  if(!class(Polys) %in% c("SpatialPolygonsDataFrame", "SpatialPolygons")) stop("Polys must be a SpatialPolygonsDataFrame")
  if(is.na(projection(Polys))) stop("Polys must be projected")

  Poly.Spdf <- spTransform(Polys, CRS=CRS("+proj=longlat +ellps=WGS84"))
  DgProj <- Polys@proj4string

  DateLine <- Poly.Spdf@bbox[1,1] < -178 & Poly.Spdf@bbox[1,2] > 178
  if(DateLine == TRUE) {print("Data crosses DateLine")}

  UDbbox <- bbox(Poly.Spdf)
  if(DateLine == TRUE)  {UDbbox[1,] <- c(-180,180)}
  BL <- floor(UDbbox[,1])                 # + (Res/2) - removed on 16 Dec 2016 because it results in some polygons outside the grid
  TR <- ceiling(UDbbox[,2])
  NRow <- ceiling(sqrt((BL[1] - TR[1])^2)/Res)
  NCol <- ceiling(sqrt((BL[2] - TR[2])^2)/Res) #+ (Res * 100)				### THIS LINE CAUSES PROBLEMS BECAUSE IT GENERATES LATITUDES >90 which will cause spTransform to fail
  Grid <- GridTopology(BL, c(Res,Res), c(NRow, NCol))
    newgrid<-SpatialGrid(Grid, proj4string = CRS("+proj=longlat + datum=wgs84"))
    spol <- as(newgrid, "SpatialPolygons")								### this seems to create an orphaned hole
    SpGridProj <- spTransform(spol, CRS=DgProj)
    GridIntersects <- over(SpGridProj, Polys)
    SpGridProj<- SpatialPolygonsDataFrame(SpGridProj, data = data.frame(ID=GridIntersects$ID, row.names=sapply(SpGridProj@polygons,function(x) x@ID)))
    SpGridProj <- subset(SpGridProj, !is.na(SpGridProj@data$ID))
  #SpGrid <- SpatialPoints(Grid, proj4string = CRS("+proj=longlat + datum=wgs84"))
  #SpdfGrid <- SpatialPointsDataFrame(SpGrid, data.frame(Longitude=SpGrid@coords[,1], Latitude=SpGrid@coords[,2]))
  #SpGridProj <- spTransform(SpdfGrid, CRS=DgProj)
  #GridIntersects <- over(SpGridProj, Polys)
  #SpGridProj@data$Intersects$ID <- GridIntersects$ID
  #SpGridProj <- subset(SpGridProj, !is.na(SpGridProj@data$Intersects$ID))   ### SpGridProj[!is.na(SpGridProj@data$Intersects$ID),] 			###
  plot(SpGridProj)

  Count <- 0
  for(i in 1:length(Polys))
    {
    TempB <- Polys[i,]
    Temp <- over(SpGridProj, TempB)[,1]     ### inserted based on Matthew Carroll's advice; MAY NEED TO SWAP arguments in 'over'?
    Temp[is.na(Temp)] <- 0
    Temp[Temp > 0] <- 1
    Count <- Count + Temp
    #Prop <- Count/i                        ### removed to improve efficiency
    }
  Prop <- Count/length(Polys)     ### removed from loop over polys as it only needs to be calculated once
  #GridIntersects$inside<-as.numeric(as.character(GridIntersects$ID))    ### this only works for numeric trip_id!!
  #GridIntersects$Prop <- 0
  #GridIntersects$Prop[!is.na(GridIntersects$inside)] <- Prop    #[,1] removed based on Matthew Carroll's advice, because fixed in L. 678
  SpGridProj@data$Prop <- Prop
  SpGridOUT <- spTransform(SpGridProj, CRS=CRS("+proj=longlat +ellps=WGS84"))   ### show output in WGS84
  SGExtent <- extent(SpGridOUT)
  RT <- raster(SGExtent, as.double(NCol), as.double(NRow))
  WgsRas <- (rasterize(x=SpGridOUT,y=RT, field = "Prop"))

  plot(WgsRas, asp=1)
  map("world", add=T, fill=T, col="darkolivegreen3")
  projection(WgsRas) <- CRS("+proj=longlat + datum=wgs84")
  return(WgsRas)
  }


## thresholdRaster     ######################################################################################################

## Phil Taylor & Mark Miller, 2012

## thresholdRaster applies a threshold to a raster, and isolates any areas above
## that threshold value. Converting the raster values to polygons is difficult and
## so this part takes some time. The function returns a SpatialPolygonsDataFrame
## containing the polygons that are above threshold, and with an attributes table
## holding each sites Maximum raster value.

## CountRas must be a raster object with values 0 - 1.
## Threshold must be a number indicating the percentage value to be used
## as the threshold.


thresholdRaster <- function(CountRas, Threshold = 10)
    {

    require(raster)
    require(maps)
    require(geosphere)

    plot(CountRas, asp=1)
    map("world", add=T, fill=T, col="darkolivegreen3")
    Threshold <- Threshold/100
    RasSites <- CountRas >= Threshold
    plot(RasSites, asp=1, col=rev(heat.colors(25)))
    map("world", add=T, fill=T, col="darkolivegreen3")

    if(length(which(getValues(CountRas) > Threshold)) < 1)
      {
      Mid <- c(bbox(CountRas)[1,1]+(as.numeric(bbox(CountRas)[1,2])- as.numeric(bbox(CountRas)[1,1]))/2,  midPoint(bbox(CountRas)[,2], bbox(CountRas)[,1])[2])
      text(Mid, "No Site Identified", cex=1.25)
      stop("No cells were above the threshold value")
      }

    Cells <- rasterToPolygons(CountRas, fun=function(x) {x>Threshold})
    DateLine <- Cells@bbox[1,1] < -178 & Cells@bbox[1,2] > 178
    if(DateLine == TRUE)    {Cells <- spTransform(Cells, CRS=DgProj)}

    Sites <- dissolve(Cells)
    ifelse(DateLine == TRUE, projection(Sites) <- DgProj, projection(Sites) <- "+proj=longlat + datum=wgs84")
    Sites <- spTransform(Sites, CRS=CRS("+proj=longlat + datum=wgs84"))
    SiteTable <- data.frame(SiteID = names(Sites), MaxPerc = round(extract(CountRas, Sites, fun=max)*100,2))
    Sites <- SpatialPolygonsDataFrame(Sites, data=SiteTable)
    print(SiteTable)
    return(Sites)
    }

dissolve <- function(Cells)
{
require(sp)
require(rgeos)
require(maptools)
CellsAv <- Cells
plot(Cells)
j <- 0
while(length(CellsAv) > 0)
  {
  j <- j + 1

  Cell1 <- CellsAv[1,]
  CellsAv <- CellsAv[-1,]
  if(length(CellsAv) < 1)
    {
    CellsMerge <- spChFIDs(Cell1, as.character(j))
    if(j == 1) {Sites <- CellsMerge} else
    Sites <- spRbind(Sites, CellsMerge)
    next
    }
  CellsNr <- which(gTouches(Cell1, CellsAv, byid=T))
  if(length(CellsNr) < 1)
    {
    CellsMerge <- spChFIDs(Cell1, as.character(j))
    if(j == 1) {Sites <- CellsMerge} else
    Sites <- spRbind(Sites, CellsMerge)
    next
    }
  CellsSel <- CellsAv[as.double(CellsNr),]
  CellsMerge <- gUnion(Cell1, CellsSel)
  CellsAv <- CellsAv[-as.double(CellsNr),]
  if(length(CellsAv) < 1)
    {
    CellsMerge <- spChFIDs(CellsMerge, as.character(j))
    if(j == 1) {Sites <- CellsMerge} else
    Sites <- spRbind(Sites, CellsMerge)
    next
    }
  CellsNr <- which(gTouches(CellsMerge, CellsAv, byid=T))
  while(length(CellsNr) > 0)
    {
    CellsSel <- CellsAv[as.double(CellsNr),]
    CellsMerge <- gUnion(CellsMerge, CellsSel)
    plot(CellsMerge, add=T, col=2)
    CellsAv <- CellsAv[-as.double(CellsNr),]
    if(length(CellsAv) < 1) break
    CellsNr <- which(gTouches(CellsMerge, CellsAv, byid=T))
    }
  CellsMerge <- spChFIDs(CellsMerge, as.character(j))
  if(j == 1) {Sites <- CellsMerge} else
  Sites <- spRbind(Sites, CellsMerge)
  plot(Sites, add=T, col=2)
  }
plot(Sites, col=names(Sites), add=T)
return(Sites)
}







