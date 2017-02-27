
########## SEABIRD SPATIAL AGGREGATION INDEX  ######################################################################################################



## STEFFEN OPPEL 2016
## This function is based on the marine IBA 'polyCount' function (Lascelles et al. 2016)
## Modified to remove raster output and include polygons with no core ranges
## Added calculation and reporting of Morisita spatial aggregation index
## Requires input from 'batchUD' function, and a resolution (in degrees)



spatInd <- function(Polys, Res = 0.1)
  {
  Poly.Spdf <- spTransform(Polys, CRS=CRS("+proj=longlat +ellps=WGS84"))
  DgProj <- Polys@proj4string

  DateLine <- Poly.Spdf@bbox[1,1] < -178 & Poly.Spdf@bbox[1,2] > 178
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
  #plot(SpGridProj)

  Count <- 0
  for(i in 1:length(Polys))
    {
    TempB <- Polys[i,]
    Temp <- over(SpGridProj, TempB)[,1]     ### inserted based on Matthew Carroll's advice
    Temp[is.na(Temp)] <- 0
    Temp[Temp > 0] <- 1
    Count <- Count + Temp
    }

	#### CALCULATE STATISTICS FOR SPATIAL AGGREGATION ####
report<-data.frame(Scale=sc,extent=length(Count),E_W=NCol,N_S=NRow)
index<- as.numeric(dispindmorisita(Count))
report$Morisita<- index[1]
report$mclu<- index[2]
report$muni<- index[3]
report$imst<- index[4]
report$pchisq<- index[5]
return(report)

  }





########## SEABIRD SPATIAL OVERLAP INDEX  ######################################################################################################

## STEFFEN OPPEL 2016
## This function is based on the marine IBA 'batchUD' function (Lascelles et al. 2016)
## Modified to use adehabitatHR and calculate overlap indices
## returns now a list with the SpatialPolygonsDataFrame and a data.frame with the overlap indices

## updated in February 2017 to include EMD index (Kranstauber et al. 2017)
## required smaller extent, grid, and specified threshold to avoid internal error 9


batchUDOL <- function(DataGroup, Scale = 50, UDLev = 50)
    {
    require(sp)
    require(maptools)
    require(rgdal)
    require(adehabitatHR)
    require(geosphere)
    require(reshape)
    require(move)

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

    UIDs <- unique(DataGroup$ID)


##### USE NEW kernelUD function to bypass loop over individuals ####

TripCoords<-SpatialPointsDataFrame(DataGroup, data=data.frame(ID=DataGroup@data$ID,TrackTime=DataGroup@data$TrackTime))		
TripCoords@data$TrackTime<-NULL
Ext <- (min(coordinates(TripCoords)[,1]) + 3 * diff(range(coordinates(TripCoords)[,1])))
if(Ext < (Scale * 1000 * 2)) {BExt <- ceiling((Scale * 1000 * 3)/(diff(range(coordinates(TripCoords)[,1]))))} else {BExt <- 5}

KDE.Surface <- adehabitatHR::kernelUD(TripCoords, h=(Scale * 1000), grid=1000, extent=BExt, same4all=T)			## NEEDS TO BE same4all=T otherwise overlap will produce rubbish output!
KDE.Sp <- adehabitatHR::getverticeshr(KDE.Surface, percent = UDLev,unin = "m", unout = "km2")	

    UIDs <- names(which(table(DataGroup$ID)>5))
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
      va90b@data<-KDE.Spdf@data                                               ### this adds the original data back into the data frame - may not work if entire polygons are removed
      
      
    #return(va90b)     ## R does not allow return of multiple objects, hence combined in list below


##### CALCULATE VARIOUS OVERLAP INDICES #####

olInd<-data.frame(method=c("HR", "PHR", "VI", "BA", "UDOI", "HD"), mean_index=0, n_comps=0)
for(ix in c("HR", "PHR", "VI", "BA", "UDOI", "HD")){
OL1<-adehabitatHR::kerneloverlaphr(KDE.Surface, method = ix,percent= UDLev, conditional = FALSE)		### set conditional=T will make overlap smaller because it sets everything to 0 outside overlap zone

## CALCULATE MEAN BUT REMOVE DIAGONAL
#str(OL1)
OL1<-as.data.frame(melt(OL1))
OL1<-OL1[!(OL1$X1==OL1$X2),]
olInd$mean_index[olInd$method==ix]<-median(OL1$value)		## mean gives crazy values because occasionally index is >100
olInd$n_comps[olInd$method==ix]<-dim(OL1)[1]
}


##### CALCULATE EARTH MOVERS DISTANCE INDEX #####
#KDE.Small <- adehabitatHR::kernelUD(TripCoords, h=(Scale * 1000), grid=100, same4all=T)	## REDUCED GRID TO LIMIT COMPUTATION TIME
#udspdf <- estUDm2spixdf(KDE.Small)
#all<-stack(udspdf)
emdthresh<-sqrt((diff(range(coordinates(TripCoords)[,1]))^2)+(diff(range(coordinates(TripCoords)[,2]))^2))	## needed to avoid internal error 9
#emdout<-emd(all, threshold=emdthresh)



for (gr in c(100,75,50,25,20)){
rm(emdout)
KDE.Small <- adehabitatHR::kernelUD(TripCoords, h=(Scale * 1000), grid=gr, same4all=T)	## REDUCED GRID TO LIMIT COMPUTATION TIME
udspdf <- estUDm2spixdf(KDE.Small)
all<-stack(udspdf)
withTimeout({try(emdout<-emd(all, threshold=emdthresh), silent=T)}, timeout=500, onTimeout="silent")
if('emdout' %in% ls()){break}
}
print(paste(sprintf("EMD succeed on grid %s",gr)))
add<-data.frame(method=c("EMD"), mean_index=mean(emdout), n_comps=length(emdout))
olInd<-rbind(olInd,add)
rm(gr,emdout)




####
outlist<-list("UDpolygons"=va90b,"OverlapIndex"=olInd)
return(outlist)


    }

