
########## SEABIRD SPATIAL AGGREGATION INDEX  ######################################################################################################
### v.3 7 March 2017 - includes sample size loop over each function derived from mIBA bootstrap function
### removed Morisita and Ripley's K index - see v.1 for emd code
### v. 4 12 Sept 2017, includes area increase by individual
### v. 5 15 Sept 2017, removed MCP and streamlined the polyCount function

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
o <- over(newgrid, Poly.Spdf, returnList=TRUE, minDimension=0)
ct <- sapply(o, nrow)
#summary(ct)
SGDF <- SpatialGridDataFrame(newgrid, data=data.frame(prop=ct/length(Polys)))
#SGDF.proj <- spTransform(SGDF, DgProj)

#spplot(SGDF, "ct", col.regions=bpy.colors(20))
RT <- raster(SGDF)
Cells <- rasterToPolygons(RT)

plot(RT)
plot(Poly.Spdf, col=2, add=T)


### CALCULATE AREAS OF IBAs
try(IBAarea10<-sum(areaPolygon(Cells[Cells@data$prop>=0.01,]))/1000000,silent=T)    	### IBA area with >10% of tracked pop
try(IBAarea125<-sum(areaPolygon(Cells[Cells@data$prop>=0.125,]))/1000000,silent=T)	### IBA area with >12.5% of tracked pop
try(IBAarea20<-sum(areaPolygon(Cells[Cells@data$prop>=0.2,]))/1000000,silent=T)	    	### IBA area with >20% of tracked pop

report<-data.frame(Scale=Res,extent=length(ct),E_W=NCol,N_S=NRow,
		IBA10=ifelse('IBAarea10' %in% ls(),IBAarea10,0),
		IBA125=ifelse('IBAarea125' %in% ls(),IBAarea125,0),
		IBA20=ifelse('IBAarea20' %in% ls(),IBAarea20,0))
rm(IBAarea10,IBAarea125,IBAarea20)



#### CALCULATE STATISTICS FOR SPATIAL AGGREGATION ####
index<- as.numeric(dispindmorisita(ct))
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
## updated in MARCH 2017 to remove EMD index but include bootstrap functionality


batchUDOL <- function(DataGroup, Scale = 10, UDLev = 50)
    {
    require(sp)
    require(maptools)
    require(rgdal)
    require(adehabitatHR)
    require(geosphere)
    require(reshape)

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
    if(is.factor(DataGroup@data$ID)==T){DataGroup@data$ID<-droplevels(DataGroup@data$ID)} ### encountered weird error when unused levels were retained (27 Feb 2017)

##### ERRORS OCCUR IF ANIMALS HAVE ONLY FEW LOCATIONS, SO WE EXCLUDE THEM ###
    UIDs <- names(which(table(DataGroup$ID)>5))
    DataGroup <- DataGroup[DataGroup$ID %in% UIDs,]


##### IF SUBSAMPLING LEAVES TOO FEW DATA STOP THE LOOP ####

if(length(UIDs)>1){



##### USE NEW kernelUD function to bypass loop over individuals ####

TripCoords<-SpatialPointsDataFrame(DataGroup, data=data.frame(ID=DataGroup@data$ID,TrackTime=DataGroup@data$TrackTime))		
TripCoords@data$TrackTime<-NULL
Ext <- (min(coordinates(TripCoords)[,1]) + 3 * diff(range(coordinates(TripCoords)[,1])))
if(Ext < (Scale * 1000 * 2)) {BExt <- ceiling((Scale * 1000 * 3)/(diff(range(coordinates(TripCoords)[,1]))))} else {BExt <- 5}
KDE.Surface <- adehabitatHR::kernelUD(TripCoords, h=(Scale * 1000), grid=1000, extent=BExt, same4all=T)			## NEEDS TO BE same4all=T otherwise overlap will produce rubbish output!

##### THIS FUNCTION CAN FAIL WHEN GRID IS TOO SMALL ###
try(KDE.Sp <- adehabitatHR::getverticeshr(KDE.Surface, percent = UDLev,unin = "m", unout = "km2"), silent=T)


##### IF GETVERTICES FAILED REPEAT UNTIL GRID IS LARGE ENOUGH ###
### this is very time consuming
### problem may occur from PTT data with very few locations in original data set
### interpolated data are all regularly spaced so kernels are impossible - need to find a way to exclude such 'tracks'

if(!('KDE.Sp' %in% ls())){
KDE.Surface <- adehabitatHR::kernelUD(TripCoords, h=(Scale * 1000), grid=1000, extent=0.2, same4all=T)			## NEEDS TO BE same4all=T otherwise overlap will produce rubbish output!
try(KDE.Sp <- adehabitatHR::getverticeshr(KDE.Surface, percent = UDLev,unin = "m", unout = "km2"), silent=T)
}


if(!('KDE.Sp' %in% ls())){
KDE.Surface <- adehabitatHR::kernelUD(TripCoords, h=(Scale * 1000), grid=500, extent=0.2, same4all=T)			## NEEDS TO BE same4all=T otherwise overlap will produce rubbish output!
try(KDE.Sp <- adehabitatHR::getverticeshr(KDE.Surface, percent = UDLev,unin = "m", unout = "km2"), silent=T)
}


if(!('KDE.Sp' %in% ls())){
KDE.Surface <- adehabitatHR::kernelUD(TripCoords, h=(Scale * 1000), grid=2000, extent=1.2, same4all=T)			## NEEDS TO BE same4all=T otherwise overlap will produce rubbish output!
try(KDE.Sp <- adehabitatHR::getverticeshr(KDE.Surface, percent = UDLev,unin = "m", unout = "km2"), silent=T)
}
	

    #UIDs <- names(which(table(DataGroup$ID)>5))
    KDE.Sp@proj4string <- DgProj
    KDE.Wgs <- spTransform(KDE.Sp, CRS=CRS("+proj=longlat +ellps=WGS84"))
    Tbl <- data.frame(Name_0 = rep(1, length(UIDs)), Name_1 = 1:length(UIDs), ID = UIDs)
    row.names(Tbl) <- UIDs
    KDE.Spdf <- SpatialPolygonsDataFrame(KDE.Sp, data=Tbl)

    #plot(KDE.Spdf, border=factor(UIDs))
      
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


##### CALCULATE BHATTACHARYS AFFINITY INDICES #####
									### removed others to speed up process
OL1<-adehabitatHR::kerneloverlaphr(KDE.Surface, method = "BA",percent= UDLev, conditional = FALSE)		### set conditional=T will make overlap smaller because it sets everything to 0 outside overlap zone

## CALCULATE MEAN BUT REMOVE DIAGONAL
OL1<-as.data.frame(melt(OL1))
OL1<-OL1[!(OL1$X1==OL1$X2),]
olInd<-median(OL1$value)		## mean gives crazy values because occasionally index is >100


#### CALCULATE AREA INCREASE FOR EACH INDIVIDUAL
## added after StatsClub meeting on 9 Sept 2017

      ## pick start individual
      #pollistnames<-as.character(va90b@data$ID)    ## extract the names of the UD polygons for all individuals
      #startINDS<-OL1$X1[which(OL1$value==max(OL1$value))][1] ## extract the names for the individuals with the greatest area overlap
      areaByID<-data.frame(ID=as.character(va90b@data$ID),area=sapply(slot(va90b, "polygons"), slot, "area"), step=1, tot.area=0)
	
	## start with the individual that has the smallest area ##
	startINDS<-as.character(areaByID$ID[which(areaByID$area==min(areaByID$area))])[1]
      startpolygon<-va90b[va90b@data$ID==startINDS,]
      remainIDs<-UIDs[!(UIDs==startINDS)]
      areaByID$tot.area[areaByID$ID==startINDS]<-areaByID$area[areaByID$ID==startINDS]
	

  ## loop over sample size of all individuals ###          
  for (al in 2:length(UIDs)){

      ## assess area increase for each added individual ###
	stepincrease<-data.frame(ID=remainIDs, area=0)
      for (lp in remainIDs) {
        selinds<-c(as.character(startINDS),lp)      
        mergedpol<-unionSpatialPolygons(va90b[va90b@data$ID %in% selinds,], IDs=selinds)
        stepincrease$area[stepincrease$ID==lp]<-area(mergedpol)    
      }

      stepINDS<-as.character(stepincrease$ID[which(stepincrease$area==min(stepincrease$area))])[1]		### add one individual that leads to smallest area increase
      areaByID$tot.area[areaByID$ID==stepINDS]<-stepincrease$area[stepincrease$ID==stepINDS]
      areaByID$step[areaByID$ID==stepINDS]<-al
	remainIDs<-remainIDs[!(remainIDs==stepINDS)]

  }
      
      
      
#### POOL THE OUTPUT
outlist<-list("UDpolygons"=va90b,"OverlapIndex"=olInd,"AreaIncrease"=areaByID)
rm(remainIDs,KDE.Sp,va90b,olInd,OL1,va90a,KDE.Surface,TripCoords,DataGroup,DataGroup.Projected,KDE.Wgs,Tbl,KDE.Spdf,vaSP, stepincrease,mergedpol,startpolygon,selinds, stepINDS)
gc(verbose = F)
return(outlist)


} 			### exit n>9 data point condition

    }			### exit function loop

