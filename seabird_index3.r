
########## SEABIRD SPATIAL AGGREGATION INDEX  ######################################################################################################
### v.3 7 March 2017 - includes sample size loop over each function derived from mIBA bootstrap function
### removed Morisita and Ripley's K index - see v.1 for emd code


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

	### CALCULATE AREA FOR MARINE IBA AND MINIMUM MANAGEMENT AREA ####
  Prop <- Count/length(Polys)     ### removed from loop over polys as it only needs to be calculated once
  SpGridProj@data$Prop <- Prop
  RT <- raster(SpGridProj, ncols=as.double(NCol), nrows=as.double(NRow), vals=Prop)
	
	IBAarea10<-length(Prop[Prop>0.1])*(res(RT)[1]/1000)*(res(RT)[2]/1000)		    ### IBA area with >10% of tracked pop 
	IBAarea125<-length(Prop[Prop>0.125])*(res(RT)[1]/1000)*(res(RT)[2]/1000)		   ### IBA area with >12.5% of tracked pop
	IBAarea20<-length(Prop[Prop>0.2])*(res(RT)[1]/1000)*(res(RT)[2]/1000)		    ### IBA area with >20% of tracked pop 

Cells <- rasterToPolygons(RT)
IBAcells<-Cells[Cells@data$layer>0.1,]

if(dim(coordinates(IBAcells))[1]>4){
mcp10<-as.numeric(mcp.area(SpatialPoints(coordinates(IBAcells)), percent = 100,unin = "m",unout = "km2", plotit = FALSE))}else{mcp10<-IBAarea10}
IBAcells<-Cells[Cells@data$layer>0.125,]
if(dim(coordinates(IBAcells))[1]>4){
mcp125<-as.numeric(mcp.area(SpatialPoints(coordinates(IBAcells)), percent = 100,unin = "m",unout = "km2", plotit = FALSE))}else{mcp125<-IBAarea125}
IBAcells<-Cells[Cells@data$layer>0.2,]
if(dim(coordinates(IBAcells))[1]>4){
mcp20<-as.numeric(mcp.area(SpatialPoints(coordinates(IBAcells)), percent = 100,unin = "m",unout = "km2", plotit = FALSE))}else{mcp20<-IBAarea20}

report<-data.frame(Scale=Res,extent=length(Count),E_W=NCol,N_S=NRow, IBA10=IBAarea10,IBA125=IBAarea125,IBA20=IBAarea20,MCP10=mcp10,MCP125=mcp125,MCP20=mcp20)

	#### CALCULATE STATISTICS FOR SPATIAL AGGREGATION ####
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

##### ERRORS OCCUR IF ANIMALS HAVE ONLY FEW LOCATIONS< SO WE EXCLUDE THEM ###
    UIDs <- names(which(table(DataGroup$ID)>5))
    DataGroup <- DataGroup[DataGroup$ID %in% UIDs,]

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

####
outlist<-list("UDpolygons"=va90b,"OverlapIndex"=olInd)
rm(KDE.Sp,va90b,olInd,OL1,va90a,KDE.Surface,TripCoords,DataGroup,DataGroup.Projected,KDE.Wgs,Tbl,KDE.Spdf,vaSP)
gc(verbose = F)

return(outlist)

    }

