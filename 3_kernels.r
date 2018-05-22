#####################################################
################ LOADING PACKAGES ###################
#####################################################

lu=function (x=x) length(unique(x))
library(rgeos)
library(rgdal)
library(sp)
library(geosphere)
library(adehabitatHR)
library(raster)



#####################################################
######### GENERAL DIRECTIONS AND FILES ##############
#####################################################

## GENERAL DIR
dir <- "C:/Users/ana.carneiro/Documents/GEF analysis"  

## PROJECTIONS
land <- readOGR(dsn="C:/Users/ana.carneiro/Documents/Baselayers/land", layer = "mundo")
DgProj <- CRS("+proj=laea +lon_0=111.686552278329 +lat_0=-88.6447576802893 +ellps=WGS84") ## based on my null grid coords

## TO SAVE KERNEL RESULTS
dir_kernels_v1 <- "C:/Users/ana.carneiro/Documents/GEF analysis/results/3_kernels/breed_stage_v1"
dir_kernels_v2 <- "C:/Users/ana.carneiro/Documents/GEF analysis/results/3_kernels/breed_stage_v2"



#####################################################
################# LOADING SPP DATA ##################
#####################################################

species <- "Southern Royal Albatross"  ## Tristan Albatross  ## "Northern Royal Albatross"

df <- read.csv(paste0(dir, "/results/2_cleaning/", species, ".csv"))
head(df)



#######################################################
################### YEAR QUARTERS #####################
# Q1(Jan-Mar), Q2(Apr-Jun), Q3(July-Sep), Q4(Oct-Dec) #
#######################################################
#######################################################

## THIS IS TO CREATE NB DISTRIBUTIONS BASED ON YEAR QUARTER
## FOR JUVENILES, USE THE WINTERING DISTRIBUTION OF THE NB BIRDS (Q2_Q3)

df$dtime <- as.POSIXct(strptime(df$dtime, "%Y-%m-%d %H:%M:%S")) 
df$quarter <- quarters(df$dtime)
str(df)
df$quarter <- factor(df$quarter)
levels(df$quarter)

df <- subset(df, !quarter %in% "QNA")
df$quarter <- factor(df$quarter)
levels(df$quarter)

num_tracks <- aggregate(df$track_id, list(df$common_name, df$site_name, df$colony_name, df$device, df$age, df$breed_stage), lu)
colnames(num_tracks) <- c("common_name", "site_name", "colony_name", "device", "age", "breed_stage", "n_tracks")
num_tracks


## EXCLUDING DATA (not to be included in the analysis)
## OR KEEP ALL DATA AND CHECK SAMPLE SIZES USED TO CREATE THE KERNEL LATER 
## Tristan Albatross (EXCLUDING)
# df <- subset(df, !device %in% "GLS" | !breed_stage %in% "chick-rearing")
# df <- subset(df, !device %in% "GLS" | !breed_stage %in% "unknown")


## CREATING NB DISTRIBUTIONS PER QUARTER (v1) AND COMBINING WINTERING PART (v2)
df$breed_stage_v1 <-  df$breed_stage
df$breed_stage_v2 <-  df$breed_stage


## DATASETS WITH NO NB DATA WILL PRODUCE AN ERROR - IGNORE
x1 <- subset(df, breed_stage == "non-breeding")
x2 <- subset(df, !breed_stage == "non-breeding")
q <- unique(as.factor(x1$quarter))
A_summer <- data.frame()
A_winter <- data.frame()

for(i in 1:nlevels(q)){
  print(q)
  a <- subset(x1, quarter == q[i])
  a$quarter <- factor(a$quarter)
  if(a$quarter[1] %in% c("Q1", "Q4")){
    a$breed_stage_v1 <- paste0(a$breed_stage[1], "_", a$quarter[1])
    a$breed_stage_v2 <- paste0(a$breed_stage[1], "_", a$quarter[1])
    A_summer <- rbind(A_summer, as.data.frame(a))} else {
      a$breed_stage_v1 <- paste0(a$breed_stage[1], "_", a$quarter[1])
      a$breed_stage_v2 <- paste0(a$breed_stage[1], "_", "Q2_Q3")
      A_winter <- rbind(A_winter, as.data.frame(a))}
}
  
NB <- rbind(A_summer, A_winter)        
df <- rbind(NB, x2)



#######################################################
############## SPECIES METADATA FOR KERNELS ###########
#######################################################

df$device_comb <- df$device
levels(df$device_comb)[levels(df$device_comb)=="GPS"] <- "GPSorPTT"
levels(df$device_comb)[levels(df$device_comb)=="PTT"] <- "GPSorPTT"

num_birds_v1 <- aggregate(df$bird_id, list(df$common_name, df$site_name, df$colony_name, df$device_comb, df$age, df$breed_stage_v1), lu)
colnames(num_birds_v1) <- c("common_name", "site_name", "colony_name", "device_comb", "age", "breed_stage_v1", "n_birds")
num_birds_v1

num_birds_v2 <- aggregate(df$bird_id, list(df$common_name, df$site_name, df$colony_name, df$device_comb, df$age, df$breed_stage_v2), lu)
colnames(num_birds_v2) <- c("common_name", "site_name", "colony_name", "device_comb", "age", "breed_stage_v2", "n_birds")
num_birds_v2
num_birds_v2_sub <- subset(num_birds_v2, grepl("^.+(_Q2_Q3)$", breed_stage_v2))

write.csv(num_birds_v1, paste0(dir, "/metadata/species/", "metadata_", num_birds_v1$common_name[1], "_breed_stage_v1.csv"), row.names = FALSE)
write.csv(num_birds_v2_sub, paste0(dir, "/metadata/species/", "metadata_", num_birds_v2_sub$common_name[1], "_breed_stage_v2.csv"), row.names = FALSE)



#######################################################
################## KERNEL ANALYSIS ####################
#######################################################

## CREATING A NULL GRID (SOUTHERN HEMISPHERE)

## COORDINATES
so <- readOGR(dsn="C:/Users/ana.carneiro/Documents/GEF analysis/metadata/SO points", layer = "SO_points")
so_proj <- spTransform(so, CRS=DgProj)

coords <- so_proj@coords

c <- min(coords[,1])   ## to check my min lon
d <- max(coords[,1])   ## to check my max lon

e <- min(coords[,2])   ## to check my min lat
f <- max(coords[,2])   ## to check my max lat

## by= 10000 means my cell size
a= seq(c, d, by=10000)
b= seq(e, f, by=10000)  
null.grid <- expand.grid(x=a,y=b)
coordinates(null.grid) <- ~x+y
gridded(null.grid) <- TRUE
class(null.grid)
plot(null.grid)


## BASED ON BREEDING STAGE V1

meta <- read.csv(paste0(dir, "/metadata/species/metadata_Southern Royal Albatross_breed_stage_v1.csv"))  ## REPLACE HERE
all_data <- df
head(all_data)
str(all_data)
all_data$breed_stage_v1 <- factor(all_data$breed_stage_v1)
all_data$bird_id <- factor(all_data$bird_id)


## REMOVING GROUPS WITH LESS THAN 5 ROWS
file_to_run <- data.frame()

for (i in 1:nrow(meta)){
  birds_row <- all_data[(all_data$common_name==meta$common_name[i]) & (all_data$colony_name==meta$colony_name[i]) & (all_data$device_comb==meta$device_comb[i]) &
                          (all_data$age==meta$age[i]) & (all_data$breed_stage_v1==meta$breed_stage_v1[i]),]
  birds_row$bird_id <- factor(birds_row$bird_id)
  id <- unique(birds_row$bird_id)
  for (b in 1:nlevels(id)){
    sub_bird <- birds_row[(birds_row$bird_id==id[b]),]
    if (nrow(sub_bird)>5){
    file_to_run <- rbind(file_to_run, as.data.frame(sub_bird))}
  }
}

all_data <- file_to_run


## KERNEL PER GROUP (BASED ON METADATA ROWS WITH V1 FILE)

for (i in 1:nrow(meta)){
  print(i)
  tracks <- all_data[(all_data$common_name==meta$common_name[i]) & (all_data$colony_name==meta$colony_name[i]) & (all_data$device_comb==meta$device_comb[i]) &
                       (all_data$age==meta$age[i]) & (all_data$breed_stage_v1==meta$breed_stage_v1[i]),]
  coordinates(tracks) <- ~coords_lon+coords_lat
  proj4string(tracks) <- DgProj
  tracks$bird_id <- factor(tracks@data$bird_id)
  tracks$device_comb <- factor(tracks$device_comb)
  
  if(tracks$device[1]=="GLS"){
    kerns_gls <- kernelUD(tracks[,11], grid = null.grid, h = 200000)
    stk_gls <- stack(estUDm2spixdf(kerns_gls))
    sum_all_gls_raw <- overlay(stk_gls, fun = mean)
    sum_all_gls_norm <- sum_all_gls_raw/max(getValues(sum_all_gls_raw))
    
    KDERasName_gls_raw = paste(dir_kernels_v1, "/raw/", meta$common_name[i], "_", meta$colony_name[i], "_", meta$device_comb[i], "_", 
                                   meta$age[i], "_", meta$breed_stage_v1[i], "_raw", ".tif", sep="")
    writeRaster(sum_all_gls_raw, filename=KDERasName_gls_raw, format="GTiff", overwrite = TRUE)
    KDERasName_gls_norm = paste(dir_kernels_v1, "/norm/", meta$common_name[i], "_", meta$colony_name[i], "_", meta$device_comb[i], "_", 
                                    meta$age[i], "_", meta$breed_stage_v1[i], "_norm", ".tif", sep="")
    writeRaster(sum_all_gls_norm, filename=KDERasName_gls_norm, format="GTiff", overwrite = TRUE)} else {
      
      kerns_gps_ptt <- kernelUD(tracks[,11], grid = null.grid, h = 50000)
      stk_gps_ptt <- stack(estUDm2spixdf(kerns_gps_ptt))
      sum_all_gps_ptt_raw <- overlay(stk_gps_ptt, fun = mean)
      sum_all_gps_ptt_norm <- sum_all_gps_ptt_raw/max(getValues(sum_all_gps_ptt_raw))
      
      KDERasName_gps_ptt_raw = paste(dir_kernels_v1, "/raw/", meta$common_name[i], "_", meta$colony_name[i], "_", meta$device_comb[i], "_", 
                                     meta$age[i], "_", meta$breed_stage_v1[i], "_raw", ".tif", sep="")
      writeRaster(sum_all_gps_ptt_raw, filename=KDERasName_gps_ptt_raw, format="GTiff", overwrite = TRUE)
      KDERasName_gps_ptt_norm = paste(dir_kernels_v1, "/norm/", meta$common_name[i], "_", meta$colony_name[i], "_", meta$device_comb[i], "_", 
                                      meta$age[i], "_", meta$breed_stage_v1[i], "_norm", ".tif", sep="")
      writeRaster(sum_all_gps_ptt_norm, filename=KDERasName_gps_ptt_norm, format="GTiff", overwrite = TRUE)}
}



## BASED ON BREEDING STAGE V2

meta <- read.csv(paste0(dir, "/metadata/species/metadata_Northern Royal Albatross_breed_stage_v2.csv"))  ## change spp here
all_data <- df
head(all_data)
str(all_data)
all_data_v2 <- subset(all_data, breed_stage_v2 %in% c("non-breeding_Q2_Q3"))
all_data <- all_data_v2
all_data$breed_stage_v2 <- factor(all_data$breed_stage_v2)
all_data$bird_id <- factor(all_data$bird_id)
all_data$colony_name <- factor(all_data$colony_name)
all_data$device_comb <- factor(all_data$device_comb)


## REMOVING GROUPS WITH LESS THAN 5 ROWS

file_to_run <- data.frame()

for (i in 1:nrow(meta)){
  birds_row <- all_data[(all_data$common_name==meta$common_name[i]) & (all_data$colony_name==meta$colony_name[i]) & (all_data$device_comb==meta$device_comb[i]) &
                          (all_data$age==meta$age[i]) & (all_data$breed_stage_v2==meta$breed_stage_v2[i]),]
  birds_row$bird_id <- factor(birds_row$bird_id)
  id <- unique(birds_row$bird_id)
  for (b in 1:nlevels(id)){
    sub_bird <- birds_row[(birds_row$bird_id==id[b]),]
    if (nrow(sub_bird)>5){
      file_to_run <- rbind(file_to_run, as.data.frame(sub_bird))}
  }
}

all_data <- file_to_run


## KERNEL PER GROUP (BASED ON METADATA ROWS WITH V2 FILE)

for (i in 1:nrow(meta)){
  print(i)
  tracks <- all_data[(all_data$common_name==meta$common_name[i]) & (all_data$colony_name==meta$colony_name[i]) & (all_data$device_comb==meta$device_comb[i]) &
                       (all_data$age==meta$age[i]) & (all_data$breed_stage_v2==meta$breed_stage_v2[i]),]
  coordinates(tracks) <- ~coords_lon+coords_lat
  proj4string(tracks) <- DgProj
  tracks$bird_id <- factor(tracks@data$bird_id)
  tracks$device_comb <- factor(tracks$device_comb)
  
  if(tracks$device[1]=="GLS"){
    kerns_gls <- kernelUD(tracks[,11], grid = null.grid, h = 200000)
    stk_gls <- stack(estUDm2spixdf(kerns_gls))
    sum_all_gls_raw <- overlay(stk_gls, fun = mean)
    sum_all_gls_norm <- sum_all_gls_raw/max(getValues(sum_all_gls_raw))
    
    KDERasName_gls_raw = paste(dir_kernels_v2, "/raw/", meta$common_name[i], "_", meta$colony_name[i], "_", meta$device_comb[i], "_", 
                               meta$age[i], "_", meta$breed_stage_v2[i], "_raw", ".tif", sep="")
    writeRaster(sum_all_gls_raw, filename=KDERasName_gls_raw, format="GTiff", overwrite = TRUE)
    KDERasName_gls_norm = paste(dir_kernels_v2, "/norm/", meta$common_name[i], "_", meta$colony_name[i], "_", meta$device_comb[i], "_", 
                                meta$age[i], "_", meta$breed_stage_v2[i], "_norm", ".tif", sep="")
    writeRaster(sum_all_gls_norm, filename=KDERasName_gls_norm, format="GTiff", overwrite = TRUE)} else {
      
      kerns_gps_ptt <- kernelUD(tracks[,11], grid = null.grid, h = 50000)
      stk_gps_ptt <- stack(estUDm2spixdf(kerns_gps_ptt))
      sum_all_gps_ptt_raw <- overlay(stk_gps_ptt, fun = mean)
      sum_all_gps_ptt_norm <- sum_all_gps_ptt_raw/max(getValues(sum_all_gps_ptt_raw))
      
      KDERasName_gps_ptt_raw = paste(dir_kernels_v2, "/raw/", meta$common_name[i], "_", meta$colony_name[i], "_", meta$device_comb[i], "_", 
                                     meta$age[i], "_", meta$breed_stage_v2[i], "_raw", ".tif", sep="")
      writeRaster(sum_all_gps_ptt_raw, filename=KDERasName_gps_ptt_raw, format="GTiff", overwrite = TRUE)
      KDERasName_gps_ptt_norm = paste(dir_kernels_v2, "/norm/", meta$common_name[i], "_", meta$colony_name[i], "_", meta$device_comb[i], "_", 
                                      meta$age[i], "_", meta$breed_stage_v2[i], "_norm", ".tif", sep="")
      writeRaster(sum_all_gps_ptt_norm, filename=KDERasName_gps_ptt_norm, format="GTiff", overwrite = TRUE)}
}




