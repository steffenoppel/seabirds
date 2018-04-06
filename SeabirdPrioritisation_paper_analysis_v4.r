############################################################################################################
####### SEABIRD SPACE USE - ANALYSES AND SUMMARIES FOR PAPER #################################
############################################################################################################

## written by steffen.oppel@rspb.org.uk in October 2017
## colony data provided by Lizzie Pearmain (BirdLife)
## last update 1 Dec 2017 after comments by John Croxall
## removed IBA figures and modified Table 1 and Fig. 2
## added new Fig. 1 (boxplot for range) and re-labelled all other figures

## reverted Fig. 3 to simple as the fig. suggested by John Croxall looks shit
## changed sequence of incubation and brood-guard based on Ana Carneiro's suggestion of sequence

## changed on 27 Dec 2017 to remove non-Atlantic species
## changed on 29 Dec 2017 to adopt Juan Masello's comments


### revised on 2 Feb 2018 after including policy comments and Jonathan Green's suggestions
### removed exploratory code that did not make it into the final manuscript draft (v7)

### finalised on 4 March 2018 after co-author comments on manuscript
### removed DG32 (pre-breeding) Great Shearwaters
### removed body mass from analysis (Ewan Wakefield and Paulo Catry)
### added BA to Table 2 (Rob Ronconi)
### added supplementary table with BA alternatives


### revised on 20 March to combine brood-guard and post-guard


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD PACKAGES AND CUSTOM SCRIPTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)
library(dplyr)
require(foreign)
library(readr)
library(data.table)
library(lme4)
library(RColorBrewer)
library(forcats)
library(tidyr)
library(ggplot2)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD SAVED RESULTS FROM DATA GROUP SPECIFIC ANALYSIS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## this analysis is based on Spatial_index_FINAL_ANALYSIS.r
## output compiled by SeabirdPrioritisation_data_aggregation_for_paper.R

setwd("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Outputv15")
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Outputv15")
load("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Outputv15\\AggIndex_compiled_output_v17.RData")
load("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Outputv15\\AggIndex_compiled_output_v17.RData")


### POOL BROOD AND POST GUARD INTO ONE STAGE
ORIG$breed_stage<-ifelse(ORIG$breed_stage %in% c("brood-guard","post-guard"),"chick-rearing",ORIG$breed_stage)
SIMUL$breed_stage<-ifelse(SIMUL$breed_stage %in% c("brood-guard","post-guard"),"chick-rearing",SIMUL$breed_stage)
TRIPS$breed_stage<-ifelse(TRIPS$breed_stage %in% c("brood-guard","post-guard"),"chick-rearing",TRIPS$breed_stage)

### use breed stage as brood-guard or other
ORIG$chick<-ifelse(ORIG$breed_stage=="chick-rearing",1,0)
SIMUL$chick<-ifelse(SIMUL$breed_stage=="chick-rearing",1,0)
TRIPS$chick<-ifelse(TRIPS$breed_stage=="chick-rearing",1,0)


setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Data")
#setwd("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation\\Data")
overview<-fread("Seabird_priority_overview_v16.csv")
head(overview)
overview$breed_stage<-ifelse(overview$breed_stage %in% c("brood-guard","post-guard"),"chick-rearing",overview$breed_stage)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# REMOVE DATA FROM OUTSIDE THE ATLANTIC OCEAN AND FROM PRE-BREEDING
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

overview<-overview %>% filter(!(scientific_name %in% c("Fregata minor","Pterodroma ultima")))
overview<-overview %>% filter(!(DataGroup==32))

ORIG<-ORIG %>% filter(DataGroup %in% overview$DataGroup)
SIMUL<-SIMUL %>% filter(DataGroup %in% overview$DataGroup)
AREA<-AREA %>% filter(DataGroup %in% overview$DataGroup)
TRIPS<-TRIPS %>% filter(DataGroup %in% overview$DataGroup)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# READ IN SamplingRate DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Data")
setwd("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation\\Data")

SR<-fread("Seabird_DataGroup_SamplingRates.csv")
overview<-merge(overview, SR[,c(1,10)], by="DataGroup", all.x=T)
head(overview)

ORIG<-ORIG %>% mutate(SamplingRate=as.numeric(overview$SamplingRate[match(DataGroup,overview$DataGroup)]))
SIMUL<-SIMUL %>% mutate(SamplingRate=as.numeric(overview$SamplingRate[match(DataGroup,overview$DataGroup)]))
AREA<-AREA %>% mutate(SamplingRate=as.numeric(overview$SamplingRate[match(DataGroup,overview$DataGroup)]))
TRIPS<-TRIPS %>% mutate(SamplingRate=as.numeric(overview$SamplingRate[match(DataGroup,overview$DataGroup)]))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# READ IN COLONY SIZE DATA AND MODIFY TO COMMON CURRENCY (pairs)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Data")
setwd("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation\\Data")

colonies<-fread("Seabird_priority_colony_sizes.csv")
colonies<-colonies %>% select(DataGroup,scientific_name,PopYearStart,PopYearEnd,PopMin,PopMax,PopUnits,MostRecent)
head(colonies)
unique(colonies$PopUnits)

colonies<-colonies %>%
	dplyr::arrange(desc(MostRecent)) %>%
	mutate(COL_SIZE=ifelse(PopUnits %in% c('breeding pairs','chicks'), PopMin, PopMin*0.5)) %>%
	group_by(DataGroup,scientific_name) %>%
	summarise(COL_SIZE=dplyr::first(COL_SIZE))

overview<-merge(overview, colonies, by=c("DataGroup","scientific_name"), all.x=T)
head(overview)


### MISSING COLONY SIZE DATA ###

misscol<-overview %>% filter(is.na(COL_SIZE))
misscol<-misscol %>% group_by(DataGroup,Family,scientific_name,site_name,colony_name,LATITUDE,LONGITUDE) %>%
	summarise(N_tracks=sum(n_tracks))
#fwrite(misscol,"ColonySizes_NEEDED.csv")




ORIG<-ORIG %>% mutate(COL_SIZE=as.numeric(overview$COL_SIZE[match(DataGroup,overview$DataGroup)]))
SIMUL<-SIMUL %>% mutate(COL_SIZE=as.numeric(overview$COL_SIZE[match(DataGroup,overview$DataGroup)]))
AREA<-AREA %>% mutate(COL_SIZE=as.numeric(overview$COL_SIZE[match(DataGroup,overview$DataGroup)]))
TRIPS<-TRIPS %>% mutate(COL_SIZE=as.numeric(overview$COL_SIZE[match(DataGroup,overview$DataGroup)]))





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SIMPLE SUMMARIES FOR TEXT IN RESULTS SECTION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\MANUSCRIPTS\\in_prep\\SeabirdPrioritisation")
setwd("C:\\STEFFEN\\MANUSCRIPTS\\in_prep\\SeabirdPrioritisation")





dim(ORIG)
length(unique(ORIG$DataGroup))  			### number of datasets analysed
length(unique(ORIG$scientific_name))  		### number of species used
dim(TRIPS)							### number of foraging trips

TRIPS$ID<-paste(TRIPS$DataGroup,TRIPS$ID, sep="_")
length(unique(TRIPS$ID))				### number of individual birds



#### ASSESS AVERAGE NUMBER OF IND PER DATASET

meanNind<-overview %>%
  filter(DataGroup %in% unique(SIMUL$DataGroup)) %>%
  summarise(average=round(mean(n_individuals),0), min=round(min(n_individuals),0), max=round(max(n_individuals),0))
meanNind



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATING TABLE 1: OVERVIEW TABLE OF USED DATASETS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TABLE1<- ORIG %>%
  group_by(Family, breed_stage) %>%
  summarise(n_datasets=length(unique(DataGroup))) %>%
  spread(key=breed_stage, value=n_datasets, fill=0)
TABLE1add<- ORIG %>%
  group_by(Family) %>%
  summarise(n_species=length(unique(scientific_name)))
TABLE1<- merge(TABLE1,TABLE1add, by="Family")

## add proportion of species for each family
## manually adjusted to ATLANTIC SPECIES only based on John Croxall comments


#head(species_list)
#TABLE1add2<-species_list %>% mutate(count=1) %>%
#  group_by(Family) %>%
#  summarise(n_species=sum(count)) %>%
#  filter(Family %in% TABLE1$Family) %>%
#  mutate(prop_family=round((TABLE1$n_species/n_species)*100,1)) %>%
#  select(Family,prop_family)
#TABLE1<- merge(TABLE1,TABLE1add2, by="Family")

n_ATLANTIC_species<-c(6,6,4,5,42,2,11,37,9,5)
TABLE1$prop_family<- (TABLE1$n_species/n_ATLANTIC_species)*100

#TABLE1<-TABLE1[,c(1,6,7,2:5)]
TABLE1<-TABLE1[,c(1,5,6,4,3,2)]

#fwrite(TABLE1,"S:\\ConSci\\DptShare\\SteffenOppel\\MANUSCRIPTS\\in_prep\\SeabirdPrioritisation\\Table1_chick.csv")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATING TABLE 2: TRIP SUMMARY ACROSS SPECIES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## modified after John Croxall suggested to include breed stage
## modified to use median based on Juan Masello's suggestion


### CHECK CERTAIN DATASETS FOR ISSUES
TRIPS %>% filter(scientific_name=="Hydroprogne caspia") %>%
select(ID,DataGroup,scientific_name,breed_stage,duration,max_dist,total_dist) %>%
arrange(duration)


### SUMMARISE TABLE

TABLE2<- TRIPS %>%
  	group_by(Family, scientific_name, breed_stage) %>%
  	summarise(median_duration=round(median(duration),ifelse(median(duration)<5,1,0)), dist_col=round(median(max_dist),0), travel_dist=round(median(total_dist),0),
		durL=round(min(duration),ifelse(median(duration)<5,1,0)), distL=round(min(max_dist),0), travelL=round(min(total_dist),0),
		durU=round(max(duration),ifelse(median(duration)<5,1,0)), distU=round(max(max_dist),0), travelU=round(max(total_dist),0)) %>%
	mutate(range=paste(dist_col," (",distL,"-",distU,")",sep="")) %>%
	mutate(travel_distance=paste(travel_dist," (",travelL,"-",travelU,")", sep="")) %>%
	mutate(trip_duration=paste(median_duration," (",durL,"-",durU,")", sep="")) %>%
	mutate(trip_duration=paste(median_duration," (",durL,"-",durU,")", sep="")) %>%
	arrange(Family, scientific_name,desc(breed_stage)) %>%
	select(Family, scientific_name,breed_stage, trip_duration,travel_distance, range)

TABLE2add<- ORIG %>%
  	group_by(Family, scientific_name, breed_stage) %>%
  	summarise(MMAm=round(median(MMA)/1000,ifelse(median(MMA)<1000,2,0)), MMAl=round(min(MMA)/1000,ifelse(median(MMA)<1000,2,0)), MMAu=round(max(MMA)/1000,ifelse(median(MMA)<1000,2,0)), BAm=round(median(BA),2),BAl=round(min(BA),2),BAu=round(max(BA),2)) %>%
	mutate(AREA=if_else(MMAl==MMAu,as.character(MMAm),paste(MMAm," (",MMAl,"-",MMAu,")",sep=""))) %>%
	mutate(BAIndex=if_else(BAl==BAu,as.character(BAm),paste(BAm," (",BAl,"-",BAu,")",sep=""))) %>%
	arrange(Family, scientific_name,desc(breed_stage)) %>%
	select(Family, scientific_name,breed_stage, AREA,BAIndex)

TABLE2 <- merge(TABLE2,TABLE2add, by =c("Family", "scientific_name","breed_stage"))

TABLE2<- TABLE2 %>% arrange(Family, scientific_name, desc(breed_stage))
TABLE2$breed_stage <-ifelse(TABLE2$breed_stage=="breeding","unspec. breeding",TABLE2$breed_stage)



#fwrite(TABLE2,"S:\\ConSci\\DptShare\\SteffenOppel\\MANUSCRIPTS\\in_prep\\SeabirdPrioritisation\\Table2_chick.csv")

ORIG[ORIG$scientific_name=="Hydroprogne caspia",]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE NEW FIGURE 1 RANGE BOXPLOT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### CREATE SORT ORDER OF FAMILIES
sortfam<-TRIPS %>%
  filter(breed_stage %in% c('chick-rearing')) %>%
  group_by(Family) %>%
  summarise(range=median(max_dist))%>%
  arrange(desc(range)) %>%
  select(Family)

getwd()

#pdf("Fig1_chick.pdf", height=8, width=8)
#jpeg("S:\\ConSci\\DptShare\\SteffenOppel\\MANUSCRIPTS\\in_prep\\SeabirdPrioritisation\\Fig1_chick.jpg", height=780, width=680, quality=100)

TRIPS %>%
  filter(breed_stage %in% c('incubation','chick-rearing')) %>%
  mutate(breed_stage=as.factor(breed_stage)) %>%
  mutate(breed_stage=factor(breed_stage, levels=levels(breed_stage)[c(2,1)])) %>%
  mutate(Family=as.factor(Family)) %>%
  mutate(Family=fct_relevel(Family,sortfam$Family)) %>%
  
  ggplot(aes(x=Family, y=log10(max_dist), width=1), size=1)+geom_boxplot()+
  facet_wrap("breed_stage", ncol=1, scales = "fixed")+
  xlab("Seabird family") +
  ylab(expression(Maximum ~ distance ~from ~colony ~(10^x ~ km))) +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text.y=element_text(size=18, color="black"),
        axis.text.x=element_text(size=14, color="black", angle=90, vjust=0.5, hjust=0.99), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

dev.off()






#### USE ONLY FIRST TRIP PER BIRD ####
head(TRIPS)

FIRSTTRIPS<- TRIPS %>%  group_by(Family, scientific_name, ID) %>%
	summarise(firsttrip=min(trip))

TABLES1<- TRIPS %>%
	filter(trip %in% FIRSTTRIPS$firsttrip) %>%
  	group_by(Family, scientific_name, breed_stage) %>%
  	summarise(median_duration=round(median(duration),ifelse(median(duration)<5,1,0)), dist_col=round(median(max_dist),0), travel_dist=round(median(total_dist),0),
		durL=round(min(duration),ifelse(median(duration)<5,1,0)), distL=round(min(max_dist),0), travelL=round(min(total_dist),0),
		durU=round(max(duration),ifelse(median(duration)<5,1,0)), distU=round(max(max_dist),0), travelU=round(max(total_dist),0)) %>%
	mutate(range=paste(dist_col," (",distL,"-",distU,")",sep="")) %>%
	mutate(travel_distance=paste(travel_dist," (",travelL,"-",travelU,")", sep="")) %>%
	mutate(trip_duration=paste(median_duration," (",durL,"-",durU,")", sep="")) %>%
	mutate(trip_duration=paste(median_duration," (",durL,"-",durU,")", sep="")) %>%
	select(Family, scientific_name, breed_stage,travel_distance, range,trip_duration)

TABLES1 <- merge(TABLES1,TABLE2add, by =c("Family", "scientific_name", "breed_stage"))
TABLES1<- TABLES1 %>% arrange(Family, scientific_name, desc(breed_stage))
TABLES1$breed_stage <-ifelse(TABLES1$breed_stage=="breeding","unspec. breeding",TABLES1$breed_stage)

#fwrite(TABLES1,"S:\\ConSci\\DptShare\\SteffenOppel\\MANUSCRIPTS\\in_prep\\SeabirdPrioritisation\\TableS2_chick.csv")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# STATISTICAL ANALYSIS OF FAMILY EFFECT FOR MAXIMUM DISTANCE TO COLONY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

head(TRIPS)

TRIPS<-TRIPS %>% 
  mutate(COL_SIZE_sc=as.numeric(scale(COL_SIZE))) %>%
  mutate(SamplingRate_sc=as.numeric(scale(SamplingRate)))

### FOR ORIGINAL DATA ###
m0range<-lmer(log(max_dist)~breed_stage+COL_SIZE_sc+SamplingRate_sc+(1|colony_name)+(1|scientific_name), data=TRIPS)
m1range<-lmer(log(max_dist)~Family+breed_stage+COL_SIZE_sc+SamplingRate_sc+(1|colony_name)+(1|scientific_name), data=TRIPS)
mORIGaov<-anova(m0range, m1range)
modsum<-summary(m1range)
modsum$coefficients



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# STATISTICAL ANALYSIS OF FAMILY EFFECT FOR BA OVERLAP INDEX
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

head(ORIG %>% arrange(BA))
head(ORIG %>% arrange(desc(BA)))

### FOR ORIGINAL DATA ###
ORIG<-ORIG %>% mutate(n_ind_sc=as.numeric(scale(n_individuals))) %>%
     mutate(COL_SIZE_sc=as.numeric(scale(COL_SIZE))) %>%
     mutate(SamplingRate_sc=as.numeric(scale(SamplingRate)))

m0BA<-glmer(BA~breed_stage+COL_SIZE_sc+SamplingRate_sc+(1|colony_name)+(1|scientific_name), data=ORIG, family=binomial, weights=n_individuals)
m1BA<-glmer(BA~breed_stage+Family+COL_SIZE_sc+SamplingRate_sc+(1|colony_name)+(1|scientific_name), data=ORIG, family=binomial, weights=n_individuals)

mBAaov<-anova(m0BA, m1BA)
modsum<-summary(m1BA)
out<-as.data.frame(modsum$coefficients)
out$parameter<-rownames(out)

out %>% arrange(desc(Estimate))





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE FIGURE 2 BA BOXPLOT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summary(ORIG$BA)

ORIG %>% arrange(BA) %>% select(scientific_name,breed_stage,BA,Scale)

### CREATE SORT ORDER OF FAMILIES
sortfam<-ORIG %>%
  filter(breed_stage %in% c('chick-rearing')) %>%
  group_by(Family) %>%
  summarise(range=median(BA))%>%
  arrange(desc(range)) %>%
  select(Family)


#pdf("Fig2_chick.pdf", height=8, width=8)
#jpeg("S:\\ConSci\\DptShare\\SteffenOppel\\MANUSCRIPTS\\in_prep\\SeabirdPrioritisation\\Fig1_chick.jpg", height=780, width=680, quality=100)

ORIG %>%
  filter(breed_stage %in% c('incubation','chick-rearing')) %>%
  mutate(breed_stage=as.factor(breed_stage)) %>%
  mutate(breed_stage=factor(breed_stage, levels=levels(breed_stage)[c(2,1)])) %>%
  mutate(Family=as.factor(Family)) %>%
  mutate(Family=fct_relevel(Family,sortfam$Family)) %>%

ggplot(aes(x=Family, y=BA, width=1), size=1)+geom_boxplot()+
  facet_wrap("breed_stage", ncol=1, scales = "fixed")+
  xlab("Seabird family") +
  ylab("Bhattacharyya's affinity index") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text.y=element_text(size=18, color="black"),
        axis.text.x=element_text(size=14, color="black", angle=90, vjust=0.5, hjust=0.99), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

dev.off()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ASSESSMENT OF REPRESENTATIVITY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
head(SIMUL)
ORIG$representativity<-0
ORIG$IncVal<-0

dgs<-unique(SIMUL$DataGroup)


for (d in dgs){

x<-SIMUL[SIMUL$DataGroup==d,]

try(M1 <-nls((x$InclusionMean ~ (a*x$SampleSize)/(1+b*x$SampleSize)), data=x, start=list(a=1,b=0.1)), silent = TRUE)
if ('M1' %in% ls()){       ### run this only if nls was successful

out<- x%>% mutate(pred=predict(M1)) %>%
	group_by(SampleSize) %>%
	summarise(mean=mean(pred), sd=sd(InclusionMean), actual=mean(InclusionMean))

Asymptote <- (summary(M1)$coefficients[1]/summary(M1)$coefficients[2])
ORIG$representativity[ORIG$DataGroup==d] <- max(out$mean)/Asymptote

### plot for testing ###
#plot(InclusionMean~SampleSize, data=x, pch=16, cex=0.2, col="darkgray", ylim=c(0,1), xlim=c(0,max(x$SampleSize)), ylab="Inclusion", xlab="SampleSize")
#  yTemp <- c(as.numeric((out$mean + 0.5*out$sd)), as.numeric(rev(out$mean  - 0.5*out$sd)))
#  xTemp <- c(1:nrow(out), nrow(out):1)
#  polygon(x=xTemp, y=yTemp, col="gray93", border=F)
#  points(InclusionMean~SampleSize, data=x, pch=16, cex=0.2, col="darkgray")
#  lines(out, lty=1,lwd=2)

}else{ORIG$representativity[ORIG$DataGroup==d] <- mean(x$InclusionMean[x$SampleSize==max(x$SampleSize)])}   ### if nls is unsuccessful then use mean output for largest sample size
ORIG$IncVal[ORIG$DataGroup==d]<-mean(x$InclusionMean[x$SampleSize==max(x$SampleSize)])
rm(M1,x,out,Asymptote)
}

ORIG$REPRESENT<-apply(ORIG[,30:31],1,max)			### take the max of either asymptote or inclusion value



### Troubleshoot negative representativity

ORIG %>% arrange(representativity)

ORIG$representativity<-ifelse(ORIG$representativity<0,0,ORIG$representativity)


### SUMMARISE PROPORTION OF REPRESENTATIVE DATASETS
### CHECK WHY THIS ADDS TO 189 and not 186


hist(ORIG$REPRESENT, breaks=c(0,0.7,0.8,0.9,1.0), plot=F)$counts/length(ORIG$representativity)

REPRES<-ORIG %>%
	mutate(REP=ifelse(representativity<0.7,ifelse(IncVal<0.5,"not_rep","rep20"),ifelse(representativity<0.8,"rep20",ifelse(representativity<0.9,"rep12.5","rep10"))))%>%
	mutate(count=1) %>%
	group_by(Family,REP) %>%
	summarise(N=sum(count)) %>%
	spread(REP,N, fill=0) %>%
	mutate(Tot_N=sum(not_rep,rep10,rep12.5,rep20)) %>%
	mutate(not_rep=round(not_rep/Tot_N,3),rep10=rep10/Tot_N,rep12.5=rep12.5/Tot_N,rep20=rep20/Tot_N)

REPREStot<-ORIG %>%
  mutate(REP=ifelse(representativity<0.7,ifelse(IncVal<0.5,"not_rep","rep20"),ifelse(representativity<0.8,"rep20",ifelse(representativity<0.9,"rep12.5","rep10"))))%>%
  mutate(count=1) %>%
  group_by(REP) %>%
  summarise(N=sum(count))
REPREStot$N[REPREStot$REP=="not_rep"]
	
ORIG %>%
	mutate(REP=ifelse(representativity<0.7,"not rep",ifelse(representativity<0.8,"rep20",ifelse(representativity<0.9,"rep12.5","rep10"))))%>%
	select(DataGroup, Family, scientific_name,colony_name,n_individuals,IncVal,REP)%>%
	arrange(Family, IncVal)


#fwrite(REPRES,"S:\\ConSci\\DptShare\\SteffenOppel\\MANUSCRIPTS\\in_prep\\SeabirdPrioritisation\\TableS4.csv")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOT FIGURE S1: REPRESENTATIVITY AGAINST SAMPLE SIZE 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

colourpalette<-c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999','#000120')



#pdf("FigS1_chick.pdf", height=6, width=8)
#jpeg("S:\\ConSci\\DptShare\\SteffenOppel\\MANUSCRIPTS\\in_prep\\SeabirdPrioritisation\\FigS1_chick.jpg", height=520, width=680, quality=100)

ORIG %>%
  mutate(BreedingStage=ifelse(chick==1,"chick-rearing","other")) %>%
  filter(Family!="") %>%
  
  ggplot(aes(x=n_individuals, y=REPRESENT,colour=Family, pch=BreedingStage))+
  geom_point(size=2.5) +
  scale_colour_manual(values = colourpalette)+
  geom_hline(aes(yintercept=0.7),colour="red", size=0.5, linetype=2) +
  xlab("number of individuals tracked") +
  ylab("Representativity") +
  labs(shape="Breeding Stage", colour="Family")+
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),
        legend.title=element_text(size=20),
        legend.text=element_text(size=16),
        legend.key = element_rect(fill = "white"),
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

dev.off()


REPCOR<-cor.test(ORIG$REPRESENT,ORIG$n_individuals, method="pearson")
REPCOR$p.value


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# STATISTICAL ANALYSIS OF FAMILY EFFECT FOR IBA SIZE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
head(ORIG)
ORIG$weight<-apply(ORIG[,27:28],1,max)
ORIG$weight<-ifelse(ORIG$weight==0,0.001,ORIG$weight)

### FOR ORIGINAL DATA ###
m0orig<-lmer(log(IBA20)~breed_stage+COL_SIZE+n_individuals+SamplingRate+(1|colony_name)+(1|scientific_name), data=ORIG[ORIG$IBA20>0,], weights=ORIG$weight[ORIG$IBA20>0])
m1orig<-lmer(log(IBA20)~breed_stage+COL_SIZE+n_individuals+Family+SamplingRate+(1|colony_name)+(1|scientific_name), data=ORIG[ORIG$IBA20>0,], weights=ORIG$weight[ORIG$IBA20>0])
mIBAaov<-anova(m0orig, m1orig)
modsum<-summary(m1orig)
out<-as.data.frame(modsum$coefficients)
out$parameter<-rownames(out)

out %>% arrange(desc(Estimate))

str(mIBAaov[2,8])





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# STATISTICAL ANALYSIS OF FAMILY EFFECT FOR MMA SIZE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
head(ORIG)

### FOR ORIGINAL DATA ###
m0orig<-lmer(log(MMA)~breed_stage+COL_SIZE+n_individuals+SamplingRate+(1|colony_name)+(1|scientific_name), data=ORIG[ORIG$IBA20>0,], weights=ORIG$weight[ORIG$IBA20>0])
m1orig<-lmer(log(MMA)~breed_stage+COL_SIZE+n_individuals+Family+SamplingRate+(1|colony_name)+(1|scientific_name), data=ORIG[ORIG$IBA20>0,], weights=ORIG$weight[ORIG$IBA20>0])
mMMAaov<-anova(m0orig, m1orig)
modsum<-summary(m1orig)
MMAout<-as.data.frame(modsum$coefficients)
MMAout$parameter<-rownames(out)

MMAout %>% arrange(desc(Estimate))

str(mMMAaov[2,8])



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATING FIGURE 3: HISTOGRAM OF IBA AREA SIZES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## changed to NOT mention IBA as concept not introduced
## included two shades for stages

head(SIMUL)

#pdf("Fig3.pdf", height=8, width=6)
#jpeg("S:\\ConSci\\DptShare\\SteffenOppel\\MANUSCRIPTS\\in_prep\\SeabirdPrioritisation\\Fig2.jpg", height=680, width=520, quality=100)

SIMUL %>%
  filter(SampleSize>5) %>%
  
  ggplot()+
  geom_histogram(aes(x=log10(IBA20),y=..count../tapply(..count..,..PANEL..,sum)[..PANEL..]),binwidth=0.2)+    
  #geom_histogram(data=subset(SIMUL,chick==1), aes(x=log10(IBA20),y=..count../tapply(..count..,..PANEL..,sum)[..PANEL..]),binwidth=0.2,fill = "darkred", alpha = 0.5)+    
  #geom_histogram(data=subset(SIMUL,chick==0), aes(x=log10(IBA20),y=..count../tapply(..count..,..PANEL..,sum)[..PANEL..]),binwidth=0.2,fill = "darkgrey", alpha = 0.5)+    
  facet_wrap("Family", ncol=2, scales = "fixed")+
  ylab("Proportion of simulated populations") +
  xlab(expression(Area ~ of ~concentrated ~use ~(10^x ~ km^{2}))) +  
  scale_x_continuous(limits=c(0,5))+
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text.y=element_text(size=18, color="black"),
        axis.text.x=element_text(size=18, color="black", vjust=0.5), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

dev.off()






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATING FIGURE 4: SCATTERPLOT BA AGAINST MMA SIZE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(RColorBrewer)
colourpalette<-c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999','#000120')

#pdf("Fig4_chick.pdf", height=6, width=8)
#jpeg("S:\\ConSci\\DptShare\\SteffenOppel\\MANUSCRIPTS\\in_prep\\SeabirdPrioritisation\\Fig4_chick.jpg", height=520, width=680, quality=100)

ORIGplot<-ORIG %>%
  mutate(BreedingStage=ifelse(chick==1,"chick-rearing","other")) %>%
  filter(Family!="")
  
  ggplot()+
  geom_histogram(data=WDPA, aes(x=log10(GIS_M_AREA),y=..count../tapply(..count..,..PANEL..,sum)[..PANEL..]),binwidth=0.4,fill='lightgrey')+ 
  geom_point(data=ORIGplot, aes(x=log10(MMA), y=BA,colour=Family, pch=BreedingStage),size=2.5) +
  labs(shape="Breeding Stage", colour="Family")+
  scale_colour_manual(values = colourpalette)+
  xlab(expression(Size ~ of ~marine ~area ~(10^x ~ km^{2}))) +  
  scale_x_continuous(limits=c(0,8))+
  ylab("Bhattacharrya's Affinity index") +
  
  #geom_segment(aes(x=log10(large5MPA),y=0,xend=log10(large5MPA), yend=1),colour="red", size=0.5, linetype=2) +

  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),
        legend.title=element_text(size=20),
        legend.text=element_text(size=16),
        legend.key = element_rect(fill = "white"),
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CORRELATION OF BA AGAINST MMA SIZE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ORIG %>%
  mutate(BreedingStage=ifelse(chick==1,"chick-rearing","other")) %>%
  filter(Family!="")

## 
BAMMACOR<-cor.test(ORIG$BA,ORIG$MMA, method="pearson")
BAMMACOR$p.value
round(BAMMACOR$estimate,2)
paste(round(BAMMACOR$conf.int,2), sep =" - ", collapse=" - ")






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SENSITIVITY OF BA TO H VALUE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


setwd("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Outputv15")
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\Outputv15")
BA_SENS<-fread("BA_sensitivity.csv")
BA_SENS<-BA_SENS[!(BA_SENS$DataGroup==32),]
BA_SENS$breed_stage<-ifelse(BA_SENS$breed_stage %in% c("brood-guard","post-guard"),"chick-rearing",BA_SENS$breed_stage)
head(BA_SENS)

## CORRELATION TEST
cor.test(BA_SENS$BA,BA_SENS$BAalt1, method="pearson")


## TABLE FOR SUPPLEMENTARY MATERIAL
TABLES4<- BA_SENS %>%
  group_by(Family, scientific_name, breed_stage) %>%
  summarise(BAm=round(median(BA),2),BAl=round(min(BA),2),BAu=round(max(BA),2),BAm1=round(median(BAalt1),2),BAl1=round(min(BAalt1),2),BAu1=round(max(BAalt1),2)) %>%
  mutate(BAIndex=if_else(BAl==BAu,as.character(BAm),paste(BAm," (",BAl,"-",BAu,")",sep=""))) %>%
  mutate(BAIndexALT=if_else(BAl1==BAu1,as.character(BAm1),paste(BAm1," (",BAl1,"-",BAu1,")",sep=""))) %>%
  arrange(Family, scientific_name, desc(breed_stage)) %>%
  select(Family, scientific_name, breed_stage, BAIndex, BAIndexALT) 

#fwrite(TABLES4,"S:\\ConSci\\DptShare\\SteffenOppel\\MANUSCRIPTS\\in_prep\\SeabirdPrioritisation\\TableS3_chick.csv")


  BA_SENS %>%

  ggplot()+
geom_point(aes(x=BA, y=BAalt1,colour=Family),size=2.5) +
  labs(colour="Family")+
  scale_colour_manual(values = colourpalette)+
  xlab("h from area restricted search") +  
  ylab("h from foraging range") +
theme(panel.background=element_rect(fill="white", colour="black"), 
      axis.text=element_text(size=18, color="black"), 
      axis.title=element_text(size=20),
      legend.title=element_text(size=20),
      legend.text=element_text(size=16),
      legend.key = element_rect(fill = "white"),
      strip.text.x=element_text(size=18, color="black"), 
      strip.background=element_rect(fill="white", colour="black"), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_blank())





##################################################################
### PRODUCE RESULTS SECTION WITH KEY TABLES AND FIGURES ###
##################################################################
#detach(packages:htmlwidgets)
#detach(name="package:htmlwidgets", unload=TRUE, character.only=TRUE)
#install.packages(c('plotly','htmlwidgets'), dependencies=T)

library(markdown)
library(rmarkdown)
library(knitr)

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\MANUSCRIPTS\\in_prep\\SeabirdPrioritisation")


### create HTML report for overall summary report
Sys.setenv(RSTUDIO_PANDOC="C:/Program Files (x86)/RStudio/bin/pandoc")

rmarkdown::render('S:\\ConSci\\DptShare\\SteffenOppel\\MANUSCRIPTS\\in_prep\\SeabirdPrioritisation\\Results_SeabirdPrioritisation_v9.Rmd',
                  output_file = "Results_v9.docx",
                  output_dir = 'S:\\ConSci\\DptShare\\SteffenOppel\\MANUSCRIPTS\\in_prep\\SeabirdPrioritisation')



Sys.setenv(RSTUDIO_PANDOC="C:\\Program Files\\RStudio\\bin\\pandoc")
rmarkdown::render('C:\\STEFFEN\\MANUSCRIPTS\\in_prep\\SeabirdPrioritisation\\Results_SeabirdPrioritisation_v8.Rmd',
                  output_file = "Results_v8.docx",
                  output_dir = 'C:\\STEFFEN\\MANUSCRIPTS\\in_prep\\SeabirdPrioritisation')

head(TABLE2)


##################################################################
### CREATE AUTHOR LIST AND EMAIL ADDRESS ###
##################################################################

setwd("C:\\STEFFEN\\MANUSCRIPTS\\in_prep\\SeabirdPrioritisation")
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\MANUSCRIPTS\\in_prep\\SeabirdPrioritisation")
authors<-fread("Author_details.csv")

paste(authors$Coauthor, sep=", ", collapse=", ")
paste(authors$Coauthor, sep="; ", collapse="; ")



paste(authors$email, sep="; ", collapse="; ")

nonrespondents <- authors %>% filter(RESPOND!="YES")
paste(nonrespondents$email, sep="; ", collapse="; ")