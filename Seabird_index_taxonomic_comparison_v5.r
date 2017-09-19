############################################################################################################
####### SEABIRD PRIORITISATION BASED ON SPATIAL AGGREGATION INDEX ANALYSIS #################################
############################################################################################################
## this analysis is based on Spatial_index_COMPLETE_REANALYSIS.r
## output compiled by Seabird_index_OUTPUT_COMPILATION.r
## developed by steffen.oppel@rspb.org.uk in December 2016

## v4 modified 23 August 2017
## based on complete new reanalysis of cleaned data from 23 Aug 2017


## v5modified on 7 Sept 2017 after inclusion of storm petrel data
## removed the less interesting figures, exploratory plots comparing lots of indices (see v1-4 for code)
## removed latitudinal variation plot (most families have very little latitudinal variation)

## modified on 13 Sept 2017 to include the output of the representativity analysis

## modified on 15 Sept 2017 after re-run with new ARS settings
## included analysis of 'clumpedness' to detect datasets where IBA area is larger than MMA (which should not be possible)

## modified 18 Sept 2019 to include AREA output

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD PACKAGES AND CUSTOM SCRIPTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)
require(foreign)
library(readr)
library(data.table)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD SAVED RESULTS FROM DATA GROUP SPECIFIC ANALYSIS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs")
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs")
load("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs\\AggIndex_compiled_output_v12.RData")

### use breed stage as brood-guard or other
ORIG$chick<-ifelse(ORIG$breed_stage=="brood-guard",1,0)
SIMUL$chick<-ifelse(SIMUL$breed_stage=="brood-guard",1,0)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SIMPLE SUMMARIES FOR TEXT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
length(unique(ORIG$DataGroup))  ### number of datasets analysed
length(unique(ORIG$scientific_name))  ### number of species used
dim(TRIPS)
length(unique(TRIPS$ID))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATING SUMMARY OUTPUT TABLE OF USED DATASETS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#detach(package:plyr)             ## if the functions below do not work then plyr is the problem
head(SIMUL)
summarytab<- SIMUL %>%
  group_by(Family, breed_stage) %>%
  summarise(n_datasets=length(unique(DataGroup))) %>%
  spread(key=breed_stage, value=n_datasets, fill=0)

samplesizetab<- overview %>%
  filter(DataGroup %in% unique(SIMUL$DataGroup)) %>%
  group_by(Family) %>%
  summarise(average=round(mean(n_individuals),0),min=min(n_individuals), max=max(n_individuals)) %>%
  select(Family, average, min, max)

tabtripSummary<- TRIPS %>%
  group_by(Family, scientific_name) %>%
  summarise(mean_duration=round(mean(duration),1), dist_col=round(mean(max_dist),0), travel_dist=round(mean(total_dist),0))


#### INTERROGATE OUTPUT THAT LOOKS VERY WRONG:

TRIPS %>% filter(total_dist>10000)


#### ASSESS NUMBER OF FORAGING TRIPS PER DATA GROUP

TRIPS %>%
  mutate(count=1) %>%
  group_by(DataGroup) %>%
  summarise(n_trips=sum(count)) %>%
  arrange(desc(n_trips))


#### ASSESS AVERAGE NUMBER OF IND PER DATASET

overview %>%
  filter(DataGroup %in% unique(SIMUL$DataGroup)) %>%
  summarise(average=round(mean(n_individuals),0))

round((dim(overview[overview$n_individuals<10])[1]/dim(overview)[1])*100,1)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE AND PLOT CLUMPEDNESS RATIO IN ORIGINAL DATA SET
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
head(ORIG)
library(RColorBrewer)

colourpalette<-c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999','#000120')




ORIG %>% 
  mutate(Clumpedness=log(IBA20)/log(MMA)) %>%
  #mutate(Scale=Scale*50) %>%
  #filter(Clumpedness>1) %>%
  filter(Family!="") %>%
  #select(DataGroup,Family,scientific_name,colony_name,breed_stage,MMA,IBA20)

  
#pdf("IBA_MMA_ratio.pdf", height=7, width=9)

  ggplot(aes(x=Family, y=Clumpedness, width=1), size=1)+geom_boxplot()+
  xlab("Seabird family") +
  ylab("IBA/MMA ratio") +
  geom_hline(aes(yintercept=1),colour="red", size=1) +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text.y=element_text(size=18, color="black"),
        axis.text.x=element_text(size=14, color="black", angle=45, vjust=0.5), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

#dev.off()






ORIG %>%
  filter(IBA10>0) %>%
  mutate(Clumpedness=IBA10/MMA) %>%
  mutate(BreedingStage=ifelse(chick==1,"brood-guard","other")) %>%
  filter(Family!="") %>%
  
  ggplot(aes(x=log(MMA), y=log(Clumpedness),colour=Family, pch=BreedingStage))+
  geom_point(size=2.5) +
  #scale_colour_brewer(type = "qual" , palette = "Set1") +
  scale_colour_manual(values = colourpalette)+
  xlab("log(size) of marine management area") +
  ylab("log(ratio) of area with concentrated activity") +
  geom_segment(aes(x=log(MPAsizes[5]),y=-4,xend=log(MPAsizes[5]), yend=1),colour="red", size=0.5, linetype=2) +
  #geom_segment(aes(x=0,y=0.5,xend=log(MPAsizes[5]), yend=0.5),colour="red", size=0.5, linetype=2) +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOT MMA SIZE AGAINST OVERLAP RATIO IN ORIGINAL DATA SET SO THAT GRADIENT IS VISIBLE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


pdf("Seabird_prioritisation_BA_IBA scatter.pdf", height=7, width=9)

ORIG %>%
  mutate(BreedingStage=ifelse(chick==1,"brood-guard","other")) %>%
  filter(Family!="") %>%
  
  ggplot(aes(x=log(MMA), y=BA,colour=Family, pch=BreedingStage))+
  geom_point(size=2.5) +
  #scale_colour_brewer(type = "qual" , palette = "Set1") +
  scale_colour_manual(values = colourpalette)+
  xlab("log(size) of marine management area") +
  ylab("Bhattacharrya's Affinity index") +
  geom_segment(aes(x=log(MPAsizes[5]),y=0,xend=log(MPAsizes[5]), yend=1),colour="red", size=0.5, linetype=2) +
  geom_segment(aes(x=0,y=0.5,xend=log(MPAsizes[5]), yend=0.5),colour="red", size=0.5, linetype=2) +
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
# CREATING HISTOGRAMS OF THE SIZE OF MANAGEMENT AREAS FOR EACH FAMILY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



head(SIMUL)

plotout<- SIMUL %>%
  #filter(InclusionMean>0.4) %>%
  filter(SampleSize>5)


pdf("Seabird_prioritisation_MMA_histogram.pdf", height=10, width=7)

#jpeg("Seabird_prioritisation_MMA_histogram.jpg", height=580, width=640, quality=100)
ggplot(plotout)+    ## colour=breed_stage looks shit
  geom_histogram(aes(x=log(MMA),y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]),binwidth=0.5)+                               ##breaks=log(MPAsizes)
  #geom_smooth(fill="lightblue", size=1.5)+
  facet_wrap("Family", ncol=2, scales = "fixed")+
  #geom_point(colour="black", size=1.5) +
  geom_vline(aes(xintercept=log(MPAsizes[5])),colour="red", size=0.5, linetype=2) +
  geom_vline(aes(xintercept=log(MPAsizes[4])),colour="red", size=0.5, linetype=2) +
  #scale_y_continuous(name="Proportion of tracked populations",labels = percent, limits = c(0,1))+
  ylab("Proportion of tracked populations") +
  xlab("log(Size) of marine management area (sq km)") +
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





pdf("Seabird_prioritisation_IBA20_histogram.pdf", height=10, width=7)
ggplot(plotout)+    ## colour=breed_stage looks shit
  geom_histogram(aes(x=log(IBA20),y=..count../tapply(..count..,..PANEL..,sum)[..PANEL..]),binwidth=0.5)+                               ##breaks=log(MPAsizes)
  #geom_smooth(fill="lightblue", size=1.5)+
  facet_wrap("Family", ncol=2, scales = "fixed")+
  #geom_point(colour="black", size=1.5) +
  geom_vline(aes(xintercept=log(MPAsizes[5])),colour="red", size=0.5, linetype=2) +
  geom_vline(aes(xintercept=log(MPAsizes[4])),colour="red", size=0.5, linetype=2) +
  ylab("Proportion of tracked populations") +
  xlab("log(Size) of marine IBA (sq km)") +
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







pdf("Seabird_prioritisation_BA_histogram.pdf", height=10, width=7)

ggplot(plotout)+ 
  geom_histogram(aes(x=BA,y=..count../tapply(..count..,..PANEL..,sum)[..PANEL..]),binwidth=0.05)+                               ##breaks=log(MPAsizes)
  facet_wrap("Family", ncol=2, scales = "fixed")+
  ylab("Proportion of tracked populations") +
  xlab("Bhattacharyya's affinity index") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text.y=element_text(size=18, color="black"),
        axis.text.x=element_text(size=15, color="black"), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())


dev.off()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TEST AND QUANTIFY THE FAMILY EFFECT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### use breed stage as brood-guard or other
ORIG$chick<-ifelse(ORIG$breed_stage=="brood-guard",1,0)
SIMUL$chick<-ifelse(SIMUL$breed_stage=="brood-guard",1,0)


library(lme4)

### FOR ORIGINAL DATA ###
m0orig<-lmer(log(IBA20)~chick+n_individuals+(1|colony_name)+(1|scientific_name), data=ORIG[ORIG$IBA20>0,])
m1orig<-lmer(log(IBA20)~chick+n_individuals+Family+(1|colony_name)+(1|scientific_name), data=ORIG[ORIG$IBA20>0,])
anova(m0orig, m1orig)
modsum<-summary(m1orig)
#modsum$coefficients


### FOR SIMULATED DATA ###
m0simul<-lmer(log(IBA20)~chick+n_individuals+SampleSize+(1|colony_name)+(1|scientific_name), data=SIMUL[SIMUL$IBA20>0,])
m1simul<-lmer(log(IBA20)~chick+n_individuals+SampleSize+Family+(1|colony_name)+(1|scientific_name), data=SIMUL[SIMUL$IBA20>0,])
anova(m0simul, m1simul)
modsum<-summary(m1simul)
#modsum$coefficients




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATING PLOTS OF SAMPLE SIZE AND REPRESENTATIVITY FOR EACH FAMILY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



pdf("Seabird_prioritisation_InclusionValue_incubation.pdf", height=10, width=7)

ggplot(SIMUL[SIMUL$breed_stage=="incubation",], aes(x=SampleSize, y=InclusionMean), size=0.1)+
  geom_point(colour="black", size=0.1) +
  geom_smooth(fill="lightblue", size=1.5, method='nls', formula = y  ~ (a*x)/(1+b*x), se=FALSE)+
  facet_wrap("Family", ncol=2, scales = "fixed")+
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




pdf("Seabird_prioritisation_InclusionValue_brood_guard.pdf", height=10, width=7)

ggplot(SIMUL[SIMUL$breed_stage=="brood-guard",], aes(x=SampleSize, y=InclusionMean), size=0.1)+
  geom_point(colour="black", size=0.1) +
  geom_smooth(fill="lightblue", size=1.5, method='nls', formula = y  ~ (a*x)/(1+b*x), se=FALSE)+
  facet_wrap("Family", ncol=2, scales = "fixed")+
  xlab("N tracked individuals") +
  xlim(0,100)+
  ylab("prop. of untracked locations in 50%UD") +
  #scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20))+
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())
dev.off()






pdf("Seabird_prioritisation_MMA_breed_stage.pdf", height=10, width=7)

ggplot(plotout, aes(x=breed_stage, y=log(MMA), width=1), size=1)+geom_boxplot()+
#geom_smooth(fill="lightblue", size=1.5)+
facet_wrap("Family", ncol=2, scales = "fixed")+
#geom_point(colour="black", size=1.5) +
geom_hline(aes(yintercept=log(MPAsizes[5])),colour="red", size=0.5) +
  xlab("Breeding Stage") +
  ylab("log(Size) of marine management area (sq km)") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text.y=element_text(size=18, color="black"),
        axis.text.x=element_text(size=14, color="black", angle=45, vjust=0.5), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())
dev.off()







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOT AREA INCREASE FOR PROPORTION OF INDIVIDUALS PROTECTED
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

head(AREA)
AREA %>% mutate(prop.ind=step/n_tracks) %>%
  mutate(MMA=ORIG$MMA[match(AREA$DataGroup,ORIG$DataGroup)]) %>%

ggplot(aes(x=prop.ind, y=log(tot.area/1000000)), size=0.1)+
  geom_point(colour="black", size=0.1) +
  geom_smooth(fill="lightblue", size=1.5, method='nls', formula = y  ~ (a*x)/(1+b*x), se=FALSE)+
  #geom_smooth(fill="lightblue", size=1.5, method='lm', se=T)+
  facet_wrap("Family", ncol=2, scales = "fixed")+
  xlab("Proportion of tracked individuals") +
  ylab("cumulative log(area) of core foraging ranges") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())









#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# USE OUTPUT FROM THE ACTUAL DATA SET AND CONTRAST VARIOUS METRICS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

summary <- ORIG %>%
  group_by(Family) %>%
  summarise(BAmean=mean(BA),BAsd=sd(BA),BAmax=max(BA),BAmin=min(BA),MMAmean=mean(MMA),MMAsd=sd(MMA),MMAmax=max(MMA),MMAmin=min(MMA)) %>%
  arrange(BAsd)
summary




pdf("Seabird_prioritisation_BA_range.pdf", height=7, width=9)

ORIG %>%
  filter(breed_stage %in% c('incubation','brood-guard')) %>%

ggplot(aes(x=Family, y=BA, width=1), size=1)+geom_boxplot()+
  #geom_smooth(fill="lightblue", size=1.5)+
  facet_wrap("breed_stage", ncol=1, scales = "fixed")+
  #geom_point(colour="black", size=1.5) +
  geom_hline(aes(yintercept=0),colour="red", size=0.5) +
  xlab("Seabird family") +
  ylab("Bhattacharyya's affinity index") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text.y=element_text(size=18, color="black"),
        axis.text.x=element_text(size=14, color="black", angle=45, vjust=0.5), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

dev.off()













##################################################################
### PRODUCE OUTPUT REPORT WITH KEY TABLES AND FIGURES ###
##################################################################
#detach(packages:htmlwidgets)
#detach(name="package:htmlwidgets", unload=TRUE, character.only=TRUE)
#install.packages(c('plotly','htmlwidgets'), dependencies=T)

library(markdown)
library(rmarkdown)
library(knitr)
library(plotly)

plotout<- SIMUL %>%
  #filter(InclusionMean>0.4) %>%
  filter(SampleSize>5)

plotout2<- SIMUL %>%
  #filter(InclusionMean>0.4) %>%
  filter(breed_stage %in% c('incubation','brood-guard')) %>%
  filter(SampleSize>5)


  ### create HTML report for overall summary report
Sys.setenv(RSTUDIO_PANDOC="C:/Program Files (x86)/RStudio/bin/pandoc")

  rmarkdown::render('S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Analysis\\SummarySeabirdPrioritisation_v3_internal_policy.Rmd',
                    output_file = "SeabirdPrioritisation_Progress_Report.html",
                    output_dir = 'S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs')
  
  
  rmarkdown::render('S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Analysis\\SummarySeabirdPrioritisation_v3_externalScientists.Rmd',
                    output_file = "SeabirdSpaceUse_Progress_Report.html",
                    output_dir = 'S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs')
  
  
  
  ##################################################################
  ### TROUBLESHOOT MYSTERIOUS OUTPUT ###
  ##################################################################
  
ORIG %>%
  filter(Family=="Procellariidae") %>%
  filter(breed_stage=="brood-guard") %>%
  select(DataGroup, scientific_name, colony_name, n_individuals,MMA)

  
  
SIMUL %>%
  filter(Family=="Procellariidae") %>%
  filter(breed_stage=="brood-guard") %>%
  select(DataGroup, scientific_name, colony_name, SampleSize,InclusionMean)


  
  

