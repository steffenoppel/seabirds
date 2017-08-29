############################################################################################################
####### SEABIRD PRIORITISATION BASED ON SPATIAL AGGREGATION INDEX ANALYSIS #################################
############################################################################################################
## this analysis is based on Seabird_spatial_agg_index_v9.r
## developed by steffen.oppel@rspb.org.uk in December 2016
## modified on 31 July 2017

## modified 23 August 2017
## based on complete new reanalysis of cleaned data from 23 Aug 2017


### IDEAS TO DEVELOP: include latitudinal, pelagic vs continental or environmental analysis
### pick species with most datasets and look at intra-specific variation compared to variation within family or across seabirds




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

setwd("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation\\Data")
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Data")

overview<-fread("Seabird_priority_overview_v2.csv")
head(overview)


setwd("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs")
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs")


TRIPS<-fread("Seabird_trip_summaries_all.csv")
TRIPS<-merge(TRIPS,overview, by="DataGroup", all.x=T)
head(TRIPS)

SpatInd_orig<-fread("SpatIndex_OrigData_all.csv")
head(SpatInd_orig)

SpatInd_simul<-fread("SpatIndex_Simulations_all.csv")
head(SpatInd_simul)


WDPA<-fread("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\WDPA_MPA_Aug2017.csv")
head(WDPA)
MPAsizes<-quantile(WDPA$GIS_M_AREA,c(0.5,0.75,0.9,0.95,1), na.rm=T)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# STANDARDISE BREEDING STAGES ACROSS DIFFERENT GROUPS AND REMOVE NON-BREEDING STAGES FROM ALL THREE DATASETS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SpatInd_simul %>%
  group_by(Family, breed_stage) %>%
  summarise(n_tracks=max(n_tracks)) %>%
  arrange(n_tracks) %>%
  print(n=40)

dim(SpatInd_simul)
SIMUL<-SpatInd_simul %>%
  filter(!(breed_stage=="pre-egg")) %>%
  filter(!(breed_stage=="pre-moult")) %>%
  mutate(breed_stage=as.character(breed_stage)) %>%
  mutate(breed_stage=ifelse(breed_stage=="chick-rearing","brood-guard",breed_stage)) %>%
  mutate(breed_stage=ifelse(breed_stage=="unknown","breeding",breed_stage)) %>%
  mutate(breed_stage=ifelse(breed_stage=="creche","post-guard",breed_stage))%>%
  select(DataGroup,Family,scientific_name,colony_name,breed_stage,Nind,Ntrips, BA,MMA,Scale,extent,E_W,N_S,IBA10,IBA125,IBA20,MCP10,MCP125,MCP20,Morisita,mclu,muni,imst)
dim(SIMUL)


dim(SpatInd_orig)
ORIG<-SpatInd_orig %>%
  filter(!(breed_stage=="pre-egg")) %>%
  filter(!(breed_stage=="pre-moult")) %>%
  mutate(breed_stage=as.character(breed_stage)) %>%
  mutate(breed_stage=ifelse(breed_stage=="chick-rearing","brood-guard",breed_stage)) %>%
  mutate(breed_stage=ifelse(breed_stage=="unknown","breeding",breed_stage)) %>%
  mutate(breed_stage=ifelse(breed_stage=="creche","post-guard",breed_stage)) %>%
  select(DataGroup,Family,scientific_name,colony_name,breed_stage,Nind,Ntrips, BA,MMA,Scale,extent,E_W,N_S,IBA10,IBA125,IBA20,MCP10,MCP125,MCP20,Morisita,mclu,muni,imst)
dim(ORIG)


TRIPS2<-TRIPS %>%
  filter(!(breed_stage=="pre-egg")) %>%
  filter(!(breed_stage=="pre-moult")) %>%
  mutate(breed_stage=as.character(breed_stage)) %>%
  mutate(breed_stage=ifelse(breed_stage=="chick-rearing","brood-guard",breed_stage)) %>%
  mutate(breed_stage=ifelse(breed_stage=="unknown","breeding",breed_stage)) %>%
  mutate(breed_stage=ifelse(breed_stage=="creche","post-guard",breed_stage)) %>%
  select(DataGroup,ID,trip,Family,scientific_name,colony_name,breed_stage,device,duration,max_dist,total_dist)
dim(TRIPS2)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATING SUMMARY OUTPUT TABLE OF USED DATASETS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
length(unique(SIMUL$DataGroup))
length(unique(SIMUL$scientific_name))
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

tabtripSummary<- TRIPS2 %>%
  group_by(Family, scientific_name) %>%
  summarise(mean_duration=round(mean(duration),1), dist_col=round(mean(max_dist),0), travel_dist=round(mean(total_dist),0))


#### INTERROGATE OUTPUT THAT LOOKS VERY WRONG:

TRIPS2 %>% filter(total_dist>10000)


#### ASSESS NUMBER OF FORAGING TRIPS PER DATA GROUP

TRIPS2 %>%
  mutate(count=1) %>%
  group_by(DataGroup) %>%
  summarise(n_trips=sum(count)) %>%
  arrange(desc(n_trips))






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATING HISTOGRAMS OF THE SIZE OF MANAGEMENT AREAS FOR EACH FAMILY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("C:\\STEFFEN\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs")
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs")


head(SIMUL)

plotout<- SIMUL %>%
  #filter(InclusionMean>0.4) %>%
  filter(SampleSize>5)


pdf("Seabird_prioritisation_MMA_histogram.pdf", height=7, width=9)

#jpeg("Seabird_prioritisation_MMA_histogram.jpg", height=580, width=640, quality=100)
ggplot(plotout)+    ## colour=breed_stage looks shit
  geom_histogram(aes(x=log(MMA),y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]),binwidth=0.5)+                               ##breaks=log(MPAsizes)
  #geom_smooth(fill="lightblue", size=1.5)+
  facet_wrap("Family", ncol=3, scales = "fixed")+
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





pdf("Seabird_prioritisation_IBA20_histogram.pdf", height=7, width=9)
ggplot(plotout)+    ## colour=breed_stage looks shit
  geom_histogram(aes(x=log(IBA20),y=..count../tapply(..count..,..PANEL..,sum)[..PANEL..]),binwidth=0.5)+                               ##breaks=log(MPAsizes)
  #geom_smooth(fill="lightblue", size=1.5)+
  facet_wrap("Family", ncol=3, scales = "fixed")+
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







pdf("Seabird_prioritisation_BA_histogram.pdf", height=7, width=9)

ggplot(plotout)+ 
  geom_histogram(aes(x=BA,y=..count../tapply(..count..,..PANEL..,sum)[..PANEL..]),binwidth=0.05)+                               ##breaks=log(MPAsizes)
  facet_wrap("Family", ncol=3, scales = "fixed")+
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


########## CHECK WHY FRIGATEBIRDS ARE SO LOW ######

SIMUL %>%
  select(DataGroup, scientific_name, breed_stage, site_name,BA,MMA, InclusionMean,SampleSize) %>%
  filter(DataGroup %in% c(57,58,60)) %>%    ## removes all non-friagetbird datasets
  filter(InclusionMean>0.4) %>%             ## removes the unrepresentative datasets
  filter(SampleSize>10)                     ## removes the Tobago data





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATING PLOTS OF SAMPLE SIZE AND REPRESENTATIVITY FOR EACH FAMILY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



pdf("Seabird_prioritisation_InclusionValue_incubation.pdf", height=7, width=9)


ggplot(SIMUL[SIMUL$breed_stage=="incubation",], aes(x=SampleSize, y=InclusionMean), size=0.1)+
  geom_point(colour="black", size=0.1) +
  geom_smooth(fill="lightblue", size=1.5, method='nls', formula = y  ~ (a*x)/(1+b*x), start=list(a=0.5, b=0), se=FALSE)+
  facet_wrap("Family", ncol=3, scales = "fixed")+
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




pdf("Seabird_prioritisation_InclusionValue_brood_guard.pdf", height=7, width=9)


ggplot(SIMUL[SIMUL$breed_stage=="brood-guard",], aes(x=SampleSize, y=InclusionMean), size=0.1)+
  geom_point(colour="black", size=0.1) +
  geom_smooth(fill="lightblue", size=1.5, method='nls', formula = y  ~ (a*x)/(1+b*x), start=list(a=0.5, b=0), se=FALSE)+
  facet_wrap("Family", ncol=3, scales = "fixed")+
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






pdf("Seabird_prioritisation_MMA_breed_stage.pdf", height=7, width=9)

ggplot(plotout, aes(x=breed_stage, y=log(MMA), width=1), size=1)+geom_boxplot()+
#geom_smooth(fill="lightblue", size=1.5)+
facet_wrap("Family", ncol=3, scales = "fixed")+
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
# USE OUTPUT FROM THE ACTUAL DATA SET AND CONTRAST VARIOUS METRICS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUMMARY<-fread("SeabirdPriority_orig_output_fixed_samplesize.csv")
head(SUMMARY)

summary <- SUMMARY %>%
  group_by(Family) %>%
  summarise(BAmean=mean(BA),BAsd=sd(BA),BAmax=max(BA),BAmin=min(BA),MMAmean=mean(MMA),MMAsd=sd(MMA),MMAmax=max(MMA),MMAmin=min(MMA)) %>%
  arrange(BAsd)
summary




pdf("Seabird_prioritisation_BA_range.pdf", height=7, width=9)

SUMMARY %>%
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








pdf("Seabird_prioritisation_Morisita.pdf", height=7, width=9)

SUMMARY %>%
  filter(breed_stage %in% c('incubation','brood-guard')) %>%
  
  ggplot(aes(x=Family, y=imst, width=1), size=1)+geom_boxplot()+
  #geom_smooth(fill="lightblue", size=1.5)+
  facet_wrap("breed_stage", ncol=1, scales = "fixed")+
  #geom_point(colour="black", size=1.5) +
  geom_hline(aes(yintercept=0.5),colour="red", size=0.5) +
  xlab("Seabird family") +
  ylab("Morisita's aggregation index") +
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





pdf("Seabird_indices.pdf", height=7, width=9)

SUMMARY %>%
  filter(breed_stage %in% c('incubation','brood-guard')) %>%
  select(DataGroup, Family, BA, MMA, IBA10, MCP10, Morisita, mclu, imst, Scale,muni) %>%
  #group_by(DataGroup, Family) %>%
  gather(key=index, value=value,BA, MMA, IBA10, MCP10, Morisita, mclu, imst, Scale,muni) %>%
  
  ggplot(aes(x=Family, y=value, width=1), size=1)+geom_boxplot()+
  #geom_smooth(fill="lightblue", size=1.5)+
  facet_wrap("index", ncol=3, scales = "free_y")+
  #geom_point(colour="black", size=1.5) +
  #geom_hline(aes(yintercept=0.5),colour="red", size=0.5) +
  xlab("Seabird family") +
  ylab("Value of spatial indices") +
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

  rmarkdown::render('A:\\RSPB\\Marine\\SeabirdPrioritisation\\Analysis\\SummarySeabirdPrioritisation_v2.Rmd',
                    output_file = "SeabirdPrioritisation_Summary_Report.html",
                    output_dir = 'A:\\RSPB\\Marine\\SeabirdPrioritisation\\Outputs')
  

