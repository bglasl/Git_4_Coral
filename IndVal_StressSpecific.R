##### Figure 6 ####
# IndVal others

load(file="PHYLOSEQ_TABLE.RData")
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(microbiome)


PHYLOSEQ_TABLE=subset_samples(PHYLOSEQ_TABLE, Treatment!="field_control")
PHYLOSEQ_TABLE
Other_PHY<-subset_taxa(PHYLOSEQ_TABLE, Genus!="D_5__Endozoicomonas")
Other_PHY #96 samples
Other_PHY_rel = transform_sample_counts(Other_PHY, function(x)100*x/sum(x))

METADATA<-as(sample_data(Other_PHY_rel), "data.frame")
METADATA<-rownames_to_column(METADATA, var="Sample_ID")

OTU_table<-as(otu_table(Other_PHY_rel), "matrix")
if(taxa_are_rows(Other_PHY_rel)){OTU_table <- t(OTU_table)}
OTU_table_df = as.data.frame(OTU_table)
OTU_table_df<-rownames_to_column(OTU_table_df, var="Sample_ID")

# join datasets
TEM<-read.csv("Grouping_zeropointfive_threshold.csv", header = TRUE, sep = ",", dec=".", strip.white=TRUE)
METAnew<-full_join(TEM,METADATA,by="Sampling_Date")
METAnew
library(dplyr)
library(tidyr)
OTU_table_IndVal<-right_join(METAnew, OTU_table_df, by="Sample_ID")

#indval

library(labdsv)
library(MASS)
library(vegan)
library(cluster)
library(indicspecies)
library(permute)
(INDVAL_OTUs_species_GC=(as.data.frame(OTU_table_IndVal[,15:5003])))

#Cluster 
(INDVAL_Groups_Origin_GC=(as.character(OTU_table_IndVal$avg_temp)))
INDVAL_Origin_GC=multipatt(INDVAL_OTUs_species_GC, INDVAL_Groups_Origin_GC, func="IndVal.g", duleg=TRUE, control=how(nperm=1000))
summary(INDVAL_Origin_GC, indvalcomp=TRUE)
head(summary.multipatt(INDVAL_Origin_GC, alpha = 0.01, indvalcomp=TRUE)) #At=0.80, Bt=0.80 
IndValdf<-data.frame(INDVAL_Origin_GC$sign)
IndValdf<-rownames_to_column(IndValdf, var="ASV_ID")
IndValdf<-IndValdf %>% filter(p.value<=0.01) %>% dplyr:: rename(low=s.low) %>% 
  dplyr:: rename(average=s.average)%>%
  dplyr:: rename(high=s.high)

IndValdf<-IndValdf%>% dplyr::select(-index,-p.value,-stat)

library(reshape)
IndVallong<-melt(IndValdf)
indicatorNUMBER<-IndVallong %>% group_by(variable,value) %>% summarise(numberofOTUs=n())

# Component ‘A’ is the probability that the surveyed
# site belongs to the target site group given the fact that the species has been found. 
# This conditional probability is called the specificity or positive predictive value of the species as indicator of the site group. 
# Component ‘B’ is the probability of finding the species in sites belonging to the site group. 
# This second conditional probability is called the fidelity or sensitivity of the species as indicator of the target site group
# capture.output(summary.multipatt(INDVAL_Origin_GC, alpha = 0.05,At=0.90, Bt=0.90,indvalcomp=TRUE),file="IndVal_SeawaterTemperature.csv")

## prep Graphs ####
# Define the taxa you want like this:
goodTaxa = (IndValdf$ASV_ID)
goodTaxa
allTaxa = taxa_names(PL_MI)
allTaxa <- allTaxa[(allTaxa %in% goodTaxa)]
ex1 = prune_taxa(allTaxa, PL_MI)
# new phyloseq object with just the taxa you kept.
ex1 #110taxa

OTUTABLE<-as(otu_table(ex1), "matrix")
OTUTABLE<-as.data.frame(OTUTABLE)
OTUTABLE<-rownames_to_column(OTUTABLE, var="ASV_ID")

test<-inner_join(IndVallong,OTUTABLE, by=c("ASV_ID"))
test1<-test%>%gather("Sample","Abundance",4:33)
test2<-test1%>% group_by(variable,value,Sample) %>% summarise(SUM=sum(Abundance))

META<-METAnew %>% dplyr::select(Sampling_Date,Sample_ID) %>% dplyr::rename(Sample=Sample_ID)
test3<-inner_join(META,test2, by=c("Sample"))

library(raster)
MeanOTUtable<-test3 %>% group_by(variable,value,Sampling_Date) %>% summarise(MEAN=mean(SUM),CV=cv(SUM),SD=sd(SUM)) %>%filter(value==1)

TEMnew<-TEM%>%dplyr::select(Sampling_Date,avg_temp) 
test4<-inner_join(TEMnew,MeanOTUtable, by=c("Sampling_Date")) %>% rowwise(.) %>% mutate(Indicator = ifelse(avg_temp %in% variable, "yes", "no"))

test4$variable = factor(test4$variable, levels=c('low','average','high'))

p<-ggplot(test4, aes(x=Sampling_Date, y=MEAN)) +
  geom_point(aes(size=CV, color=Indicator))+
  scale_size_area(max_size = 12, breaks=c(5,10,20))+
  scale_color_manual(values=c("grey","black"))+
  facet_grid(~variable)+
  theme_bw(base_size = 16, base_family = "Helvetica")+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  ylab("relative abundance [%]")+
  xlab("sampling date")
p

pdf('IndVal_Temp_seawater_new.pdf', width=15, height=10)
print(p)
graphics.off()

postscript(file = "IndVal_Temp_seawater_mew.eps", width =15, height = 10)
print(p)
dev.off()

test4$avg_temp = factor(test4$avg_temp, levels=c("high","average", "low"))

p<-ggplot(test4, aes(x=Sampling_Date, y=""))+
  theme_bw(base_size = 16, base_family = "Helvetica")+
  geom_raster(aes(fill=avg_temp))+
  scale_fill_manual(values=c("#b2182b","white","#2166ac"))
p
postscript(file = "IndVal_Temp_seawater_gradientlegend.eps", width =11, height = 1)
print(p)
dev.off()

# alluvial graph ####
TAXA<-as(tax_table(ex1),"matrix")
TAXA<-data.frame(TAXA)

all1<-test1 %>%filter(value==1) %>% inner_join(TAXA,test1, by=c("ASV_ID")) %>% inner_join(META,test1, by=c("Sample")) 

all2<-inner_join(TEMnew,all1, by=c("Sampling_Date")) %>% rowwise(.) %>% mutate(Indicator = ifelse(avg_temp %in% variable, "yes", "no")) %>% filter(Indicator=="yes")
all3<-all2 %>% group_by(variable,ASV_ID) %>% summarise(MEAN=mean(Abundance))
all3
all4<-all3 %>% filter(variable=="high"|variable=="low") %>% filter(MEAN>0)
all4<-inner_join(TAXA,all4, by=c("ASV_ID"))
all4
all4a<-all3 %>% filter(MEAN>0)
all4a<-inner_join(TAXA,all4a, by=c("ASV_ID"))
write.csv(all4a, "Temp_alluvial.csv")

p<-ggplot(all4, aes(x=genus, y=sqrt(MEAN), color=phylum)) + geom_point(size=6) + 
  facet_wrap(~variable , ncol = 1)+
  scale_color_brewer(palette = "Set2")+
  theme_bw(base_size = 16, base_family = "Helvetica")+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  scale_color_manual(values=c("#b0b7c1","#d64242","#dda42a","#9dc663","#44baba","#7b4ec4"))+
  ylab("sqrt(rel.abundance)")
p

pdf('IndVal_Temp_seawater_TAXONOMY.pdf', width=8, height=10)
print(p)
graphics.off()

postscript(file = "IndVal_Temp_seawater_TAXONOMY.eps", width =8, height = 12)
print(p)
dev.off()