##### Figure 6 ####
# IndVal 

load(file="PHYLOSEQ_TABLE.RData")
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(microbiome)


PHYLOSEQ_TABLE=subset_samples(PHYLOSEQ_TABLE, Treatment!="field_control")
PHYLOSEQ_TABLE
Other_PHY<-subset_samples(PHYLOSEQ_TABLE, SamplingTimepoint!="t1") 
Other_PHY_rel = transform_sample_counts(Other_PHY, function(x)100*x/sum(x))

METADATA<-as(sample_data(Other_PHY_rel), "data.frame")
METADATA<-rownames_to_column(METADATA, var="Sample_ID")

OTU_table<-as(otu_table(Other_PHY_rel), "matrix")
if(taxa_are_rows(Other_PHY_rel)){OTU_table <- t(OTU_table)}
OTU_table_df = as.data.frame(OTU_table)
OTU_table_df<-rownames_to_column(OTU_table_df, var="Sample_ID")

# join datasets
library(dplyr)
library(tidyr)
OTU_table_IndVal<-right_join(METADATA, OTU_table_df, by="Sample_ID")

#indval
library(labdsv)
library(MASS)
library(vegan)
library(cluster)
library(indicspecies)
library(permute)
(INDVAL_OTUs_species_GC=(as.data.frame(OTU_table_IndVal[,14:4637])))

#Cluster 
(INDVAL_Groups_Origin_GC=(as.character(OTU_table_IndVal$Treatment)))
INDVAL_Origin_GC=multipatt(INDVAL_OTUs_species_GC, INDVAL_Groups_Origin_GC, func="IndVal.g", duleg=FALSE, control=how(nperm=10000))
summary(INDVAL_Origin_GC, indvalcomp=TRUE)
head(summary.multipatt(INDVAL_Origin_GC, alpha = 0.05, indvalcomp=TRUE))
IndValdf<-data.frame(INDVAL_Origin_GC$sign)
IndValdf<-rownames_to_column(IndValdf, var="ASV_ID")
IndValdf<-IndValdf %>% filter(p.value<=0.05) %>% dplyr:: rename(control=s.control) %>% 
  dplyr:: rename(`cumulative stress`=s.stress_2100)%>%
  dplyr:: rename(`acute stress`=s.stress_ambient)

IndValdf<-IndValdf%>% dplyr::select(-index,-p.value,-stat)

library(reshape)
IndVallong<-melt(IndValdf)
indicatorNUMBER<-IndVallong %>% group_by(variable,value) %>% summarise(numberofOTUs=n())

# Component ‘A’ is the probability that the surveyed
# site belongs to the target site group given the fact that the species has been found. 
# This conditional probability is called the specificity or positive predictive value of the species as indicator of the site group. 
# Component ‘B’ is the probability of finding the species in sites belonging to the site group. 
# This second conditional probability is called the fidelity or sensitivity of the species as indicator of the target site group
capture.output(summary.multipatt(INDVAL_Origin_GC, alpha = 0.05,indvalcomp=TRUE),file="IndVal_Treatment.csv")

## prep Graphs ####
# Define the taxa you want like this:
goodTaxa = (IndValdf$ASV_ID)
goodTaxa
allTaxa = taxa_names(Other_PHY_rel)
allTaxa <- allTaxa[(allTaxa %in% goodTaxa)]
ex1 = prune_taxa(allTaxa, Other_PHY_rel)
# new phyloseq object with just the taxa you kept.
ex1 #12 taxa

OTUTABLE<-as(otu_table(ex1), "matrix")
OTUTABLE<-as.data.frame(OTUTABLE)
OTUTABLE<-rownames_to_column(OTUTABLE, var="ASV_ID")

test<-inner_join(IndVallong,OTUTABLE, by=c("ASV_ID"))
test1<-test%>%gather("Sample_ID","Abundance",4:74)

# graph ####
TAXA<-as(tax_table(ex1),"matrix")
TAXA<-data.frame(TAXA)
TAXA<-rownames_to_column(TAXA, var="ASV_ID")


all1<-test1 %>%filter(value==1) %>% inner_join(TAXA,test1, by=c("ASV_ID")) %>% inner_join(METADATA, test1, by=c("Sample_ID")) 

all3<-all1 %>% group_by(variable,ASV_ID) %>% summarise(MEAN=mean(Abundance))
all3
all4<-inner_join(TAXA,all3, by=c("ASV_ID"))
all4
all4$variable<-factor(all4$variable, levels=c("control", "acute stress", "cumulative stress"))


p<-ggplot(all4, aes(x=Family, y=(MEAN), color=Genus)) + geom_point(size=6) + 
  facet_grid(~variable)+
  theme_bw(base_size = 16, base_family = "Helvetica")+
  coord_flip()+
  ylab("mean rel.abundance [%]")
p

pdf('IndVal_Treatment_TAXONOMY.pdf', width=10, height=8)
print(p)
graphics.off()

postscript(file = "IndVal_Treatment_TAXONOMY.eps", width =12, height = 6)
print(p)
dev.off()
