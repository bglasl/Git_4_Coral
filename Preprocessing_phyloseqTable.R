#############################
#### Coral Experiment ######
#############################

library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

### load data ####
# load OTU table
OTU_table<-read.csv('otu_table.csv', header = TRUE, sep = ",", dec = ".", strip.white = TRUE)
head(OTU_table)
str(OTU_table)

OTU_table_T<-read.csv('otu_table3506_T.csv', header = TRUE, sep = ",", dec = ".", strip.white = TRUE)
head(OTU_table_T) # OTUs in rows
str(OTU_table_T)

METADATA<-read.csv("sample_metadata_new2.csv", header = TRUE, sep = ",")
head(METADATA)
levels(METADATA$Genotype)

TAXADATA<-read.csv("taxonomy.csv", header = TRUE, sep = ",")
head(TAXADATA)

# rooted tree
TREE<-read_tree("tree.nwk")
TREE<-phy_tree(TREE)

# join datasets
library(dplyr)
library(tidyr)
OTU_table_new<-right_join(METADATA, OTU_table, by="Sample_ID")
levels(OTU_table_new$Treatment)

### prep Phyloseq Table ####
# Taxa table
TAXA_TABLE<-as.data.frame(TAXADATA[,1])
rownames(TAXADATA)<-TAXADATA[,1]
TAXA_TABLE_split<-TAXADATA %>% separate(Taxon, c("Domain","Phylum","Class","Order","Family","Genus","Species"), sep=";", remove=FALSE) 
TAXA_TABLE_matrix<-as.matrix(TAXA_TABLE_split)
TAXA_TABLE<-tax_table(TAXA_TABLE_matrix) # tax_table=matrix

# OTU table
rownames(OTU_table_T)<-OTU_table_T[,1]
OTU_table_T<-data.frame(OTU_table_T[,-1])

OTU_TABLE<-otu_table(OTU_table_T, taxa_are_rows = TRUE)
sample_names(OTU_TABLE)

# Sample data
SAMPLEDATA<-OTU_table_new[,1:13]
rownames(SAMPLEDATA)<-OTU_table_new$Sample_ID
SAMPLEDATA<-data.frame(SAMPLEDATA[,-1])
levels(SAMPLEDATA$Genotype)

SAMPLEDATA<-sample_data(SAMPLEDATA) #Sample_data = data.frame
sample_names(SAMPLEDATA)
SAMPLEDATA$Treatment<-factor(SAMPLEDATA$Treatment, levels=c("field_control", "control", "stress_ambient", "stress_2100"))
SAMPLEDATA$SamplingTimepoint<-factor(SAMPLEDATA$SamplingTimepoint, levels=c("t0","t1","t3","t4","t5"))
SAMPLEDATA$Treatment_Timepoint<-factor(SAMPLEDATA$Treatment_Timepoint, levels=c("field control t0",
                                                                                "control t1","control t3", "control t4","control t5",
                                                                                "ambient stress t1","ambient stress t3", "ambient stress t4","ambient stress t5",
                                                                                "2100 stress t1","2100 stress t3", "2100 stress t4","2100 stress t5"))


#combine Sampledata, Taxadata, OTU table to phyloseq table
PHYLOSEQ_TABLE_count<-phyloseq(OTU_TABLE, SAMPLEDATA, TAXA_TABLE, TREE)
PHYLOSEQ_TABLE_count
save(PHYLOSEQ_TABLE_count, file="PHYLOSEQ_TABLE_Count.RData")

# transform sample counts to relative abundance
PHYLOSEQ_TABLE = transform_sample_counts(PHYLOSEQ_TABLE_count, function(x)100*x/sum(x))
otu_table(PHYLOSEQ_TABLE)
PHYLOSEQ_TABLE
save(PHYLOSEQ_TABLE, file="PHYLOSEQ_TABLE.RData")

sample_data(PHYLOSEQ_TABLE)
