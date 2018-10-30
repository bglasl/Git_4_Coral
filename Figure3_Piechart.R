##### Piechart Figure_3 ####
# Piechart

load(file="PHYLOSEQ_TABLE.RData")
library(phyloseq)

PHYLOSEQ_TABLE_control=subset_samples(PHYLOSEQ_TABLE, Treatment!="field_control")
PHYLOSEQ_TABLE_control

#### PIECHART ####
OTU<-as.data.frame(otu_table(PHYLOSEQ_TABLE_control))
library(data.table)
setDT(OTU, keep.rownames = TRUE)[]
OTU<-dplyr::rename(OTU, Feature_ID=rn)

TAXA<-as.data.frame(tax_table(PHYLOSEQ_TABLE_control))

library(dplyr)
library(tidyr)
Piecharttable<-right_join(OTU, TAXA, by="Feature_ID")

Endos_all<-Piecharttable %>% filter(Genus=="D_5__Endozoicomonas") %>%
  tidyr::gather("Sample","Abundance", 2:97) %>% 
  group_by(Sample) %>% 
  summarise(MEAN_Endo=sum(Abundance)) %>% 
  summarise(MEAN_Endo_total=mean(MEAN_Endo), SD=sd(MEAN_Endo))
Endos_all

Others_all<-Piecharttable %>% filter(Genus!="D_5__Endozoicomonas") %>%
  tidyr::gather("Sample","Abundance", 2:97) %>% 
  group_by(Sample) %>% 
  summarise(MEAN_Endo=sum(Abundance)) %>% 
  summarise(MEAN_Endo_total=mean(MEAN_Endo), SD=sd(MEAN_Endo))
Others_all

Endos<-Piecharttable %>% filter(Genus=="D_5__Endozoicomonas") %>%
  filter(Feature_ID!="c0e6c66cfd0673747e59936d6c050407") %>% 
  tidyr::gather("Sample","Abundance", 2:97) %>% 
  group_by(Sample) %>% 
  summarise(MEAN_Endo=sum(Abundance)) %>% 
  summarise(MEAN_Endo_total=mean(MEAN_Endo), SD=sd(MEAN_Endo))
Endos

Core<-Piecharttable %>%
  filter(Feature_ID=="c0e6c66cfd0673747e59936d6c050407") %>% 
  tidyr::gather("Sample","Abundance", 2:97) %>% 
  group_by(Sample) %>% 
  summarise(MEAN_Endo=sum(Abundance)) %>% 
  summarise(MEAN_Endo_total=mean(MEAN_Endo),SD=sd(MEAN_Endo))

#mean Endo abundance (minus core) = 28.98157
#mean Core = 18.91282
#mean Others =(100-(28.98157+18.91282)
Endo<-data_frame(28.98157,18.91282,(100-(28.98157+18.91282)))
NEW<-t(Endo)
NEW<-as.data.frame(NEW)
setDT(NEW, keep.rownames = TRUE)

library(ggplot2)
plot<-ggplot(NEW, aes(x = 1, weight = V1, fill = rn)) +
  geom_bar(width = 1, position = "stack") +
  coord_polar(theta = "y")+
  theme_minimal(base_size = 16, base_family = "Helvetica")+
  scale_fill_manual(values=c("#aaaaaa", "#a6cee3","#1f78b4"), 
                    name="",labels=c("others","core","Endozoicomonas"))+
  theme(legend.position="bottom")+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank())

plot
pdf('Figure_2a', width=6, height=6)
print(plot)
graphics.off()


#### Figure 2b-c ####
# all samples - general trends

setwd("~/Documents/2_Coral_Microbiome/Stats_2/Coral_Microbiome_analysis")
load(file="PHYLOSEQ_TABLE.RData")
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)


PHYLOSEQ_TABLE_control=subset_samples(PHYLOSEQ_TABLE, Treatment!="field_control")
PHYLOSEQ_TABLE_control


#### NMDS ####
set.seed(200)
ORDCCA<-ordinate(PHYLOSEQ_TABLE_control,"NMDS","bray")
ORDCCA # stress 0.163385

p0 = plot_ordination(Endo_Control_PHY_rel, ORDCCA, color="Genotype", shape="Treatment")+
  scale_color_brewer(palette = "Paired")+
  theme_classic(base_size = 16, base_family = "Helvetica")+
  geom_point(size=3)+
  facet_wrap(~Genotype)
p0

postscript(file = "NMDS_all.eps", width =10, height = 10)
print(p0)
dev.off()

pdf('NMDS_all.pdf', width=10, height=10)
print(p0)
graphics.off()

#### ADONIS & BETA DISPERSION ####
# Adonis = permutational Multivariate Analysis of variance using distance matrices
### using adonis2 ### - this is the way to go!
df=as(sample_data(PHYLOSEQ_TABLE_control), 'data.frame')
d=phyloseq::distance(PHYLOSEQ_TABLE_control,'bray')

ADONIS<-adonis(d~Genotype, data=df, permutations = 10000, method = "bray")
ADONIS #sign -> sign -  community varies significanlty between coral host genotypes
capture.output(ADONIS,file="Adonis_allGenotype.doc")

ADONIS<-adonis(d~Treatment, data=df, permutations = 10000, method = "bray")
ADONIS #sign -> sign -  community varies significanlty between coral host genotypes
capture.output(ADONIS,file="Adonis_allGenotype.doc")


perm<-how(nperm=10000)
setBlocks(perm)<-with(df, Genotype)
ADONIS2<-adonis2(d~Treatment+SamplingTimepoint, data=df, permutations = perm, method = "bray")
ADONIS2 #not sign -> no shift in the microbial community between treatments within a genotype
capture.output(ADONIS2,file="Adonis2_allGenotypeblocked.doc")

# graph Betadisper ####
library(ggplot2)
Genotype_group<-get_variable(PHYLOSEQ_TABLE_control, "Genotype")
BETADISP<-betadisper(d, Genotype_group,type=c("centroid"))
anova(BETADISP) # varies between Genotype significant <2.2e-16

Disp<-as.data.frame(BETADISP$distances) 
Sample_ID<-rownames(Disp)
Disp<-cbind(Sample_ID, Disp)
Disp_new<-right_join(METADATA, Disp, by="Sample_ID")
Disp_new<-Disp_new%>% group_by(Genotype) %>% arrange(desc(BETADISP$distances), .by_group=TRUE)
Disp_new$Treatment<-factor(Disp_new$Treatment, levels=c("control", "stress_ambient", "stress_2100"))

save(Disp_new, file="Disp_new.RData")

plot<-ggplot(Disp_new, aes(y=(BETADISP$distances),x=SamplingTimepoint))+
  geom_point(aes(color=Genotype, shape=Treatment), size=3)+
  geom_line(aes(group=Treatment, linetype=Treatment, color=Genotype))+
  scale_color_brewer(palette = "Paired")+
  theme_classic(base_size = 16, base_family = "Helvetica")+
  facet_wrap(~Genotype)+
  guides(color=guide_legend(title="Genotype"))+
  labs(x= "day", y = "distance to group centroids")+
  scale_x_discrete(breaks=c("t1","t3","t4","t5"), labels=c("1","10","14","19"))

plot

pdf('Betadisper_all.pdf', width=10, height=10)
print(plot)
graphics.off()

setEPS()
postscript('Betadisper_all.eps', width=10, height=10)
plot
dev.off()


df=as(sample_data(PHYLOSEQ_TABLE_control), 'data.frame')
df=df%>%dplyr::mutate(Genotype_Treatment=paste(Genotype,Treatment,sep="_"))
d=phyloseq::distance(PHYLOSEQ_TABLE_control,'bray')


Host_Treatment_group<-get_variable(df, "Genotype_Treatment")
BETADISP<-betadisper(d, Host_Treatment_group,type=c("centroid"))
anova(BETADISP) #significant 0.006481 **
boxplot(BETADISP)
(BETADISP)
permutest(BETADISP, pairwise = TRUE, permutations = 10000)
BETADISP_HSD<-TukeyHSD(BETADISP)
BETADISP_HSD


Treatment_group<-get_variable(PHYLOSEQ_TABLE_control, "Treatment")
BETADISP<-betadisper(d, Treatment_group,type=c("centroid"))
anova(BETADISP)

#Mantel test ####
# transfrom metadata z-score
TEMP<-read.csv("sample_metadata_new2.csv", header = TRUE, sep=",", dec=".", strip.white = TRUE)
str(TEMP)
TEMP<-TEMP %>% filter(Treatment != "field_control")
TEMP
TEMP_stand<-decostand(TEMP[,c(9,10,12,13)], method = "normalize", na.rm=TRUE)
TEMP_stand

# upload z-score metadata
Sample_ID<-TEMP$Sample_ID
newMETA<-cbind(TEMP_stand, Sample_ID)
str(newMETA)
rownames(newMETA)<-newMETA$Sample_ID

META<-inner_join(TEMP, newMETA, by="Sample_ID")
rownames(META)<-META$Sample_ID
METADATA<-sample_data(META)

sample_data(PHYLOSEQ_TABLE_control)<-METADATA
sample_data(PHYLOSEQ_TABLE_control)

library(dplyr)
conMETA<-sample_data(PHYLOSEQ_TABLE_control) %>% select(PSYield.y, Chla.y, Zoox.y, Protein.y)

bray.dist.otu<-phyloseq::distance(PHYLOSEQ_TABLE_control,'bray') #using rarefied species data

bray.dist.host<-vegdist(conMETA, method = "bray") #using scaled host health proxy data

bray.mantel<-mantel(bray.dist.otu, bray.dist.host, permutations = 10000, method="pearson")

bray.mantel

# dbRDA ####
ORDCCA<-ordinate(PHYLOSEQ_TABLE_control,"CAP",formula = ~Genotype+Treatment+ Zoox+Protein+Chla+PSYield)
ORDCCA # explains 32 % of observed variation

library(vegan)
ANOVA<-anova.cca(ORDCCA, by="term", permutations = 1000) # 
ANOVA # Genotype sig
capture.output(ANOVA,file="ANOVA_CCA_GenotypeHostHealth.doc")


