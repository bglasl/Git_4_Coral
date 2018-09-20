#######################
#### Endozoicomonas ####

library(ggplot2)
library(phyloseq)

#load data
load(file="PHYLOSEQ_TABLE.RData")
load(file="PHYLOSEQ_TABLE_Count.RData")

#remove field_control samples
PHYLOSEQ_TABLE_control=subset_samples(PHYLOSEQ_TABLE, Treatment!="field_control")
PHYLOSEQ_TABLE_control
PHYLOSEQ_TABLE_control_count=subset_samples(PHYLOSEQ_TABLE_count, Treatment!="field_control")
PHYLOSEQ_TABLE_control_count

#only keep Endozoicmonas 
Endo_Control_PHY_rel<-subset_taxa(PHYLOSEQ_TABLE_control, Genus=="D_5__Endozoicomonas")
Endo_Control_PHY_rel
Endo_Control_PHY_count<-subset_taxa(PHYLOSEQ_TABLE_control_count, Genus=="D_5__Endozoicomonas")
Endo_Control_PHY_count

#### frequency ####
Endo_freq<-psmelt(Endo_Control_PHY_rel)
head(Endo_freq)

ANOVA<-aov(Abundance~Treatment*SamplingTimepoint, data=Endo_freq)
summary(ANOVA)
capture.output(summary(ANOVA), file="Endo_freq_TreatTime.doc")

ANOVA<-aov(Abundance~Genotype, data=Endo_freq)
summary(ANOVA)
capture.output(summary(ANOVA), file="Endo_freq_Genotype.doc")

ANOVA<-aov(Abundance~Treatment+Error(Genotype), data=Endo_freq)
summary(ANOVA)
capture.output(summary(ANOVA), file="Endo_freq_withinGenotype.doc")

#### Alpha diversity ####
Endo_Control_PHY_count #133ASVs

plot<-plot_richness(Endo_Control_PHY_count, x="Treatment", measures=c("Shannon","Observed"))+
  geom_boxplot()+
  theme_classic()
plot

pdf('ALPHA_ENDO.pdf', width=6, height=6)
print(plot)
graphics.off()

plot<-plot_richness(Endo_Control_PHY_count, x="Genotype", measures=c("Shannon","Observed"))+
  geom_boxplot()+
  theme_classic()
plot

pdf('ALPHA_ENDOGENOTYPE.pdf', width=6, height=6)
print(plot)
graphics.off()

RICHNESS<-estimate_richness(Endo_Control_PHY_count, split=TRUE, measures=c("Shannon","Observed"))
RICHNESS<-as.data.frame(RICHNESS)
RICHNESS<-tibble::rownames_to_column(as.data.frame(RICHNESS), var="Sample_ID")
Richness_table<-right_join(METADATA, RICHNESS, by="Sample_ID")

# treatment nor SamplingTimepoint has no overall effect alpha diversity (Endozoicomonas)
model1<-aov((Shannon)~SamplingTimepoint+Error(Treatment), data=Richness_table)
summary(model1) #not significant
TukeyHSD(model1)

# Effect of Genotype on the Endozoicomoas diversity
model1<-aov((Shannon)~Genotype, data=Richness_table)
summary(model1) #p=0.0452 *  significant
TukeyHSD(model1)

# Treatment also not signficant when tested wihtin a genotype
model2<-aov((Shannon)~Treatment*SamplingTimepoint+Error(Genotype), data=Richness_table)
summary(model2)
TukeyHSD(model2)

# Richness
# treatment nor SamplingTimepoint has no overall effect alpha diversity (Endozoicomonas)
model1<-aov((Observed)~SamplingTimepoint+Error(Treatment), data=Richness_table)
summary(model1) #not significant
TukeyHSD(model1)

# Effect of Genotype on the Endozoicomoas diversity
model1<-aov((Observed)~Genotype, data=Richness_table)
summary(model1) #sign p=1.14e-08 ***
TukeyHSD(model1)

# Treatment & SamplingTimepoint signficant when tested wihtin a genotype
model2<-aov((Observed)~Treatment*SamplingTimepoint+Error(Genotype), data=Richness_table)
summary(model2) #Treatment & SamplingTimepoint  significantly affect Richness within a genotype (p=0.02531 *  & p=0.00942 **, respectively)

plot_richness(CORE_count, x="Treatment", measures=c("Observed"))+
  geom_boxplot()+
  facet_grid(.~Genotype)

#### ADONIS ####
# Adonis = permutational Multivariate Analysis of variance using distance matrices
### using adonis2 ### - this is the way to go!
df=as(sample_data(Endo_Control_PHY_rel), 'data.frame')
d=phyloseq::distance(Endo_Control_PHY_rel,'bray')

ADONIS<-adonis(d~Genotype, data=df, permutations = 10000, method = "bray")
ADONIS #sign -> sign - endozoicomonas community varies significanlty between coral host genotypes
capture.output(ADONIS,file="Adonis_EndozoicomonasGenotypeControl.doc")

perm<-how(nperm=10000)
setBlocks(perm)<-with(df, Genotype)
ADONIS2<-adonis2(d~Treatment, data=df, permutations = perm, method = "bray")
ADONIS2 #not sign -> no shift in the microbial community between treatments within a genotype
capture.output(ADONIS2,file="Adonis2_TimeeffectEndoGenotypeControl.doc")

#### CCA ####
ORDCCA<-ordinate(Endo_Control_PHY_rel,"CCA", formula = ~Genotype)
plot_ordination(Endo_Control_PHY_rel, ORDCCA,color = "Genotype")
ORDCCA # Genotype explains 27.23% of observed variability
library(vegan)
anova.cca(ORDCCA,by="term", permutations = 10000) # Genotype also highly sign p=0.000999 ***

p0 = plot_ordination(Endo_Control_PHY_rel, ORDCCA, color="Genotype")+
  scale_color_brewer(palette = "Paired")+
  theme_classic(base_size = 16, base_family = "Helvetica")+
  geom_hline(yintercept=c(0,0), linetype="dotted")+
  geom_vline(xintercept=c(0,0), linetype="dotted")+
  geom_point(size=3)
p0

postscript(file = "CCAplot_EndoGenotype.eps", width =6, height = 5)
print(p0)
dev.off()

pdf('CCAplot_EndoGenotype.pdf', width=6, height=5)
print(p0)
graphics.off()

#### Endo Tree ####
Endo_Control_PHY_rel = transform_sample_counts(PHYLOSEQ_TABLE_control_count, function(x)100*x/sum(x))
head(sample_data(Endo_Control_PHY_rel))

Endo_Control_PHY_rel<-subset_taxa(PHYLOSEQ_TABLE_control, Genus=="D_5__Endozoicomonas")
test<-prune_taxa(names(sort(taxa_sums(Endo_Control_PHY_rel), TRUE))[1:12], Endo_Control_PHY_rel)
tax_table(test)
test<-subset_taxa(test, Feature_ID!="52236970d9a33a65914b7ee46274639d")
stressASV<-subset_taxa(Endo_Control_PHY_rel, Feature_ID=="ffd7161776070abe4c1c2add0317f58d"|
                         Feature_ID=="7644ba5e228f6c7e4cab13c973eebd3c"|
                         Feature_ID=="e09c9c87223e713ca3e185b78a49e825"|
                         Feature_ID=="1ca72d4dd6b29e683338f1ba895b5bc1"|
                         Feature_ID=="3964b537a5b584afc40192b91d5283b3"|
                         Feature_ID=="61a0606c0017caac996f520569863b69"|
                         Feature_ID=="b772fe52a7ceca98c424d9c16da7c826"|
                         Feature_ID=="fd123868368f46e30f2b2f7a2fa45861"|
                         Feature_ID=="52236970d9a33a65914b7ee46274639d"|
                         Feature_ID=="c0e6c66cfd0673747e59936d6c050407"|
                         Feature_ID=="8024ea926f3fdd804332b645f98fc7f4"|
                         Feature_ID=="5fb983e92b5f365c10ce1f7d494eb7e1"|
                         Feature_ID=="561f809577daa2932b264f72e2ed49fe")





stressASV
newtree<-merge_phyloseq(test, stressASV)
newtree

test = transform_sample_counts(stressASV, function(x)100*x/sum(x))
test
variable1 = as.character(get_variable(test, "Treatment"))
variable2 = as.character(get_variable(test, "Genotype"))
sample_data(test)$Tre_Gen <- mapply(paste0, variable1, variable2, 
                                    collapse = "_")

mergedGP = merge_samples(test, "Genotype")
print(mergedGP)

sample_data(mergedGP)$Genotype <- factor(sample_names(mergedGP))
sample_data(mergedGP)$Genotype <- as.character(sample_data(mergedGP)$Genotype)
sample_data(mergedGP)$Treatment <- as.character(sample_data(mergedGP)$Treatment)

str(sample_data(mergedGP))
mergedGP = transform_sample_counts(mergedGP, function(x)100*x/sum(x))

library(ggtree)
tree<-plot_tree(mergedGP, ladderize = "left", label.tips="taxa_names", color = "Genotype", size="abundance",  nodelabf = nodeplotblank, 
                base.spacing = 0.04)+
  scale_color_brewer(palette = "Paired")
  theme_tree2()+
  theme(legend.position = "right")+
  scale_size_area(max_size = 8)
tree

tree<-plot_tree(mergedGP, ladderize = "right", color = "Genotype", size="abundance",  nodelabf = nodeplotblank, 
                base.spacing = 0.05)+
  scale_color_brewer(palette = "Paired")+
  theme_tree2()+
  theme(legend.position = "right")+
  scale_size_area(max_size = 6)
tree

postscript(file = "Tree_EndoGenotype.eps", width =5, height =5 )
print(tree)
dev.off()

pdf('Tree_EndoGenotype.pdf', width=5, height=5)
print(tree)
graphics.off()


#### dbRDA ####
Endo_Control_PHY_rel = transform_sample_counts(PHYLOSEQ_TABLE_control_count, function(x)100*x/sum(x))
head(sample_data(Endo_Control_PHY_rel))

# transfrom metadata z-score
TEMP<-read.csv("sample_metadata_new2.csv", header = TRUE, sep=",", dec=".", strip.white = TRUE)
str(TEMP)
TEMP<-TEMP %>% filter(Treatment != "field_control")
TEMP
TEMP_stand<-decostand(TEMP[,c(9,10,12,13)], method = "standardize", na.rm=TRUE)
TEMP_stand

# upload z-score metadata
Sample_ID<-TEMP$Sample_ID
newMETA<-cbind(TEMP_stand, Sample_ID)
str(newMETA)
rownames(newMETA)<-newMETA$Sample_ID

META<-inner_join(TEMP, newMETA, by="Sample_ID")
rownames(META)<-META$Sample_ID
METADATA<-sample_data(META)

sample_data(Endo_Control_PHY_rel)<-METADATA
sample_data(Endo_Control_PHY_rel)


# dbRDA
ORDCCA<-ordinate(Endo_Control_PHY_rel,"CAP",formula = ~ Genotype+Treatment+PSYield.y+Chla.y+Zoox.y+Protein.y)
ORDCCA # explains  32.346% of observed variation

library(vegan)
ANOVA<-anova.cca(ORDCCA, by="term", permutations = 10000) # 
ANOVA #sig
capture.output(ANOVA,file="Endo_Genotype_RDA.doc")

ORDCCA<-ordinate(Endo_Control_PHY_rel,"CAP",formula = ~ Genotype)
ORDCCA # explains 26.368% of observed variation

p0 = plot_ordination(Endo_Control_PHY_rel, ORDCCA, color="Genotype")+
  scale_color_brewer(palette = "Paired")+
  theme_classic(base_size = 16, base_family = "Helvetica")+
  geom_hline(yintercept=c(0,0), linetype="dotted")+
  geom_vline(xintercept=c(0,0), linetype="dotted")+
  geom_point(size=3)+
  xlab("db RDA1 [23.5%]")+
  ylab("db RDA2 [4.3%]")
p0

postscript(file = "RDAplot_EndoGenotype.eps", width =6, height = 5)
print(p0)
dev.off()

pdf('RDAplot_EndoGenotype.pdf', width=6, height=5)
print(p0)
graphics.off()

### VARPAT VennDiagramm####
library(vegan)
varspec<-(as(otu_table(Endo_Control_PHY_rel), "matrix"))
if(taxa_are_rows(Endo_Control_PHY_rel)){varspec <- t(varspec)}
varspecBray<-vegdist(varspec, method = "bray")
varchem<-data.frame(sample_data(Endo_Control_PHY_rel))
mod <- varpart(varspecBray, ~ Genotype, ~  Treatment+SamplingTimepoint, ~PSYield.y+Chla.y+Zoox.y+Protein.y,
               data=varchem)
mod
showvarparts(3)
plot(mod)


#### Figure 3c - Endozoicomonas stability/Genotype ####
Endo<-subset_taxa(PHYLOSEQ_TABLE_control, Genus=="D_5__Endozoicomonas")

TEMP<-read.csv("sample_metadata_new2.csv", header = TRUE, sep=",", dec=".", strip.white = TRUE)
str(TEMP)
TEMP<-TEMP %>% filter(Treatment!="field_control")
rownames(TEMP)<-TEMP$Sample_ID
METADATA<-sample_data(TEMP)

sample_data(Endo)<-METADATA
sample_data(Endo)

Endomelt<-psmelt(Endo)
head(Endomelt)
Endomelt$Treatment_Timepoint<-factor(Endomelt$Treatment_Timepoint, levels=c("control t1", "control t3","control t4","control t5",
                                                                            "ambient stress t1", "ambient stress t3", "ambient stress t4","ambient stress t5",
                                                                            "2100 stress t1","2100 stress t3","2100 stress t4","2100 stress t5"))

Endomelt$Treatment<-factor(Endomelt$Treatment, levels=c("control","stress_ambient","stress_2100"))

relabundplot<-ggplot(Endomelt, aes(Treatment_Timepoint, Abundance))+
  geom_col(aes(color=Treatment, fill=Treatment))+
  scale_color_manual(values=c("lightgrey","grey","darkgrey"))+
  scale_fill_manual(values=c("lightgrey","grey","darkgrey"))+
  theme_classic(base_family = "Helvetica", base_size = 16)+
  facet_wrap(~Genotype, ncol=3)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),text = element_text(size=14), legend.position = "none")+
  scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
  coord_cartesian(ylim=c(0,100))+
  labs(y = "rel.abundance [%]")

relabundplot


### Oligo 100 percent ####
sample_data(Endo)
variable1 = as.character(get_variable(Endo, "Treatment"))
variable2 = as.character(get_variable(Endo, "SamplingTimepoint"))
sample_data(Endo)$Treat_Time <- mapply(paste0, variable1, variable2, 
                                       collapse = "_")
mergedEndo<-merge_samples(Endo, "Treat_Time")

# others
othersotus=names(sort(taxa_sums(Endo), TRUE)[12:133])
othersotus

## merge others
mergedEndonew<-merge_taxa(Endo, othersotus)
tax_table(mergedEndonew)

phy <- transform_sample_counts(mergedEndonew, function(x)100*x/sum(x))

sample_data(phy)$Treatment_Timepoint<-factor(sample_data(phy)$Treatment_Timepoint, levels=c("control t1", "control t3","control t4","control t5",
                                                                                            "ambient stress t1", "ambient stress t3", "ambient stress t4","ambient stress t5",
                                                                                            "2100 stress t1","2100 stress t3","2100 stress t4","2100 stress t5"))

sample_data(phy)$Treatment<-factor(sample_data(phy)$Treatment, levels=c("control","stress_ambient","stress_2100"))

oligoplot<-plot_bar(phy, x="Treatment_Timepoint", fill="Feature_ID")+
  geom_bar(aes(fill=Feature_ID, color=Feature_ID), stat="identity", position="stack")+
  theme_classic(base_family = "Helvetica", base_size = 16)+
  scale_fill_brewer(palette="Spectral", na.value="grey")+
  scale_colour_brewer(palette = "Spectral", na.value="grey")+
  facet_wrap(~Genotype)+
  theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank(), axis.text.x = element_text(angle=90,vjust = 0.5, hjust=1))+
  theme( text = element_text(size=14))+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_continuous(labels=scaleFUN, expand = c(0, 0))+
  labs(y = "rel.abundance [%]", x="")+
  theme(legend.position="bottom")+
  guides(fill=guide_legend(nrow=4,byrow=TRUE))

oligoplot
relabundplot

library("cowplot")
combplot<-ggdraw() +
  draw_plot(relabundplot, x = 0, y = 0.60, width = 1, height = .40)+ 
  draw_plot(oligoplot, x = 0, y = 0, width = 1, height = .60) 
draw_plot_label(label = c("B", "C"), size = 15,
                x = c(0, 0), y = c(1, 0.7))

combplot


pdf('Endototal.pdf', width=8, height=10)
print(combplot)
graphics.off()  
postscript(file = "Endototal.eps", width =8, height = 20)
print(combplot)
dev.off()

# genotype A ####
EndoA<-subset_samples(Endo, Genotype=="A")

Endomelt<-psmelt(EndoA)
head(Endomelt)
Endomelt$Treatment_Timepoint<-factor(Endomelt$Treatment_Timepoint, levels=c("control t1", "control t3","control t4","control t5",
                                                                            "ambient stress t1", "ambient stress t3", "ambient stress t4","ambient stress t5",
                                                                            "2100 stress t1","2100 stress t3","2100 stress t4","2100 stress t5"))

Endomelt$Treatment<-factor(Endomelt$Treatment, levels=c("control","stress_ambient","stress_2100"))

relabundplot<-ggplot(Endomelt, aes(Treatment_Timepoint, Abundance))+
  geom_col(aes(color=Treatment, fill=Treatment))+
  scale_color_manual(values=c("lightgrey","grey","darkgrey"))+
  scale_fill_manual(values=c("lightgrey","grey","darkgrey"))+
  theme_classic(base_family = "Helvetica", base_size = 5)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),text = element_text(size=14), legend.position = "none")+
  scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
  coord_cartesian(ylim=c(0,100))+
  labs(y = "rel.abundance [%]")

relabundplot


### Oligo 100 percent ####
sample_data(Endo)
variable1 = as.character(get_variable(Endo, "Treatment"))
variable2 = as.character(get_variable(Endo, "SamplingTimepoint"))
sample_data(Endo)$Treat_Time <- mapply(paste0, variable1, variable2, 
                                       collapse = "_")
mergedEndo<-merge_samples(Endo, "Treat_Time")

# others
othersotus=names(sort(taxa_sums(Endo), TRUE)[12:133])
othersotus

## merge others
mergedEndonew<-merge_taxa(EndoA, othersotus)
tax_table(mergedEndonew)

phy <- transform_sample_counts(mergedEndonew, function(x)100*x/sum(x))

sample_data(phy)$Treatment_Timepoint<-factor(sample_data(phy)$Treatment_Timepoint, levels=c("control t1", "control t3","control t4","control t5",
                                                                                            "ambient stress t1", "ambient stress t3", "ambient stress t4","ambient stress t5",
                                                                                            "2100 stress t1","2100 stress t3","2100 stress t4","2100 stress t5"))

sample_data(phy)$Treatment<-factor(sample_data(phy)$Treatment, levels=c("control","stress_ambient","stress_2100"))

oligoplot<-plot_bar(phy, x="Treatment_Timepoint", fill="Feature_ID")+
  geom_bar(aes(fill=Feature_ID, color=Feature_ID), stat="identity", position="stack")+
  theme_classic(base_family = "Helvetica", base_size = 5)+
  scale_fill_brewer(palette="Spectral", na.value="grey")+
  scale_colour_brewer(palette = "Spectral", na.value="grey")+
  theme(legend.position="none", legend.direction="horizontal", legend.title = element_blank(), axis.text.x = element_text(angle=90,vjust = 0.5, hjust=1))+
  theme(strip.text.x = element_blank(), text = element_text(size=14))+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_continuous(labels=scaleFUN, expand = c(0, 0))+
  labs(y = "rel.abundance [%]", x="")+
  guides(fill=guide_legend(nrow=4,byrow=TRUE))

oligoplot
relabundplot

library("cowplot")
combplot<-ggdraw() +
  draw_plot(relabundplot, x = 0, y = 0.55, width = 1, height = .40)+ 
  draw_plot(oligoplot, x = 0, y = 0, width = 1, height = .60) 
draw_plot_label(label = c("B", "C"), size = 15,
                x = c(0, 0), y = c(1, 0.7))

combplot


pdf('EndoA.pdf', width=5, height=5)
print(combplot)
graphics.off()  

# gentoype B ####
EndoA<-subset_samples(Endo, Genotype=="B")

Endomelt<-psmelt(EndoA)
head(Endomelt)
Endomelt$Treatment_Timepoint<-factor(Endomelt$Treatment_Timepoint, levels=c("control t1", "control t3","control t4","control t5",
                                                                            "ambient stress t1", "ambient stress t3", "ambient stress t4","ambient stress t5",
                                                                            "2100 stress t1","2100 stress t3","2100 stress t4","2100 stress t5"))

Endomelt$Treatment<-factor(Endomelt$Treatment, levels=c("control","stress_ambient","stress_2100"))

relabundplot<-ggplot(Endomelt, aes(Treatment_Timepoint, Abundance))+
  geom_col(aes(color=Treatment, fill=Treatment))+
  scale_color_manual(values=c("lightgrey","grey","darkgrey"))+
  scale_fill_manual(values=c("lightgrey","grey","darkgrey"))+
  theme_classic(base_family = "Helvetica", base_size = 5)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),text = element_text(size=14), legend.position = "none")+
  scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
  coord_cartesian(ylim=c(0,100))+
  labs(y = "rel.abundance [%]")

relabundplot


### Oligo 100 percent ####
sample_data(Endo)
variable1 = as.character(get_variable(Endo, "Treatment"))
variable2 = as.character(get_variable(Endo, "SamplingTimepoint"))
sample_data(Endo)$Treat_Time <- mapply(paste0, variable1, variable2, 
                                       collapse = "_")
mergedEndo<-merge_samples(Endo, "Treat_Time")

# others
othersotus=names(sort(taxa_sums(Endo), TRUE)[12:133])
othersotus

## merge others
mergedEndonew<-merge_taxa(EndoA, othersotus)
tax_table(mergedEndonew)

phy <- transform_sample_counts(mergedEndonew, function(x)100*x/sum(x))

sample_data(phy)$Treatment_Timepoint<-factor(sample_data(phy)$Treatment_Timepoint, levels=c("control t1", "control t3","control t4","control t5",
                                                                                            "ambient stress t1", "ambient stress t3", "ambient stress t4","ambient stress t5",
                                                                                            "2100 stress t1","2100 stress t3","2100 stress t4","2100 stress t5"))

sample_data(phy)$Treatment<-factor(sample_data(phy)$Treatment, levels=c("control","stress_ambient","stress_2100"))

oligoplot<-plot_bar(phy, x="Treatment_Timepoint", fill="Feature_ID")+
  geom_bar(aes(fill=Feature_ID, color=Feature_ID), stat="identity", position="stack")+
  theme_classic(base_family = "Helvetica", base_size = 5)+
  scale_fill_brewer(palette="Spectral", na.value="grey")+
  scale_colour_brewer(palette = "Spectral", na.value="grey")+
  theme(legend.position="none", legend.direction="horizontal", legend.title = element_blank(), axis.text.x = element_text(angle=90,vjust = 0.5, hjust=1))+
  theme(strip.text.x = element_blank(), text = element_text(size=14))+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_continuous(labels=scaleFUN, expand = c(0, 0))+
  labs(y = "rel.abundance [%]", x="")+
  guides(fill=guide_legend(nrow=4,byrow=TRUE))

oligoplot
relabundplot

library("cowplot")
combplot<-ggdraw() +
  draw_plot(relabundplot, x = 0, y = 0.55, width = 1, height = .40)+ 
  draw_plot(oligoplot, x = 0, y = 0, width = 1, height = .60) 
draw_plot_label(label = c("B", "C"), size = 15,
                x = c(0, 0), y = c(1, 0.7))

combplot


pdf('EndoB.pdf', width=5, height=5)
print(combplot)
graphics.off() 

# genotype C ####
EndoA<-subset_samples(Endo, Genotype=="C")

Endomelt<-psmelt(EndoA)
head(Endomelt)
Endomelt$Treatment_Timepoint<-factor(Endomelt$Treatment_Timepoint, levels=c("control t1", "control t3","control t4","control t5",
                                                                            "ambient stress t1", "ambient stress t3", "ambient stress t4","ambient stress t5",
                                                                            "2100 stress t1","2100 stress t3","2100 stress t4","2100 stress t5"))

Endomelt$Treatment<-factor(Endomelt$Treatment, levels=c("control","stress_ambient","stress_2100"))

relabundplot<-ggplot(Endomelt, aes(Treatment_Timepoint, Abundance))+
  geom_col(aes(color=Treatment, fill=Treatment))+
  scale_color_manual(values=c("lightgrey","grey","darkgrey"))+
  scale_fill_manual(values=c("lightgrey","grey","darkgrey"))+
  theme_classic(base_family = "Helvetica", base_size = 5)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),text = element_text(size=14), legend.position = "none")+
  scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
  coord_cartesian(ylim=c(0,100))+
  labs(y = "rel.abundance [%]")

relabundplot


### Oligo 100 percent ####
sample_data(Endo)
variable1 = as.character(get_variable(Endo, "Treatment"))
variable2 = as.character(get_variable(Endo, "SamplingTimepoint"))
sample_data(Endo)$Treat_Time <- mapply(paste0, variable1, variable2, 
                                       collapse = "_")
mergedEndo<-merge_samples(Endo, "Treat_Time")

# others
othersotus=names(sort(taxa_sums(Endo), TRUE)[12:133])
othersotus

## merge others
mergedEndonew<-merge_taxa(EndoA, othersotus)
tax_table(mergedEndonew)

phy <- transform_sample_counts(mergedEndonew, function(x)100*x/sum(x))

sample_data(phy)$Treatment_Timepoint<-factor(sample_data(phy)$Treatment_Timepoint, levels=c("control t1", "control t3","control t4","control t5",
                                                                                            "ambient stress t1", "ambient stress t3", "ambient stress t4","ambient stress t5",
                                                                                            "2100 stress t1","2100 stress t3","2100 stress t4","2100 stress t5"))

sample_data(phy)$Treatment<-factor(sample_data(phy)$Treatment, levels=c("control","stress_ambient","stress_2100"))

oligoplot<-plot_bar(phy, x="Treatment_Timepoint", fill="Feature_ID")+
  geom_bar(aes(fill=Feature_ID, color=Feature_ID), stat="identity", position="stack")+
  theme_classic(base_family = "Helvetica", base_size = 5)+
  scale_fill_brewer(palette="Spectral", na.value="grey")+
  scale_colour_brewer(palette = "Spectral", na.value="grey")+
  theme(legend.position="none", legend.direction="horizontal", legend.title = element_blank(), axis.text.x = element_text(angle=90,vjust = 0.5, hjust=1))+
  theme(strip.text.x = element_blank(), text = element_text(size=14))+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_continuous(labels=scaleFUN, expand = c(0, 0))+
  labs(y = "rel.abundance [%]", x="")+
  guides(fill=guide_legend(nrow=4,byrow=TRUE))

oligoplot
relabundplot

library("cowplot")
combplot<-ggdraw() +
  draw_plot(relabundplot, x = 0, y = 0.55, width = 1, height = .45)+ 
  draw_plot(oligoplot, x = 0, y = 0, width = 1, height = .45) 
draw_plot_label(label = c("B", "C"), size = 15,
                x = c(0, 0), y = c(1, 0.7))

combplot


pdf('EndoC.pdf', width=5, height=5)
print(combplot)
graphics.off() 








#### PIECHART ####
OTU<-as.data.frame(otu_table(PHYLOSEQ_TABLE_control))
library(data.table)
setDT(OTU, keep.rownames = TRUE)[]
OTU<-dplyr::rename(OTU, Feature_ID=rn)

TAXA<-as.data.frame(tax_table(PHYLOSEQ_TABLE_control))

library(dplyr)
library(tidyr)
Piecharttable<-right_join(OTU, TAXA, by="Feature_ID")

Endos<-Piecharttable %>% filter(Genus=="D_5__Endozoicomonas") %>%
  filter(Feature_ID!="c0e6c66cfd0673747e59936d6c050407") %>% 
  tidyr::gather("Sample","Abundance", 2:97) %>% 
  group_by(Sample) %>% 
  summarise(MEAN_Endo=sum(Abundance)) %>% 
  summarise(MEAN_Endo_total=mean(MEAN_Endo))

Core<-Piecharttable %>%
  filter(Feature_ID=="c0e6c66cfd0673747e59936d6c050407") %>% 
  tidyr::gather("Sample","Abundance", 2:97) %>% 
  group_by(Sample) %>% 
  summarise(MEAN_Endo=sum(Abundance)) %>% 
  summarise(MEAN_Endo_total=mean(MEAN_Endo))


#### Boxplot - rel.abund/Genotype ####
Endos<-Piecharttable %>% filter(Genus=="D_5__Endozoicomonas") %>% 
  filter(Feature_ID!="c0e6c66cfd0673747e59936d6c050407") %>% 
  tidyr::gather("SAMPLE","Abundance", 2:97) %>% 
  group_by(SAMPLE) %>%summarise(Endozoicomonas=sum(Abundance))
Core<-Piecharttable %>% filter(Feature_ID=="c0e6c66cfd0673747e59936d6c050407") %>% 
  tidyr::gather("SAMPLE","Abundance", 2:97) %>% 
  group_by(SAMPLE) %>%summarise(core=sum(Abundance))

SAMPLEDATA<-data.frame(sample_data(PHYLOSEQ_TABLE_control))
Endo_table<-inner_join(SAMPLEDATA,Endos, by="SAMPLE")
Endo_table<-inner_join(Endo_table, Core, by="SAMPLE")
Endo_table_melt<-Endo_table %>% select(-PSYield, -Chla, -Zoox, -Protein) 
Endo_table_melt<-melt(Endo_table_melt)

plot1<-ggplot(Endo_table_melt, aes(x=variable, y=value, fill=variable))+
  facet_wrap(~Genotype)+
  geom_boxplot()+
  theme_classic(base_size = 16, base_family = "Helvetica")+
  scale_fill_brewer(palette = "Paired")+
  ylab("rel.abundance [%]")+
  xlab("")+
  theme(legend.position = "none")+
  scale_fill_manual(values=c("#1f78b4","#a6cee3"))
plot1

pdf('EndoCore_Abundance_Genotype.pdf', width=6, height=6)
print(plot1)
graphics.off()

ENDOANOVA<-aov(value~Genotype+variable, data=Endo_table_melt)
summary(ENDOANOVA) #sig 1.03e-05 ***

ENDOposthoc<-TukeyHSD(ENDOANOVA)
ENDOposthoc


#### indval Genotype ####
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(ggpubr)
options(expressions = 500000)

#prep data for indval
METADATA<-as(sample_data(Endo_Control_PHY_rel), "data.frame")
METADATA<-rownames_to_column(METADATA, var="Sample_ID")

OTU_table<-as(otu_table(Endo_Control_PHY_rel), "matrix")
if(taxa_are_rows(Endo_Control_PHY_rel)){OTU_table <- t(OTU_table)}
OTU_table_df = as.data.frame(OTU_table)
OTU_table_df<-rownames_to_column(OTU_table_df, var="Sample_ID")

# join datasets
library(dplyr)
library(tidyr)
OTU_table_IndVal<-right_join(META, OTU_table_df, by="Sample_ID")

#indval
library(labdsv)
library(MASS)
library(vegan)
library(cluster)
library(indicspecies)
library(permute)
(INDVAL_OTUs_species_GC=(as.data.frame(OTU_table_IndVal[,14:146])))

#Cluster 
(INDVAL_Groups_Origin_GC=(as.character(OTU_table_IndVal$Genotype)))
INDVAL_Origin_GC=multipatt(INDVAL_OTUs_species_GC, INDVAL_Groups_Origin_GC, func="IndVal.g", duleg=FALSE, control=how(nperm=10000))
summary(INDVAL_Origin_GC, indvalcomp=FALSE)
head(summary.multipatt(INDVAL_Origin_GC, alpha = 0.05, indvalcomp=FALSE)) #At=0.80, Bt=0.80 
IndValdf<-data.frame(INDVAL_Origin_GC$sign)
IndValdf<-rownames_to_column(IndValdf, var="Feature_ID")
IndValdf<-IndValdf %>% filter(p.value<=0.05) %>% dplyr:: rename(A=s.A) %>% 
  dplyr:: rename(B=s.B)%>%
  dplyr:: rename(D=s.D)%>%
  dplyr:: rename(E=s.E)%>%
  dplyr:: rename(`F`=s.F)%>%
  dplyr:: rename(H=s.H)%>%
  dplyr:: rename(I=s.I)%>%
  dplyr:: rename(J=s.J)
capture.output(IndValdf,file="IndVal_GenotypeControlEndo.doc")


## bubble plot ##
Endodf<-psmelt(Endo_Control_PHY_rel)
#### SE function #### 
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

### calculate Mean, SE, SD, CI ####
SummarySE<-summarySE(Endodf, measurevar="Abundance", groupvars=c("Genotype","OTU"))

# Bubble Plot ####
library(ggplot2)
plot_Indicator<-ggplot(SummarySE[which(SummarySE$Abundance>0),],
                       aes(x=Genotype,y=OTU, fill="Genotype"))+
  geom_point(aes(size=Abundance, fill=Genotype), shape=21)+
  scale_fill_brewer(palette = "Paired")+
  scale_size_area(max_size = 6)+
  theme(panel.grid=element_blank())+
  theme_classic(base_family = "Helvetica", base_size = 10)+
  guides(fill=guide_legend(title="host genotype"), size=guide_legend(title = "mean rel.abundance [%]"))+
  xlab("host genotype")+
  ylab("Endozoicmonas ASV_ID")

plot_Indicator

postscript(file = "Bubble_EndoGenotype.eps", width =10, height = 10)
print(plot_Indicator)
dev.off()

pdf('Bubble_EndoGenotype.pdf', width=10, height=10)
print(plot_Indicator)
graphics.off()


# ANOSIM####
##Genotype ##
library(vegan)
Genotype_group<-get_variable(??, "Genotype")
HOSTANOSIM<-anosim(distance(??,"bray"), Genotype_group)
HOSTANOSIM$signif #p=0.001
HOSTANOSIM$statistic #R=0.3695626

Treatment_group<-get_variable(test, "Treatment")
HOSTANOSIM<-anosim(distance(test,"bray"), Treatment_group)
HOSTANOSIM$signif #p=0.001
HOSTANOSIM$statistic #R=0.3695626



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

sample_data(test)<-METADATA
sample_data(test)

library(dplyr)
conMETA<-sample_data(test) %>% dplyr::filter(Treatment=="control")%>%select(PSYield.y, Chla.y, Zoox.y, Protein.y)

bray.dist.otu<-phyloseq::distance(test,'bray') #using rarefied species data

bray.dist.host<-vegdist(conMETA, method = "bray") #using scaled host health proxy data

bray.mantel<-mantel(bray.dist.otu, bray.dist.host, permutations = 10000, method="pearson")

bray.mantel