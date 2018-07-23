##### Figure 2a-b ####
# Host health proxy response

load(file="PHYLOSEQ_TABLE.RData")
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(vegan)


PHYLOSEQ_TABLE_control=subset_samples(PHYLOSEQ_TABLE, Treatment!="field_control")
PHYLOSEQ_TABLE_control

#### z score ####
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

sample_data(PHYLOSEQ_TABLE_control)<-METADATA
sample_data(PHYLOSEQ_TABLE_control)

#### re-arrange data ####
longTemp<-META%>% select(SamplingTimepoint, Protein.y, Zoox.y, Chla.y, PSYield.y, Treatment) %>% tidyr::gather(key=environmental_data, value=zscore, -SamplingTimepoint, -Treatment)
longTempMEAN<-longTemp %>% group_by(SamplingTimepoint, Treatment,environmental_data) %>% summarise(MEAN=mean(zscore), SD=sd(zscore))
longTempMEAN$Treatment<-factor(longTempMEAN$Treatment, levels=c("field_control", "control", "stress_ambient", "stress_2100"))

### plot graph ####
library(ggplot2)
plot<-ggplot(longTempMEAN, aes(x=SamplingTimepoint, y=`MEAN`, colour=environmental_data))+
  theme_classic(base_size = 16, base_family = "Helvetica")+
  geom_hline(yintercept=c(0,0), linetype="dotted")+
  geom_errorbar(aes(ymax=(MEAN+SD), ymin=(MEAN-SD)), colour="black", width=0.2)+
  geom_point(size=3)+
  geom_line(aes(group=Treatment))+
  scale_color_brewer(palette = "Set1")+
  facet_wrap(~Treatment+environmental_data)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_wrap(~Treatment+environmental_data)+
  ylab("z-scores")+
  xlab("days")

plot  
pdf('Fig_1b.pdf', width=10, height=8)
print(plot)
graphics.off()


postscript(file = "Fig_1b.eps", width =10, height = 8)
print(plot)
dev.off()

### ANOVA ####
ANOVA<-aov(zscore~Treatment*SamplingTimepoint+Error(environmental_data), data=longTemp)
summary(ANOVA)

CHLA<-longTemp %>% filter(environmental_data=="Chla.y")
ANOVACHla<-aov(zscore~Treatment, data=CHLA)
summary(ANOVACHla) # not sig
capture.output(summary(ANOVACHla),file="Chla_Treatment_ANOVA.doc")

PROTEIN<-longTemp %>% filter(environmental_data=="Protein.y")
ANOVAPROTEIN<-aov(zscore~Treatment, data=PROTEIN)
summary(ANOVAPROTEIN) # Treatment significant #0.00468 **
capture.output(summary(ANOVAPROTEIN),file="Protein_Treatment_ANOVA.doc")

TukeyHSD(ANOVAPROTEIN) # sig. lower prot conc in stress 2100 treatment
capture.output(TukeyHSD(ANOVAPROTEIN),file="Protein_Treatment_Tukey.doc")

ZOOX<-longTemp %>% filter(environmental_data=="Zoox.y")
ANOVAZOOX<-aov(zscore~Treatment, data=ZOOX)
summary(ANOVAZOOX) # Treat 0.000479 ***

TukeyHSD(ANOVAZOOX) # sig. lower Zoox conc in stress 2100 treatment
capture.output(TukeyHSD(ANOVAZOOX),file="Zoox_Treatment_Tukey.doc")

PSYield<-longTemp %>% filter(environmental_data=="PSYield.y")
ANOVAPSYield<-aov(zscore~Treatment, data=PSYield)
summary(ANOVAPSYield) # Treatment significant 2.71e-05 ***

TukeyHSD(ANOVAPSYield) # sig. higher PS in stress 2100 treatment
capture.output(TukeyHSD(ANOVAPSYield),file="PSyield_Treatment_Tukey.doc")


#### Genotype ####
longTemp<-META%>% select(Genotype, Protein.y, Zoox.y, Chla.y, PSYield.y, Treatment) %>% tidyr::gather(key=environmental_data, value=zscore, -Genotype, -Treatment)
longTempMEAN<-longTemp %>% group_by(Genotype, Treatment,environmental_data) %>% summarise(MEAN=mean(zscore), SD=sd(zscore))
longTempMEAN$Treatment<-factor(longTempMEAN$Treatment, levels=c("field_control", "control", "stress_ambient", "stress_2100"))

library(ggplot2)
plot<-ggplot(longTempMEAN, aes(x=Genotype, y=`MEAN`, colour=environmental_data))+
  theme_classic(base_size = 16, base_family = "Helvetica")+
  geom_hline(yintercept=c(0,0), linetype="dotted")+
  geom_errorbar(aes(ymax=(MEAN+SD), ymin=(MEAN-SD)), colour="black", width=0.2)+
  geom_point(size=3)+
  scale_color_brewer(palette = "Set1")+
  facet_wrap(~Treatment+environmental_data)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_wrap(~Treatment+environmental_data)+
  ylab("z-scores")+
  xlab("genotype")
plot  
pdf('Fig_1c.pdf', width=10, height=8)
print(plot)
graphics.off()


postscript(file = "Fig_1c.eps", width =10, height = 8)
print(plot)
dev.off()

### ANOVA ####
ANOVA<-aov(zscore~Genotype+Error(Treatment), data=longTemp)
summary(ANOVA) # 2.06e-06 ***

CHLA<-longTemp %>% filter(environmental_data=="Chla.y")
ANOVACHla<-aov(zscore~Genotype, data=CHLA)
summary(ANOVACHla) # 0.00864 **

TukeyHSD(ANOVACHla) # 

PROTEIN<-longTemp %>% filter(environmental_data=="Protein.y")
ANOVAPROTEIN<-aov(zscore~Genotype, data=PROTEIN)
summary(ANOVAPROTEIN) # Genotype not sign # 0.0906 .


ZOOX<-longTemp %>% filter(environmental_data=="Zoox.y")
ANOVAZOOX<-aov(zscore~Genotype, data=ZOOX)
summary(ANOVAZOOX) # Genotype 0.000186 ***

TukeyHSD(ANOVAZOOX) # 

PSYield<-longTemp %>% filter(environmental_data=="PSYield.y")
ANOVAPSYield<-aov(zscore~Genotype, data=PSYield)
summary(ANOVAPSYield) # Treatment significant 0.0106 *