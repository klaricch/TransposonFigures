#!/usr/bin/R
# this script
# 1) plots total transposons vs strain per insertions, references, and absences per class
# 2) plot histogram of transposons per strain
# USE: TE_vs_DR.R

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(grid)
setwd("/Users/kristen/Documents/transposon_figure_data/data")
summarydata <- read.table("T_kin_C_matrix_full.txt",header=TRUE)
#remove ZERO_new traits
summarydata<-subset(summarydata,!grepl('^ZERO_new', summarydata$trait))
summarydata<-subset(summarydata,!grepl('^coverage', summarydata$trait))
#clean trait names
summarydata$trait <- gsub("_C$" ,"",summarydata$trait)
summarydata$trait <- gsub("^ONE_new" ,"new",summarydata$trait)

classdata<- read.table("CtCp_all_nonredundant.txt",header=TRUE)
names(classdata)<-c("chr","start","end","TE","orientation","method","strain","class")

# add te class info to summarydata(new_TRANS_end_tes will be removed)
classdata$family<- stringr::str_split_fixed(classdata$TE, regex("_(non-)?reference"),2)[,1]
classdata$family<- paste(stringr::str_split_fixed(classdata$family, "_",4)[,3],stringr::str_split_fixed(classdata$family, "_",4)[,4],sep="_")
classdata$family <- gsub("_$" ,"",classdata$family)
classdata$family <- gsub("_non-reference(.*)$" ,"",classdata$family)
classdata<-mutate(classdata, trait=paste(method,"TRANS",family,sep="_"))
class_subset <- classdata %>% distinct(trait) %>% select(trait,class)
summarydata <-merge(summarydata, class_subset, by="trait")

#names(summarydata)
summarydata<-gather(summarydata, "sample","value",2:(ncol(summarydata)-1))
#summarydata<-rename(summarydata,total_tes=value)
tail(summarydata)

#new column that specifies what caller was used
summarydata$method<- stringr::str_split_fixed(summarydata$trait, "_TRANS_",2)[,1]
#new column that specifies TE family
summarydata$transposon<- stringr::str_split_fixed(summarydata$trait, "_TRANS_",2)[,2]
#summarydata<-filter(summarydata,transposon=="total")



families<-summarydata


#reformat the data
summarydata <- summarydata %>% group_by(sample,class, method) %>% summarise(XX=sum(value,na.rm=TRUE))
hist_data<-summarydata
total_absence<-filter(summarydata,method=="absent")
total_reference<-filter(summarydata,method=="reference")
total_insertion<-filter(summarydata,method=="new")





#SCATTER
initial_merge <- merge(total_absence,total_reference, by=c("sample","class"))

final_merge <- merge(initial_merge,total_insertion, by=c("sample","class"))

names(final_merge)<-c("sample", "class", "method.x","total_absences",	"method.y", "total_references",	"method.z","total_insertions")


#remove points in "unknown"classification
final_merge<-filter(final_merge, class !="unknown")


#revalue classes
final_merge$class <- factor(final_merge$class,
                            levels = c("dnatransposon", "retrotransposon","unknown"),
                            labels = c("DNA Transposon", "Retrotransposon", "Unknown"))
####
retro<-filter(final_merge, class=="Retrotransposon")
dna<-filter(final_merge, class=="DNA Transposon")
R_max_insertions<-max(retro$total_insertions)
R_max_absences<-max(retro$total_absences)
R_max_references<-max(retro$total_references)
D_max_insertions<-max(dna$total_insertions)
D_max_absences<-max(dna$total_absences)
D_max_references<-max(dna$total_references)


#R_max_references R_max_insertions R_max_absences
#D_max_insertions D_max_references D_max_absences


#1 ABSENCE vs INSERTION
# ADD IN CLASS
#spearman correlation

R_correlation<-cor.test(retro$total_absences, retro$total_insertions,method="spearman",exact=FALSE)
R_rho<-round(R_correlation$estimate,3)
D_correlation<-cor.test(dna$total_absences, dna$total_insertions,method="spearman",exact=FALSE)
D_rho<-round(D_correlation$estimate,3)

Rla <- paste("italic(rho) == ", R_rho)
Dla <- paste("italic(rho) == ", D_rho)
m <- ggplot(final_merge, aes(x=total_insertions, y=total_absences))+
  facet_wrap(~class,ncol=3,scale="free")
m <- m + geom_point(size=1.25) +
  geom_smooth(method="lm",se=FALSE,col="red")+
  geom_abline(slope=1,linetype="dashed",colour="gray52")+
  geom_text(data = subset(final_merge, class=="DNA Transposon"), x=.2*D_max_insertions, y=.9*D_max_insertions,label=Dla,parse=TRUE, colour="red",size=2.5)+
  geom_text(data = subset(final_merge, class=="Retrotransposon"), x=.2*R_max_insertions, y=.9*R_max_insertions,label=Rla,parse=TRUE, colour="red",size=2.5)+
  geom_point(data = subset(final_merge, class=="Retrotransposon"), aes(x=R_max_insertions,y=R_max_insertions),alpha=0)+
  geom_point(data = subset(final_merge, class=="DNA Transposon"), aes(x=D_max_insertions,y=D_max_insertions),alpha=0)+
  geom_point(aes(x=0,y=0),alpha=0)+
  
  theme(strip.text.x = element_text(size = 9, colour = "black"),
        strip.background = element_blank(),
        legend.position=c(.90,0.75),
        legend.background = element_rect(fill=FALSE),
        legend.text=element_text(size=9),
        panel.background = element_rect(fill = "white"),
        panel.margin = unit(1, "lines"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=9),
        axis.text.x = element_text(colour = "black",size=9),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9))+
  guides(fill=FALSE) +
  labs(x = "Insertion Events", y = "Absence Events")
m
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="Absence_vs_Insertion_DR.tiff",
       dpi=300,
       width=7.5,
       height=3.5,
       units="in")

#2 ABSENCE vs REFERENCE
#spearman correlation
R_correlation<-cor.test(retro$total_absences, retro$total_references,method="spearman",exact=FALSE)
R_rho<-round(R_correlation$estimate,3)
D_correlation<-cor.test(dna$total_absences, dna$total_references,method="spearman",exact=FALSE)
D_rho<-round(D_correlation$estimate,3)

Rla <- paste("italic(rho) == ", R_rho)
Dla <- paste("italic(rho) == ", D_rho)




m <- ggplot(final_merge, aes(x=total_references, y=total_absences))+
    facet_wrap(~class,ncol=3,scale="free")
m <- m + geom_point(size=1.25) +
  geom_smooth(method="lm",se=FALSE,col="red")+
  geom_abline(slope=1,linetype="dashed",colour="gray52")+
  geom_text(data = subset(final_merge, class=="DNA Transposon"), x=.2*D_max_references, y=.9*D_max_references,label=Dla,parse=TRUE, colour="red",size=2.5)+
  geom_text(data = subset(final_merge, class=="Retrotransposon"), x=.2*R_max_references, y=.9*R_max_references,label=Rla,parse=TRUE, colour="red",size=2.5)+
  geom_point(data = subset(final_merge, class=="Retrotransposon"), aes(x=R_max_references,y=R_max_references),alpha=0)+
  geom_point(data = subset(final_merge, class=="DNA Transposon"), aes(x=D_max_references,y=D_max_references),alpha=0)+
  geom_point(aes(x=0,y=0),alpha=0)+
  
  theme(strip.text.x = element_text(size = 9, colour = "black"),
        strip.background = element_blank(),
        legend.position=c(.90,0.75),
        legend.background = element_rect(fill=FALSE),
        legend.text=element_text(size=9),
        panel.background = element_rect(fill = "white"),
        panel.margin = unit(1, "lines"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=9),
        axis.text.x = element_text(colour = "black",size=9),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9))+
  guides(fill=FALSE) +
  labs(x = "Reference Events", y = "Absence Events")
m
ggsave(filename="Absence_vs_Reference_DR.tiff",
       dpi=300,
       width=7.5,
       height=3.5,
       units="in")

#3 INSERTION vs REFERENCE
#spearman correlation
R_correlation<-cor.test(retro$total_reference, retro$total_insertions,method="spearman",exact=FALSE)
R_rho<-round(R_correlation$estimate,3)
D_correlation<-cor.test(dna$total_reference, dna$total_insertions,method="spearman",exact=FALSE)
D_rho<-round(D_correlation$estimate,3)

Rla <- paste("italic(rho) == ", R_rho)
Dla <- paste("italic(rho) == ", D_rho)
m <- ggplot(final_merge, aes(x=total_insertions, y=total_references))+
  facet_wrap(~class,ncol=3,scale="free")
m <- m + geom_point(size=1.25) +
  geom_smooth(method="lm",se=FALSE,col="red")+
  geom_abline(slope=1,linetype="dashed",colour="gray52")+
  geom_text(data = subset(final_merge, class=="DNA Transposon"), x=.2*D_max_insertions, y=.9*D_max_insertions,label=Dla,parse=TRUE, colour="red",size=2.5)+
  geom_text(data = subset(final_merge, class=="Retrotransposon"), x=.2*R_max_insertions, y=.9*R_max_references,label=Rla,parse=TRUE, colour="red",size=2.5)+
  geom_point(data = subset(final_merge, class=="Retrotransposon"), aes(x=R_max_references,y=R_max_references),alpha=0)+
  geom_point(data = subset(final_merge, class=="DNA Transposon"), aes(x=D_max_insertions,y=D_max_insertions),alpha=0)+
  geom_point(aes(x=0,y=0),alpha=0)+
  
  theme(strip.text.x = element_text(size = 9, colour = "black"),
        strip.background = element_blank(),
        legend.position=c(.90,0.75),
        legend.background = element_rect(fill=FALSE),
        legend.text=element_text(size=9),
        panel.background = element_rect(fill = "white"),
        panel.margin = unit(1, "lines"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=9),
        axis.text.x = element_text(colour = "black",size=9),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9))+
  guides(fill=FALSE) +
  labs(x = "Insertion Events", y = "Reference Events")
m
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="Insertion_vs_Reference_DR.tiff",
       dpi=300,
       width=7.5,
       height=3.5,
       units="in")


#CHEXK T HIS  OVER

#FAMILY
thing<-select(families, -trait)
thing<-spread(thing,method,value)

AR<-filter(thing, absent!="NA" & reference !="NA")
AI<-filter(thing, absent!="NA" & new !="NA")
RI<-filter(thing, reference!="NA" & new !="NA")

#calculate rho values
AR_corr<-AR %>% group_by(class,transposon) %>% summarise(cor=round((cor.test(absent,reference,method="spearman",exact=FALSE))$estimate,3))
AI_corr<-AI %>% group_by(class,transposon) %>% summarise(cor=round((cor.test(absent,new,method="spearman",exact=FALSE))$estimate,3))
RI_corr<-RI %>% group_by(class,transposon) %>% summarise(cor=round((cor.test(reference,new,method="spearman",exact=FALSE))$estimate,3))

#calculate mean and sd of the correlation rho values
AR_mean<-AR_corr %>% group_by(class) %>% summarise(mean=mean(cor,na.rm=TRUE),SD=sd(cor,na.rm=TRUE))
AI_mean<-AI_corr %>% group_by(class) %>% summarise(mean=mean(cor,na.rm=TRUE),SD=sd(cor,na.rm=TRUE))
RI_mean<-RI_corr %>% group_by(class) %>% summarise(mean=mean(cor,na.rm=TRUE),SD=sd(cor,na.rm=TRUE))

#
merged_AR<-merge(AR_corr, AR_mean, by="class")
merged_AI<-merge(AI_corr, AI_mean, by="class")
merged_RI<-merge(RI_corr, RI_mean, by="class")

merged_AR$comparison<-"AR"
merged_AI$comparison<-"AI"
merged_RI$comparison<-"RI"

all<-rbind(merged_AR,merged_AI,merged_RI)

all<-mutate(all, outL=ifelse(abs(cor-mean)>SD, "OUTLIER", "n"))
outliers<-filter(all,outL=="OUTLIER")
setwd("/Users/kristen/Documents/transposon_figure_data/data")
save(outliers,file="outlier_slopes.Rda")






#Histogram Versions
names(summarydata)
head(summarydata)
method_names <- list(
  'new'="Insertion",
  'reference'="Reference",
  'absent'="Absence"  
)

method_labeller <- function(variable,value){
  return(method_names[value])
}

hist_data$method <- factor(hist_data$method,
                       levels = c("new", "reference", "absent"),
                       labels = c("Insertion", "Reference", "Absence"))



m <- ggplot(hist_data, aes(x=XX,fill=class))
m <-m + geom_bar(binwidth=5)+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))+
  facet_wrap(~method,scale="free_y")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 9, colour = "black",face="bold"),
        panel.margin = unit(.25, "lines"),
        panel.border = element_rect(fill=NA,colour = "black"),
        panel.background = element_blank(),
        axis.title=element_text(size=9,face="bold"),
        axis.text.y = element_text(colour = "black",size=9),
        axis.text.x = element_text(colour = "black",size=9),
        axis.ticks =element_line(colour = "black"),
        legend.title=element_blank(),
        legend.position="none")+
  scale_fill_manual(values = c('dnatransposon' = "navy", "retrotransposon"="brown3","unknown"="goldenrod"))+
  labs(x="Transposons Per Strain", y="Count")
m

setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="Histogram_TE_per_Strain.tiff",
       dpi=300,
       width=7.5,
       height=3.5,
       units="in")





