#!/usr/bin/R
# this script
# 1) plots no. transposon events vs each other for each possible pairing
# 2) plots total transposons vs strain per insertions, references, and absences
# USE: Te_totals_distribution.R

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(cowplot)
library(grid)
setwd("/Users/kristen/Documents/transposon_figure_data/data")
summarydata <- read.table("T_kin_C_matrix_full.txt",header=TRUE)
#remove ZERO_new traits
summarydata<-subset(summarydata,!grepl('^ZERO_new', summarydata$trait))
summarydata<-subset(summarydata,!grepl('^coverage', summarydata$trait))
#clean trait names
summarydata$trait <- gsub("_C$" ,"",summarydata$trait)
summarydata$trait <- gsub("^ONE_new" ,"new",summarydata$trait)

#new column that specifies what caller was used
summarydata$method<- stringr::str_split_fixed(summarydata$trait, "_TRANS_",2)[,1]
#new column that specifies TE family
summarydata$transposon<- stringr::str_split_fixed(summarydata$trait, "_TRANS_",2)[,2]
summarydata<-filter(summarydata,transposon=="total") # this will get total ins,ref,abs calls, NOT total DNA, Retro, Unknonwn
unique(summarydata$transposon)

#names(summarydata)
summarydata<-gather(summarydata, "sample","value",2:(ncol(summarydata)-2))
summarydata<-rename(summarydata,total_tes=value)

#reformat the data
total_absence<-filter(summarydata,method=="absent")
total_reference<-filter(summarydata,method=="reference")
total_insertion<-filter(summarydata,method=="new")

#SCATTER
final_merge<- Reduce(function(x, y) merge(x, y, all=TRUE,by="sample"), list(total_absence, total_reference, total_insertion))
names(final_merge)<-c("sample", "trait.x",	"method.x",	"transposon.x",	"total_absences",	"trait.y",	"method.y",	"transposon.y",	"total_references",	"trait",	"method",	"transposon",	"total_insertions")


#1 ABSENCE vs INSERTION
#spearman correlation
correlation<-cor.test(final_merge$total_absences, final_merge$total_insertions,method="spearman",exact=FALSE)
rho<-round(correlation$estimate,3)

max_insertions<-max(final_merge$total_insertions)
max_absences<-max(final_merge$total_absences)

la <- paste("italic(rho) == ", rho)
m1 <- ggplot(final_merge, aes(x=total_insertions, y=total_absences))
m1 <- m1 + geom_point(size=1.25) + xlim(0,max_insertions)+ ylim(0,max_insertions)+
  geom_smooth(method="lm",se=FALSE,col="red")+
  geom_abline(slope=1,linetype="dashed",colour="gray52")+
  annotate("text", x=.2*max_insertions, y=.9*max_insertions,label=la,parse=TRUE, colour="red",size=2.5)+
  theme(strip.text.x = element_text(size = 9, colour = "black"),
        strip.background = element_blank(),
        legend.position=c(.90,0.75),
        legend.background = element_rect(fill=FALSE),
        legend.text=element_text(size=9),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=9),
        axis.text.x = element_text(colour = "black",size=9),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9))+
  guides(fill=FALSE) +
  labs(x = "Insertion Events", y = "Absence Events")
m1
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="Absence_vs_Insertion.tiff",
       dpi=300,
       width=4,
       height=4,
       units="in")

#2 ABSENCE vs REFERENCE
#spearman correlation
correlation<-cor.test(final_merge$total_absences, final_merge$total_references,method="spearman",exact=FALSE)
rho<-round(correlation$estimate,3)

max_references<-max(final_merge$total_references)
max_absences<-max(final_merge$total_absences)

la <- paste("italic(rho) == ", rho)
m2 <- ggplot(final_merge, aes(x=total_references, y=total_absences))
m2 <- m2 + geom_point(size=1.25) + xlim(0,max_references)+ ylim(0,max_references)+
  geom_smooth(method="lm",se=FALSE,col="red")+
  geom_abline(slope=1,linetype="dashed",colour="gray52")+
  annotate("text", x=.2*max_references, y=.9*max_references,label=la,parse=TRUE, colour="red",size=2.5)+
  theme(strip.text.x = element_text(size = 9, colour = "black"),
        strip.background = element_blank(),
        legend.position=c(.90,0.75),
        legend.background = element_rect(fill=FALSE),
        legend.text=element_text(size=9),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=9),
        axis.text.x = element_text(colour = "black",size=9),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9))+
  guides(fill=FALSE) +
  labs(x = "Reference Events", y = "Absence Events")
m2
ggsave(filename="Absence_vs_Reference.tiff",
       dpi=300,
       width=4,
       height=4,
       units="in")



#3 INSERTION vs REFERENCE
#spearman correlation
correlation<-cor.test(final_merge$total_insertions, final_merge$total_references,method="spearman",exact=FALSE)
rho<-round(correlation$estimate,3)

max_references<-max(final_merge$total_references)
max_insertions<-max(final_merge$total_insertions)
la <- paste("italic(rho) == ", rho)

max_references<-max(final_merge$total_references)
m3 <- ggplot(final_merge, aes(x=total_references, y=total_insertions))
m3 <- m3 + geom_point(size=1.25) + xlim(0,max_references)+ ylim(0,max_references)+
  geom_smooth(method="lm",se=FALSE,col="red")+
  geom_abline(slope=1,linetype="dashed",colour="gray52")+
  annotate("text", x=.2*max_references, y=.9*max_references,label=la,parse=TRUE, colour="red",size=2.5)+
  theme(strip.text.x = element_text(size = 9, colour = "black"),
        strip.background = element_blank(),
        legend.position=c(.90,0.75),
        legend.background = element_rect(fill=FALSE),
        legend.text=element_text(size=9),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=9),
        axis.text.x = element_text(colour = "black",size=9),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9))+
  guides(fill=FALSE) +
  labs(x = "Reference Events", y = "Insertion Events")

ggsave(filename="Insertion_vs_Reference.tiff",
       dpi=300,
       width=4,
       height=4,
       units="in")

plot_grid(m1, m2, m3,ncol=3,labels=c('A', 'B','C'))
ggsave(filename="All_vs_All.tiff",
       dpi=300,
       width=7.5,
       height=2.25,
       units="in")

#TRANSPOSONS vs STRAINS
names(summarydata)
#INSERTIONS
insertions<-summarydata[summarydata$method=="new",]
insertions<-(insertions[ order(insertions$total_tes), ])
#plot(insertions$total_tes~insertions$sample)
#pdf(file = "insertions_per_strain.pdf")
m1 <- ggplot(insertions, aes(x=reorder(insertions$sample,insertions$total_tes), y=insertions$total_tes)) 
m1<- m1 + geom_point(size=.75) +aes(group=1)+
  theme(axis.text.x = element_text(color="black",size=8,angle=90,hjust=1),
        axis.text.y = element_text(color="black",size=8),
        axis.title = element_text(color="black",size=9),
        axis.ticks =element_line(colour = "black"))+
  labs(x="", y="Number of Insertion Events")
m1
ggsave(filename="Insertions_per_Strain.tiff",
       dpi=300,
       width=7.5,
       height=10,
       units="in")

#ABSENCES
absences<-summarydata[summarydata$method=="absent",]
absences<-(absences[ order(absences$total_tes), ])
#plot(absences$total_tes~absences$sample)
#pdf(file = "absences_per_strain.pdf")
m2 <- ggplot(absences, aes(x=reorder(absences$sample,absences$total_tes), y=absences$total_tes)) 
m2<- m2 + geom_point(size=.75) +aes(group=1)+
  theme(axis.text.x = element_text(color="black",size=8,angle=90,hjust=1),
        axis.text.y = element_text(color="black",size=8),
        axis.title = element_text(color="black",size=9),
        axis.ticks =element_line(colour = "black"))+
  labs(x="", y="Number of Absence Events")
m2
ggsave(filename="Absences_per_Strain.tiff",
       dpi=300,
       width=7.5,
       height=10,
       units="in")

#REFERENCE
references<-summarydata[summarydata$method=="reference",]
references<-(references[ order(references$total_tes), ])
#plot(references$total_tes~references$sample)
#pdf(file = "references_per_strain.pdf")
m3 <- ggplot(references, aes(x=reorder(references$sample,references$total_tes), y=references$total_tes)) 
m3<- m3 + geom_point(size=.75) +aes(group=1)+
  theme(axis.text.x = element_text(color="black",size=8,angle=90,hjust=1),
        axis.text.y = element_text(color="black",size=8),
        axis.title = element_text(color="black",size=9),
        axis.ticks =element_line(colour = "black"))+
  labs(x="", y="Number of Reference Events")
m3
ggsave(filename="References_per_Strain.tiff",
       dpi=300,
       width=7.5,
       height=10,
       units="in")

plot_grid(m1, m3, m2,ncol=1,labels=c('A', 'B','C'))
ggsave(filename="All_per_Strain.tiff",
       dpi=300,
       width=7.5,
       height=10,
       units="in")

test<-arrange(total_insertion, total_tes)
min(absences$total_tes)
max(absences$total_tes)
mean(absences$total_tes)
min(insertions$total_tes)
max(insertions$total_tes)
mean(insertions$total_tes)
min(references$total_tes)
max(references$total_tes)
mean(references$total_tes)