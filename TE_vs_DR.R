#!/usr/bin/R
# this script
# 1) plots total transposons vs strain per insertions, references, and absences per class
# 2) plot histogram of transposons per strain
# USE: TE_vs_DR.R
# NOTE: y max values are hard coded

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(grid)
library(MASS)
library(cowplot)
select <- dplyr::select
setwd("/Users/kristen/Documents/transposon_figure_data/data")
summarydata <- read.table("T_kin_C_matrix_full.txt",header=TRUE)
#remove ZERO_new traits
summarydata<-subset(summarydata,!grepl('^ZERO_new', summarydata$trait))
summarydata<-subset(summarydata,!grepl('^coverage', summarydata$trait))
#clean trait names
summarydata$trait <- gsub("_C$" ,"",summarydata$trait)
summarydata$trait <- gsub("^ONE_new" ,"new",summarydata$trait)

classdata<- read.table("CtCp_all_nonredundant.txt")
names(classdata)<-c("chr","start","end","TE","orientation","method","strain","class")

# add te class info to summarydata(new_TRANS_end_tes will be removed)
classdata$family<- stringr::str_split_fixed(classdata$TE, regex("_(non-)?reference"),2)[,1]
classdata$family<- paste(stringr::str_split_fixed(classdata$family, "_",4)[,3],stringr::str_split_fixed(classdata$family, "_",4)[,4],sep="_")
classdata$family <- gsub("_$" ,"",classdata$family)
classdata$family <- gsub("_non-reference(.*)$" ,"",classdata$family)
classdata<-mutate(classdata, trait=paste(method,"TRANS",family,sep="_"))
class_subset <- classdata %>% distinct(family,class,.keep_all=TRUE) %>% select(family,class)
summarydata$family<- stringr::str_split_fixed(summarydata$trait, "_TRANS_",2)[,2]
summarydata <-merge(summarydata, class_subset, by="family")
summarydata<-select(summarydata, -family)
unique(class_subset$method)


#summarydata <-merge(summarydata, class_subset, by="trait")

#nrow(distinct(classdata,method,TE))
#x<-merge(summarydata, class_subset, by="family")
#Z<-filter(summarydata, !(summarydata$trait %in% x$trait))
#unique(summarydata$trait)
#tail(classdata)
#unique(classdata$family)
#test1<-filter(classdata, family=="CEMUDR2")
#names(summarydata)

summarydata$trait
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
  stat_smooth(se=FALSE,colour="red",method=function(formula,data,weights=weight) rlm(formula,
                                                                                     data,
                                                                                     weights=weight,
                                                                                     method="MM"),fullrange=TRUE)+
  geom_abline(slope=1,linetype="dashed",colour="gray52")+
  geom_text(data = subset(final_merge, class=="DNA Transposon"), x=.15*D_max_insertions, y=.9*D_max_insertions,label=Dla,parse=TRUE, colour="red",size=4)+
  geom_text(data = subset(final_merge, class=="Retrotransposon"), x=.15*R_max_insertions, y=.9*R_max_insertions,label=Rla,parse=TRUE, colour="red",size=4)+
  geom_point(data = subset(final_merge, class=="Retrotransposon"), aes(x=R_max_insertions,y=R_max_insertions),alpha=0)+
  geom_point(data = subset(final_merge, class=="DNA Transposon"), aes(x=D_max_insertions,y=D_max_insertions),alpha=0)+
  geom_point(aes(x=0,y=0),alpha=0)+
  
  theme(strip.text.x = element_text(size = 11, colour = "black",face="bold"),
        strip.background = element_blank(),
        legend.position=c(.90,0.75),
        legend.background = element_rect(fill=FALSE),
        legend.text=element_text(size=11),
        panel.background = element_rect(fill = "white"),
        panel.spacing= unit(1, "lines"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=11),
        axis.text.x = element_text(colour = "black",size=11),
        axis.line=element_line(linetype="solid"),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.title=element_text(size=11,face="bold"))+
  guides(fill=FALSE) +
  labs(x = "Insertion Sites", y = "Active Reference Sites")
m
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="Absence_vs_Insertion_DR.tiff", dpi=350, width=6, height=4, units="in")
ggsave(filename="Absence_vs_Insertion_DR.png", dpi=350, width=7.5, height=3.5, units="in")
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
  stat_smooth(se=FALSE,colour="red",method=function(formula,data,weights=weight) rlm(formula,
                                                                                     data,
                                                                                     weights=weight,
                                                                                     method="MM"),fullrange=TRUE)+
  geom_abline(slope=1,linetype="dashed",colour="gray52")+
  geom_text(data = subset(final_merge, class=="DNA Transposon"), x=.1*D_max_references, y=.9*D_max_references,label=Dla,parse=TRUE, colour="red",size=4)+
  geom_text(data = subset(final_merge, class=="Retrotransposon"), x=.1*R_max_references, y=.9*R_max_references,label=Rla,parse=TRUE, colour="red",size=4)+
  geom_point(data = subset(final_merge, class=="Retrotransposon"), aes(x=R_max_references,y=R_max_references),alpha=0)+
  geom_point(data = subset(final_merge, class=="DNA Transposon"), aes(x=D_max_references,y=D_max_references),alpha=0)+
  geom_point(aes(x=0,y=0),alpha=0)+
  
  theme(strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        strip.background = element_blank(),
        legend.position=c(.90,0.75),
        legend.background = element_rect(fill=FALSE),
        legend.text=element_text(size=9),
        panel.background = element_rect(fill = "white"),
        panel.spacing = unit(1, "lines"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=9),
        axis.text.x = element_text(colour = "black",size=9),
        axis.line=element_line(linetype="solid"),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.title=element_text(size=9,face="bold"))+
  guides(fill=FALSE) +
  labs(x = "Reference Sites", y = "Active Reference Sites")
m
ggsave(filename="Absence_vs_Reference_DR.tiff", dpi=300, width=7.5, height=3.5, units="in")
ggsave(filename="Absence_vs_Reference_DR.png", dpi=300, width=7.5, height=3.5, units="in")

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
  stat_smooth(se=FALSE,colour="red",method=function(formula,data,weights=weight) rlm(formula,
                                                                                     data,
                                                                                     weights=weight,
                                                                                     method="MM"),fullrange=TRUE)+
  geom_abline(slope=1,linetype="dashed",colour="gray52")+
  geom_text(data = subset(final_merge, class=="DNA Transposon"), x=.1*D_max_insertions, y=.9*D_max_insertions,label=Dla,parse=TRUE, colour="red",size=4)+
  geom_text(data = subset(final_merge, class=="Retrotransposon"), x=.1*R_max_insertions, y=.9*R_max_references,label=Rla,parse=TRUE, colour="red",size=4)+
  geom_point(data = subset(final_merge, class=="Retrotransposon"), aes(x=R_max_references,y=R_max_references),alpha=0)+
  geom_point(data = subset(final_merge, class=="DNA Transposon"), aes(x=D_max_insertions,y=D_max_insertions),alpha=0)+
  geom_point(aes(x=0,y=0),alpha=0)+
  
  theme(strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        strip.background = element_blank(),
        legend.position=c(.90,0.75),
        legend.background = element_rect(fill=FALSE),
        legend.text=element_text(size=9),
        panel.background = element_rect(fill = "white"),
        panel.spacing = unit(1, "lines"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=9),
        axis.text.x = element_text(colour = "black",size=9),
        axis.line=element_line(linetype="solid"),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.title=element_text(size=9,face="bold"))+
  guides(fill=FALSE) +
  labs(x = "Insertion Sites", y = "Reference Sites")
m
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="Insertion_vs_Reference_DR.tiff", dpi=300, width=7.5, height=3.5, units="in")
ggsave(filename="Insertion_vs_Reference_DR.png", dpi=300, width=7.5, height=3.5, units="in")


#CHEXK T HIS  OVER
# NO LONGER NEED THE REST OF THE BELOW
#FAMILY
sl<-select(families, -trait)
sl<-spread(sl,method,value)

sl<-rename(sl, New = new) #Rename the column "new"
AR<-filter(sl, absent!="NA" & reference !="NA")
AI<-filter(sl, absent!="NA" & New !="NA")
RI<-filter(sl, reference!="NA" & New !="NA")

#calculate rho values
AR_corr<-AR %>% group_by(class,transposon) %>% summarise(cor=round((cor.test(absent,reference,method="spearman",exact=FALSE))$estimate,3))
AI_corr<-AI %>% group_by(class,transposon) %>% summarise(cor=round((cor.test(absent,New,method="spearman",exact=FALSE))$estimate,3))
RI_corr<-RI %>% group_by(class,transposon) %>% summarise(cor=round((cor.test(reference,New,method="spearman",exact=FALSE))$estimate,3))

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

all<-mutate(all, outL=ifelse(abs(cor-mean)>(2*SD), "OUTLIER", "n"))
outliers<-filter(all,outL=="OUTLIER")
setwd("/Users/kristen/Documents/transposon_figure_data/data")
save(outliers,file="outlier_slopes.Rda")






#Histogram Versions



test<-hist_data
hist_data<-spread(hist_data,method,XX)
hist_data<-mutate(hist_data, all=reference+new)
hist_data<-select(hist_data,new,absent,all)
hist_data<-gather(hist_data, "method","value",new,absent,all)
hist_data$method <- factor(hist_data$method,
                       levels = c("new", "absent", "all"),
                       labels = c("Insertion Sites", "Active Reference Sites", "All Transposon Sites"))


#dim(hist_data)
#sub_hist<-hist_data %>% group_by(method) %>% summarize(MM=max(XX))
#hist_data<-merge(hist_data,sub_hist,by="method")

unique(hist_data$method)
ins_data<-filter(hist_data, method=="Insertion Sites")
m <- ggplot(ins_data, aes(x=value,fill=class))
m <-m + geom_histogram(binwidth=5)+
  scale_y_continuous(expand = c(0,0),limits=c(0,180)) + scale_x_continuous(expand = c(0,0))+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 11, colour = "black",face="bold"),
        panel.spacing = unit(.25, "lines"),
        panel.border = element_rect(fill=NA,colour = "black"),
        panel.background = element_blank(),
        axis.title=element_text(size=11,face="bold"),
        axis.text.y = element_text(colour = "black",size=11),
        axis.text.x = element_text(colour = "black",size=11),
        axis.ticks =element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        legend.title=element_blank(),
        legend.position="none")+
  scale_fill_manual(values = c('dnatransposon' = "navy", "retrotransposon"="brown3","unknown"="goldenrod"))+
  labs(x="Insertion Sites Per Strain", y="Count")
m
  

ar_data<-filter(hist_data, method=="Active Reference Sites")
m2 <- ggplot(ar_data, aes(x=value,fill=class))
m2 <-m2 + geom_histogram(binwidth=5)+
  scale_y_continuous(expand = c(0,0),limits=c(0,180)) + scale_x_continuous(expand = c(0,0))+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 11, colour = "black",face="bold"),
        panel.spacing = unit(.25, "lines"),
        panel.border = element_rect(fill=NA,colour = "black"),
        panel.background = element_blank(),
        axis.title=element_text(size=11,face="bold"),
        axis.text.y = element_text(colour = "black",size=11),
        axis.text.x = element_text(colour = "black",size=11),
        axis.ticks =element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        legend.title=element_blank(),
        legend.position="none")+
  scale_fill_manual(values = c('dnatransposon' = "navy", "retrotransposon"="brown3","unknown"="goldenrod"))+
  labs(x="Active Reference Sites Per Strain", y="Count")
m2
#m2<- m2 + geom_point(data = subset(ar_data, method=="Active Reference Sites"),aes(y=160),alpha=0)+
##pull out y max of panels
#panel1<-filter(ggplot_build(m)$data[[1]],PANEL==1)
#max1<-max(panel1$y)
#max1
#panel2<-filter(ggplot_build(m)$data[[1]],PANEL==2)
#max2<-max(panel2$y)
#max2
#panel3<-filter(ggplot_build(m)$data[[1]],PANEL==3)
#max3<-max(panel3$y)
#max3

#geom_point(data = subset(hist_data, method=="Active Reference Sites"),aes(y=1.05*max2),alpha=0)
 # geom_point(data = subset(hist_data, method=="All Transposon Sites"),aes(y=1.05*max3),alpha=0)
m




#allele freq

summarydata <- read.table("kin_matrix_full.txt",header=TRUE)
#count number of "1" scores
no_cols=ncol(summarydata)
summarydata<-mutate(summarydata, freq=rowSums(summarydata[2:no_cols],na.rm = TRUE))
#count number of non-NA value
summarydata<-mutate(summarydata,total=rowSums(!is.na(summarydata[2:no_cols])))
#calcualte allele frequency
#summarydata<-mutate(summarydata,AF_pre=round((freq/total),digits=3))
summarydata<-mutate(summarydata,AF_pre=(freq/total))

summarydata<-mutate(summarydata,call=ifelse(grepl('NR',trait),"Novel Site", "Reference Site"))
summarydata<-mutate(summarydata,AF=ifelse(AF_pre>=.50, abs(AF_pre-1), AF_pre))


sort(unique(summarydata$AF))
sort(unique(summarydata$AF))
summarydata$AF
summarydata<-filter(summarydata, AF!=0.000) #double check
summarydata<-select(summarydata, AF)

a <- ggplot(summarydata, aes(x=AF))
a <- a + geom_histogram(binwidth=.05) +
  # facet_grid(call ~ .,scale="free_y")+
  theme(panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill="white"),
        axis.ticks = element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=11),
        axis.text.x = element_text(colour = "black",size=11),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=11,face="bold"))+
  guides(fill=FALSE) +
  labs(x="Allele Frequency", y="Count")+
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0))
a


a_all<-plot_grid(m,m2,a,ncol=3,labels=c("A","B","C"))
a_all
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="Histogram_TE_per_Strain_AF.tiff", dpi=350, width=6.75, height=3.375, units="in")
ggsave(filename="Histogram_TE_per_Strain_AF.png", dpi=300, width=7.5, height=3.5, units="in")

#ggsave(filename="Histogram_TE_per_Strain.tiff", dpi=300, width=7.5, height=3.5, units="in")
#ggsave(filename="Histogram_TE_per_Strain.png", dpi=300, width=7.5, height=3.5, units="in")

#alternative
alt_fig<-filter(hist_data, method=="All Transposon Sites")
m3 <- ggplot(alt_fig, aes(x=value,fill=class))
m3 <-m3 + geom_histogram(binwidth=5)+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))+
  #facet_wrap(~method,scale="free")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 9, colour = "black",face="bold"),
        panel.spacing = unit(.25, "lines"),
        panel.border = element_rect(fill=NA,colour = "black"),
        panel.background = element_blank(),
        axis.title=element_text(size=9,face="bold"),
        axis.text.y = element_text(colour = "black",size=9),
        axis.text.x = element_text(colour = "black",size=9),
        axis.ticks =element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        legend.title=element_blank(),
        legend.position="none")+
  scale_fill_manual(values = c('dnatransposon' = "navy", "retrotransposon"="brown3","unknown"="goldenrod"))+
  labs(x="All Transposon Sites Per Strain", y="Count")
m3

ggsave(filename="Histogram_TE_per_Strain_All.png", dpi=300, width=7.5, height=3.5, units="in")


b_all<-plot_grid(m3,a,ncol=2,labels=c("A","B"))
b_all
ggsave(b_all,filename="Histogram_TE_per_Strain_All_b.png", dpi=300, width=7.5, height=3.5, units="in")


b_supp<-plot_grid(m,m2,ncol=2,labels=c("A","B"))
b_supp
ggsave(b_supp,filename="Histogram_TE_per_Strain_All_supp3b.png", dpi=300, width=7.5, height=3.5, units="in")
