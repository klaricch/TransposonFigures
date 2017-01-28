#!/usr/bin/R
# this script graphs the total absences, insertions, and references per transposon family                   
# USE: family_freq.R

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(grid)
library(stringr)

setwd("/Users/kristen/Documents/transposon_figure_data/data")
summarydata <- read.table("T_kin_C_matrix_full.txt",header=TRUE)
#remove ZERO_new traits
summarydata<-subset(summarydata,!grepl('^ZERO_new', summarydata$trait))
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
class_subset<-distinct(classdata, family,class,.keep_all=TRUE)
class_subset<-dplyr::select(class_subset,family,class)
#class_subset <- classdata %>% distinct(family,class) %>% select(family,class)
summarydata$family<- stringr::str_split_fixed(summarydata$trait, "_TRANS_",2)[,2]
summarydata <-merge(summarydata, class_subset, by="family")
summarydata<-dplyr::select(summarydata, -family)
#revalue classes
summarydata$class <- factor(summarydata$class,
                            levels = c("dnatransposon", "retrotransposon","unknown"),
                            labels = c("DNA Transposon", "Retrotransposon", "Unknown"))

#double check removed total, but shouls have been removed in the merge
summarydata<-filter(summarydata,!grepl('total', trait))
summarydata <- filter(summarydata, trait != "coverage" )

no_cols<-ncol(summarydata)-1
print(summarydata[,no_cols])
summarydata<-mutate(summarydata, TOTAL=rowSums(summarydata[2:no_cols],na.rm = TRUE))

#new column that specifies what caller was used
summarydata$caller<- stringr::str_split_fixed(summarydata$trait, "_TRANS_",2)[,1]
#new column that specifies TE family
summarydata$transposon<- stringr::str_split_fixed(summarydata$trait, "_TRANS_",2)[,2]

#revalue
summarydata$caller<- factor(summarydata$caller,
                          levels = c("new","reference","absent"),
                          labels = c("Insertion","Reference","Absence"))


summarydata$transposon<-gsub("_CE$","",summarydata$transposon)
summarydata$transposon<-gsub("WBTransposon","WBT",summarydata$transposon)

#colnames(summarydata)
summarydata<-dplyr::select(summarydata,trait,class,TOTAL,caller,transposon)

########
summarydata<-summarydata
summarydata<-dplyr::select(summarydata,-trait)
summarydata<-spread(summarydata,caller,TOTAL)
summarydata$Reference[is.na(summarydata$Reference)] <- 0
summarydata<-mutate(summarydata,All=Insertion+Reference)
summarydata<-dplyr::select(summarydata,-Reference)
summarydata<-gather(summarydata, "caller","TOTAL",Insertion,Absence,All)

####
summarydata$caller = factor(summarydata$caller, levels=c('Insertion','Absence','All'),
                            labels = c("Insertion Sites", "Active Reference Sites","All Transposon Sites"))

m <- ggplot(summarydata,aes(y=transposon,x=TOTAL))
m <- m + geom_point(size=1.25,aes(color=class))+
  facet_grid(.~caller, scale="free") +
  theme(strip.background = element_rect(fill="white"),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        panel.spacing = unit(.6, "lines"),
        panel.spacing.y=unit(.50,"cm"),
        plot.margin=unit(c(.1,.1,0,.1), "cm"),
        panel.border = element_rect(fill=NA,colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey87"),
        panel.grid.minor = element_line(colour = "grey87"),
        axis.ticks =element_line(colour = "black"),
        axis.title=element_text(size=9),
        axis.title.x=element_text(face="bold"),
        axis.text.y = element_text(colour = "black",size=5),
        axis.text.x = element_text(color="black",size=8),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        legend.title=element_blank(),
        legend.background = element_rect(fill=FALSE),
        legend.key=element_rect(fill=NA),
        legend.position="none",
        legend.text=element_text(size=9))+
  scale_color_manual(values = c("DNA Transposon" = "navy", "Retrotransposon"="brown3","Unknown"="darkgoldenrod2"))+
  labs(y="", x="Total Sites")
m
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="Family_Frequency.tiff", dpi=300, width=7.5, height=10, units="in")
ggsave(filename="Family_Frequency.png", dpi=300, width=7.5, height=10, units="in")


setwd("/Users/kristen/Documents/transposon_figure_data/data")

total_means<-summarydata %>% group_by(caller,class) %>% summarise(mean=mean(TOTAL,na.rm=TRUE),SD=sd(TOTAL, na.rm=TRUE))
total_means<- mutate(total_means, id = paste(class,caller,sep="_"))
summarydata<- mutate(summarydata, id = paste(class,caller,sep="_"))

merged<-merge(summarydata, total_means, by="id")
merged<-mutate(merged, outL=ifelse(abs(TOTAL-mean)>SD, "OUTLIER", "n"))
outliers<-filter(merged,outL=="OUTLIER")
save(outliers,file="outlier_total_events_per_family.Rda")


setwd("/Users/kristen/Documents/transposon_figure_data/figures")
outlier_table<-dplyr::select(outliers,transposon, caller.x,class.x,TOTAL,mean,SD)
outlier_table$mean<-signif(outlier_table$mean,4)
outlier_table$SD<-signif(outlier_table$SD,4)
outlier_table<-arrange(outlier_table,transposon,TOTAL,caller.x)
colnames(outlier_table)<-c("Transposon", "Site Type","Class","Total Number","Mean","SD")
write.table(outlier_table, file="Outlier_Table.txt",sep="\t",quote=FALSE,row.names=FALSE)

names(summarydata)
test<-distinct(summarydata,caller, transposon)
RR<-filter(summarydata,caller=="Reference")
AA<-filter(summarydata,caller=="Active References")
NN<-filter(summarydata,caller=="Insertions")
CC<-filter(summarydata,caller=="Reference"|caller=="Active References")
length(unique(summarydata$transposon))
(unique(summarydata$caller))
length(unique(RR$transposon))
length(unique(AA$transposon))
length(unique(NN$transposon))
length(unique(CC$transposon))
test<-distinct(CC,caller, transposon)

length(unique(CC$trait))


test<-filter(merged,transposon=="MARINER2")
test<-filter(merged,caller.y=="Insertions")


