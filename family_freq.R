#!/usr/bin/R
# this script graphs the total absences, insertions, and references per transposon family                   
# USE: family_freq.R

library(dplyr)
library(stringr)
library(ggplot2)
library(grid)
library(stringr)

setwd("/Users/kristen/Documents/transposon_figure_data")
summarydata <- read.table("T_Full_Results.txt",header=TRUE)
classdata<- read.table("CtCp_all_nonredundant.txt",header=TRUE)
names(classdata)<-c("chr","start","end","TE","support","orientation","method","strain","class")

# add te class info to summarydata(new_TRANS_end_tes will be removed)
classdata$id<- stringr::str_split_fixed(classdata$TE, regex("_(non-)?reference"),2)[,1]
classdata<-mutate(classdata, trait=paste(method,"TRANS",id,sep="_"))
class_subset <- classdata %>% distinct(trait) %>% select(trait,class)
summarydata <-merge(summarydata, class_subset, by="trait")

#revalue classes
summarydata$class <- factor(summarydata$class,
                            levels = c("dnatransposon", "retrotransposon","unknown"),
                            labels = c("DNA Transposon", "Retrotransposon", "Unknown"))

#remove total counts and coverage
summarydata <- filter(summarydata, trait != "absent_TRANS_total" )
summarydata <- filter(summarydata, trait != "new_TRANS_total" )
summarydata <- filter(summarydata, trait != "reference_TRANS_total" )
summarydata <- filter(summarydata, trait != "coverage" )
summarydata<-mutate(summarydata, TOTAL=rowSums(summarydata[2:125]))

#new column that specifies what caller was used
summarydata$caller<- stringr::str_split_fixed(summarydata$trait, "_TRANS_",2)[,1]
#new column that specifies TE family
summarydata$transposon<- stringr::str_split_fixed(summarydata$trait, "_TRANS_",2)[,2]

#revalue
summarydata$caller<- factor(summarydata$caller,
                          levels = c("absent", "new","reference"),
                          labels = c("Absence", "Insertion", "Reference"))

m <- ggplot(summarydata,aes(y=transposon,x=TOTAL))
m <- m + geom_point(size=1.25,aes(color=class))+
  facet_grid(.~caller, scale="free") +
  theme(strip.background = element_rect(fill="white"),
        strip.text.y = element_text(size = 8, colour = "black",face="bold"),
        panel.margin = unit(.6, "lines"),
        panel.border = element_rect(fill=NA,colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey87"),
        panel.grid.minor = element_line(colour = "grey87"),
        axis.ticks =element_line(colour = "black"),
        axis.title=element_text(size=9),
        axis.text.y = element_text(colour = "black",size=5),
        axis.text.x = element_text(color="black",size=8),
        legend.title=element_blank(),
        legend.background = element_rect(fill=FALSE),
        legend.key=element_rect(fill=NA),
        legend.text=element_text(size=9))+
  scale_color_manual(values = c("navy", "brown3", "darkgoldenrod2"))+
  labs(y="Transposon", x="Total Transposition Events")
m
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="Family_Frequency.tiff",
       dpi=300,
       width=7.5,
       height=10,
       units="in")

