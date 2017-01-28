#!/usr/bin/R
# this script
# 1) plots the TELOCATE reference read support vs the TEMP absence read support for all contradictory calls
# 2) plots a histogram of the differences in read support between the TELOCATE reference caller and TEMP absence caller for all contradictory calls
# USE: contradictory_calls.R

library(tidyr)
library(dplyr)
library(ggplot2)
library(gtable)

setwd("/Users/kristen/Documents/transposon_figure_data/data")
summarydata <- read.table("contradictory_calls.txt")
names(summarydata)<-c("strain","chr","start","end","TE","TELOCATE_support","orient", "TE2","TEMP_support","orient2")
summarydata$difference<-abs(summarydata$TEMP_support-summarydata$TELOCATE_support)

#measure degree of difference
summarydata<-mutate(summarydata,rel=((TEMP_support*TELOCATE_support)/(TEMP_support+TELOCATE_support)))

#SCATTER
a <- ggplot(data = summarydata, aes(x = TEMP_support,y=TELOCATE_support,color=ifelse(TEMP_support>TELOCATE_support,"green4","purple4")))+scale_color_identity()
a <- a + geom_point(size=1,position="jitter",alpha=.5)+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        strip.text.y = element_text(size = 9, colour = "black",face="bold"),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=9),
        axis.text.x = element_text(colour = "black",size=9),
        axis.line=element_line(linetype="solid"),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.title=element_text(size=9,face="bold"),
        legend.position=c(.90,0.85),
        legend.background = element_rect(fill=FALSE),
        legend.key=element_rect(fill=NA),
        legend.text=element_text(size=9))+
   scale_y_continuous(expand = c(0,0),limits=c(0,max(summarydata$TELOCATE_support*1.075)))+
  scale_x_continuous(expand = c(0,0),limits=c(0,max(summarydata$TEMP_support*1.075)))+
  labs(x="Number of Reads Supporting Absence Call", y="Number of Reads Supporting Reference Call")
a
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="Contradictory_Calls.tiff",
       dpi=300,
       width=7.5,
       height=3.5,
       units="in")

ggsave(filename="Contradictory_Calls.png",
       dpi=300,
       width=7.5,
       height=3.5,
       units="in")

#HISTOGRAM
(summarydata$difference)
hh <- summarydata[summarydata$difference==0,]

mean(summarydata$difference)
median(summarydata$difference)
m <- ggplot(summarydata, aes(x=difference))
m <- m + geom_histogram(binwidth=2) +
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=9),
        axis.text.x = element_text(colour = "black",size=9),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9))+
  guides(fill=FALSE) +
  labs(x="Difference in Read Support", y="Count")+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))
m
ggsave(filename="Histogram_Contradictory_Calls.tiff",dpi=300, width=7.5,height=3.5,units="in")
ggsave(filename="Histogram_Contradictory_Calls.png",dpi=300, width=7.5,height=3.5,units="in")


#CER1 SCATTER
setwd("/Users/kristen/Documents/transposon_figure_data/data")
comparison<- read.table("cer_comparison.txt",na.strings="NAA")
names(comparison)<-c("strain","TEcaller","paper","outcome")
cer1<-summarydata %>% filter(grepl("CER1_reference", TE))
#add in "CORRECT/INCORRECT calls"
combo <-merge(comparison, cer1, by="strain")


test<-filter(combo, outcome!="NA")
#filtered
hun<-round((nrow(filter(combo, rel>100))/nrow(combo))*100,digits=2)
thirty<-round((nrow(filter(combo, rel>30))/nrow(combo))*100,digits=2)
twenty<-round((nrow(filter(combo, rel>20))/nrow(combo))*100,digits=2)
a <- ggplot(data = combo, aes(x = TEMP_support,y=TELOCATE_support,color=ifelse(rel>100,"pink",ifelse(rel>30,"orange",ifelse(rel>20,"red","black")))))+scale_color_identity()
a <- a + geom_point(size=1,position="jitter",alpha=.75)+
  annotate("text", x=400, y=225, colour="pink",size=3,label=paste(hun,"%",sep=''))+
  annotate("text", x=400, y=200, colour="orange",size=3,label=paste(thirty,"%",sep=''))+
  annotate("text", x=400, y=175, colour="red",size=3,label=paste(twenty,"%",sep=''))+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        strip.text.y = element_text(size = 9, colour = "black",face="bold"),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=9),
        axis.text.x = element_text(colour = "black",size=9),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9,face="bold"),
        legend.position=('none'))+
  labs(x="Number of Reads Supporting Absence Call", y="Number of Reads Supporting Reference Call")
a
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="CER1_Contradictory_Calls_Filtered.tiff",
       dpi=300,
       width=7.5,
       height=3.5,
       units="in")
ggsave(filename="CER1_Contradictory_Calls_Filtered.png",
       dpi=300,
       width=7.5,
       height=3.5,
       units="in")
nrow(filter(combo,outcome=="NA"))
nrow(filter(combo,outcome=="CORRECT"))
nrow(filter(combo,outcome=="INCORRECT"))
#nonfiltered
combo<-arrange(combo,desc(outcome))
a <- ggplot(data = combo, aes(x = TEMP_support,y=TELOCATE_support))
a<- a+ geom_point(aes(colour = outcome))+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        strip.text.y = element_text(size = 9, colour = "black",face="bold"),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=9,face="bold"),
        axis.text.x = element_text(colour = "black",size=9,face="bold"),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9,face="bold"),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        legend.position=('none'))+
  labs(x="Number of Reads Supporting Absence Call", y="Number of Reads Supporting Reference Call")+
  scale_colour_manual(values = c("NA"="grey85","INCORRECT"="darkorange", "CORRECT"="purple3"))
a
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="CER1_Contradictory_Calls.tiff",
       dpi=300,
       width=7.5,
       height=3.5,
       units="in")
ggsave(filename="CER1_Contradictory_Calls.png",
       dpi=300,
       width=7.5,
       height=3.5,
       units="in")
#ifelse(outcome=="CORRECT",'black', ifelse(outcome=="INCORRECT",'pink','green')))

#SCATTER REMOVALS
hun<-round((nrow(filter(summarydata, rel>100))/nrow(summarydata))*100,digits=2)
eighty<-round((nrow(filter(summarydata, rel>80))/nrow(summarydata))*100,digits=2)
thirty<-round((nrow(filter(summarydata, rel>30))/nrow(summarydata))*100,digits=2)
twenty<-round((nrow(filter(summarydata, rel>20))/nrow(summarydata))*100,digits=2)
a <- ggplot(data = summarydata, aes(x = TEMP_support,y=TELOCATE_support,color=ifelse(rel>100,"pink",ifelse(rel>80,"green",ifelse(rel>30,"orange",ifelse(rel>20, "red","black"))))))+scale_color_identity()
a <- a + geom_point(size=1,position="jitter",alpha=.75)+
  annotate("text", x=1500, y=250, colour="pink",size=3,label=paste(hun,"%",sep=''))+
  annotate("text", x=1500, y=225, colour="green",size=3,label=paste(eighty,"%",sep=''))+
  annotate("text", x=1500, y=200, colour="orange",size=3,label=paste(thirty,"%",sep=''))+
  annotate("text", x=1500, y=175, colour="red",size=3,label=paste(twenty,"%",sep=''))+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        strip.text.y = element_text(size = 9, colour = "black",face="bold"),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=9,face="bold"),
        axis.text.x = element_text(colour = "black",size=9,face="bold"),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9),
        legend.position=c(.90,0.85),
        legend.background = element_rect(fill=FALSE),
        legend.key=element_rect(fill=NA),
        legend.text=element_text(size=9))+
  labs(x="Number of Reads Supporting Absence Call", y="Number of Reads Supporting Reference Call")
a

setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="Contradictory_Calls_Removals.tiff",
       dpi=300,
       width=7.5,
       height=3.5,
       units="in")
ggsave(filename="Contradictory_Calls_Removals.png",
       dpi=300,
       width=7.5,
       height=3.5,
       units="in")