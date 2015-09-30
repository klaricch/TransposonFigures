#!/usr/bin/R
# this script
# 1) plots the TELOCATE reference read support vs the TEMP absence read support for all contradictory calls
# 2) plots a histogram of the differences in read support between the TELOCATE reference caller and TEMP absence caller for all contradictory calls
# USE: contradictory_calls.R

library(tidyr)
library(dplyr)
library(ggplot2)
library(gtable)

##REMOVE BELOW LATER
setwd("/Users/kristen/Documents/transposon_figure_data")
summarydata <- read.table("contradictory_calls.txt")
names(summarydata)<-c("strain","chr","start","end","TE","TELOCATE_support","orient", "TE2","TEMP_support","orient2")
summarydata$difference<-abs(summarydata$TEMP_support-summarydata$TELOCATE_support)

#SCATTER
a <- ggplot(data = summarydata, aes(x = TEMP_support,y=TELOCATE_support,color=ifelse(TEMP_support>TELOCATE_support,"green4","purple4")))+scale_color_identity()
a <- a + geom_point(size=1,position="jitter")+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        strip.text.y = element_text(size = 9, colour = "black",face="bold"),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=9),
        axis.text.x = element_text(colour = "black",size=9),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9),
        legend.position=c(.90,0.85),
        legend.background = element_rect(fill=FALSE),
        legend.key=element_rect(fill=NA),
        legend.text=element_text(size=9))+
  labs(x="Absence Call Read Support", y="Reference Call Read Support")
a
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="Contradictory_Calls.tiff",
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
m <- m + geom_bar(binwidth=2) +
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
