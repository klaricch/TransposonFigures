#!/usr/bin/R
# this script plots a histogram of percentage of NAs per TE position

library(dplyr)
library(tidyr)
library(ggplot2)
setwd("/Users/kristen/Documents/transposon_figure_data/data")
nas<-read.table("NA_counts_at_positions.txt")
colnames(nas)<-c("trait", "NA_count","fraction")

#m <- geom_histogram(nas)
m <- ggplot(nas, aes(x=fraction))
m <- m + geom_histogram(binwidth=.05) 
m <- m +theme(
  panel.background = element_rect(fill = "white"),
       # strip.background = element_rect(fill="white"),
        axis.ticks = element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=9),
        axis.text.x = element_text(colour = "black",size=9),
        axis.line.x=element_line(linetype="solid",colour="black"),
        axis.line.y=element_line(linetype="solid",colour="black"),
        #axis.line=element_line(linetype="solid",colour="black"),
        axis.title=element_text(size=9,face="bold"))+ 
  geom_vline(aes(xintercept=.1),color="grey60",linetype="dashed")+
  scale_x_continuous(expand = c(0,0))

max_y<-ggplot_build(m)$panel$ranges[[1]]$y.range
max_y<-max_y[2]
m<- m + scale_y_continuous(expand = c(0,0),limits=c(0,max_y*1.075))
m
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="NAs_Per_Pos.tiff",dpi=300, width=7.5,height=3.5,units="in")
ggsave(filename="NAs_Per_Pos.png",dpi=300, width=7.5,height=3.5,units="in")

test<-filter(nas, fraction>.1)



######################################################
######################################################
######################################################

