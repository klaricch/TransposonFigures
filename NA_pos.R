#!/usr/bin/R
# this script plots a histogram of percentage of NAs per TE position

library(dplyr)
library(tidyr)
library(ggplot2)
setwd("/Users/kristen/Documents/transposon_figure_data/data")
nas<-read.table("NA_counts_at_positions.txt")
colnames(nas)<-c("trait", "NA_count","fraction")
#str(nas)
nas<-mutate(nas,fraction=NA_count/152)
#nas<-mutate(nas,fraction=NA_count)
m <- geom_histogram(nas)
m <- ggplot(nas, aes(x=fraction))
m <- m + geom_histogram(binwidth=.05) 
m <- m +theme(panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill="white"),
        axis.ticks = element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=9),
        axis.text.x = element_text(colour = "black",size=9),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9,face="bold"))

m
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="NAs_Per_Pos.tiff",dpi=300, width=7.5,height=3.5,units="in")

test<-filter(nas, fraction>.1)



######################################################
######################################################
######################################################

