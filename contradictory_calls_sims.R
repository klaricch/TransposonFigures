#!/usr/bin/R
# this script
# 1) plots the TELOCATE reference read support vs the TEMP absence read support for all contradictory calls in the simulations
# 2) plots a histogram of the differences in read support between the TELOCATE reference caller and TEMP absence caller for all contradictory calls in the simulations
# NOTE: this script is for the simulation data
# USE: contradictory_calls.R

library(tidyr)
library(dplyr)
library(ggplot2)
library(gtable)
library(stringr)
library(cowplot)
library(gridExtra)


setwd("/Users/kristen/Documents/transposon_figure_data/contradictory_calls_RSV_sims")

file_list=c("RSV_simulations_Aug31_half_contra_call_TPFD.txt",
            "RSV_simulations_Aug31_quarter_contra_call_TPFD.txt",
            "RSV_simulations_Aug31_tenth_contra_call_TPFD.txt",
            "RSV_simulations_Aug31_full_contra_call_TPFD.txt")

# loop through the 4 coverage levels of the simulations
for (i in file_list){
  print(i)
summarydata <- read.table(i)

#set column names
names(summarydata)<-c("strain","chr","start","end","TE","TEMP_support","orient", "TE2","TELOCATE_support","orient2", "TEMP", "TELOCATE")
summarydata$difference<-abs(summarydata$TEMP_support-summarydata$TELOCATE_support)

file<-i
fileID<-str_match(file, 'RSV_simulations_Aug31_(.*)_contra_call_TPFD.txt')[2]

#set title names
if (fileID=="full"){ID<-"130x Depth of Coverage"}
if (fileID=="half"){ID<-"65x Depth of Coverage"}
if (fileID=="quarter"){ID<-"32.5x Depth of Coverage"}
if (fileID=="tenth"){ID<-"13x Depth of Coverage"}

#assign call of TP or FD corresponding to the call with most read support
summarydata<-mutate(summarydata, TPFD=ifelse(TELOCATE_support>TEMP_support,as.character(TELOCATE),as.character(TEMP))) #use as.character otehrwise interpret string as integer

#SCATTER
a <- ggplot(data = summarydata, aes(x = TEMP_support,y=TELOCATE_support))
a <- a + geom_point(size=1.5,position="jitter", aes(color=ifelse(TPFD=="TP","slateblue1", "darkorange")))+scale_color_identity()+
  #geom_smooth(method=lm,se=FALSE)+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 8, colour = "black",face="bold"),
        strip.text.y = element_text(size = 8, colour = "black",face="bold"),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=8),
        axis.text.x = element_text(colour = "black",size=8),
        axis.title = element_text(colour = "black",size=8),
        axis.line=element_line(linetype="solid"),
        plot.title = element_text(size=9),
        legend.text=element_text(size=8))+
  labs(x="Absence Call Read Support", y="Reference Call Read Support", title=ID)

fileID
filename <- paste(ID,"contradictory_calls.tiff", sep ="_")
ggsave(filename,
       dpi=300,
       width=7.5,
       height=3.5,
       units="in")

#HISTOGRAM
(summarydata$difference)
hh <- summarydata[summarydata$difference==0,]
print(hh)

mean(summarydata$difference)
median(summarydata$difference)
m <- ggplot(summarydata, aes(x=difference))
m <- m + geom_bar(binwidth=1) +
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=8),
        axis.text.x = element_text(colour = "black",size=8),
        axis.title = element_text(colour = "black",size=8),
        axis.line=element_line(linetype="solid"),
        plot.title = element_text(size=9),
        axis.title=element_text(size=8))+
  guides(fill=FALSE) +
  labs(x="Difference in Read Support", y="Count", title=ID)+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))
m
filename <- paste(ID,"histogram_contradictory_calls.tiff", sep ="_")
ggsave(filename,dpi=300, width=7.5,height=3.5,units="in")


if (fileID=="full"){first<-a;m_first<-m}
if (fileID=="half"){second<-a;m_second<-m}
if (fileID=="quarter"){third<-a;m_third<-m}
if (fileID=="tenth"){fourth<-a;m_fourth<-m}

}

a_all<-plot_grid(first,second,third,fourth)
m_all<-plot_grid(m_first,m_second,m_third,m_fourth)

setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(a_all,filename="combined_contradictory_scatter.tiff",dpi=300, width=7.5,height=5,units="in")
ggsave(m_all,filename="combined_contradictory_histogram.tiff",dpi=300, width=7.5,height=5,units="in")
