#!/usr/bin/R
# this script:
# 1) plots total TE counts vs coverage
# 2) plots coverage levels per strain
# USE: coverage.R

setwd("/Users/kristen/Documents/transposon_figure_data/data")
summarydata <- read.table("coverage_and_te_counts.txt",header=TRUE)
library(tidyr)
library(ggplot2)
library(dplyr)

#data<-summarydata %>%
#  gather(method, total_tes, absence:insertion:reference)

data<-summarydata %>%
  gather(method, total_tes, absence:reference)

#3X_REFACET
a <- ggplot(data = data, aes(x = coverage,y=total_tes, fill=method))
#alternative for color setting below
a <- a + geom_point(size=1,aes( color=ifelse(method=="reference", 'slateblue1',ifelse(method=="insertion",'turquoise3','darkorange')))) +scale_color_identity()+
#a <- a + geom_point(size=1) +scale_color_identity()+
  geom_smooth(method=lm,se=FALSE,aes( color=ifelse(method=="reference", 'slateblue1',ifelse(method=="insertion",'turquoise3','darkorange'))))+
  facet_wrap(~method,ncol=3,scale="free_y")+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        strip.text.y = element_text(size = 9, colour = "black",face="bold"),
        legend.position=('none'),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9))+
  labs(x = "Depth of Coverage", y="Number of Sites")
a
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="Totals_vs_Coverage.tiff",dpi=300, width=7.5,height=3.5,units="in")
ggsave(filename="Totals_vs_Coverage.png",dpi=300, width=7.5,height=3.5,units="in")


########################################################################################
########################################################################################
########################################################################################

# Coverage per Strain
names(data)
m <- ggplot(data=data,aes(x=reorder(sample,coverage),y=coverage))
m <- m + geom_point(size=.75)+
  theme(axis.text.x =element_text(angle = 90, hjust = 1, color = "black", size =5),
        axis.text.y=element_text(angle = 90, hjust = 1, color = "black", size = 9),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.ticks=element_line(colour = "black"))+
  labs(x="Strain",y="Depth of Coverage")
m
ggsave(filename="Coverage_per_Strain.tiff", dpi=300, width=7.5, height=3.5, units="in")
ggsave(filename="Coverage_per_Strain.png", dpi=300, width=7.5, height=3.5, units="in")

mean(data$coverage)
median(data$coverage)
