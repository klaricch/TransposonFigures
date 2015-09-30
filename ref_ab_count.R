#!/usr/bin/R
# this script plots the number of absence vs number of reference calls per strain
# USE: ref_ab_count.R
# NOTE: no longer plotting this

library(ggplot2)
library(grid)

setwd("/Users/kristen/Desktop/Fig_p")
summarydata <- read.table("total_ref_ab_counts.txt")
names(summarydata)<-c("ID","reference_count","absent_count")

#SCATTER
a <- ggplot(data = summarydata, aes(x = reference_count,y=absent_count))
a <- a + geom_point(size=1, position="jitter")+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        strip.text.y = element_text(size = 9, colour = "black",face="bold"),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9),
        legend.position=c(.90,0.85),
        legend.background = element_rect(fill=FALSE),
        legend.key=element_rect(fill=NA),
        legend.text=element_text(size=9))+
  labs(x = "Number of Reference Calls", y="Number of Absence Calls")

a
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="ref_abs_counts_jitter.tiff",dpi=300, width=7.5,height=3.5,units="in")
