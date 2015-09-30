#!/usr/bin/R
# this script graphs the total absences, insertions, and references per transposon family                   
# USE: family_freq.R

library(dplyr)
library(stringr)
library(ggplot2)
library(grid)

setwd("/Users/kristen/Documents/transposon_figure_data")
summarydata <- read.table("T_Full_Results.txt",header=TRUE)
names(summarydata)

ncol(summarydata)
nrow(summarydata)
#remove total coutns and coverage
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

m <- ggplot(summarydata,aes(x=transposon,y=TOTAL))
m <- m + geom_point(size=1)+
  facet_grid(caller~., scale="free_y") +
  theme(strip.background = element_rect(fill="white"),
        strip.text.y = element_text(size = 8, colour = "black",face="bold"),
        panel.margin = unit(.6, "lines"),
        panel.border = element_rect(fill=NA,colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey87"),
        panel.grid.minor = element_line(colour = "grey87"),
        axis.title=element_text(size=8),
        axis.text.y = element_text(colour = "black",size=8),
        axis.text.x = element_text(angle = 90, hjust = 1, color="black",size=3))+
  labs(x="Transposon", y="Total Transposition Events")
m
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="Family_Frequency.tiff",
       dpi=300,
       width=7.5,
       height=5,
       units="in")