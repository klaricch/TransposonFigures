#!/usr/bin/R


library(dplyr)
library(ggplot2)
library(grid)
library(stringr)

setwd("/Users/kristen/Documents/transposon_figure_data/data")
summarydata <- read.table("kin_matrix_full.txt",header=TRUE)
#count numebr of "1" scores
no_cols=ncol(summarydata)
summarydata<-mutate(summarydata, freq=rowSums(summarydata[2:no_cols],na.rm = TRUE))
#count number of non-NA value
summarydata<-mutate(summarydata,total=rowSums(!is.na(summarydata[2:no_cols])))
#calcualte allele frequency
summarydata<-mutate(summarydata,AF_pre=round((freq/total),digits=3))
summarydata<-mutate(summarydata,call=ifelse(grepl('non-reference',trait),"Novel Site", "Reference Site"))


summarydata<-mutate(summarydata,AF=ifelse(AF_pre>=.50, abs(AF_pre-1), AF_pre))
summarydata$AF_pre
tail(summarydata$AF)
summarydata$AF

m <- ggplot(summarydata, aes(x=AF))
m <- m + geom_bar(binwidth=.05) +
  facet_grid(call ~ .,scale="free_y")+
  theme(panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill="white"),
        axis.ticks = element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=9),
        axis.text.x = element_text(colour = "black",size=9),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9,face="bold"))+
  
  guides(fill=FALSE) +
  labs(x="Allele Frequency", y="Count")+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))
m

setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="Allele Frequency.tiff",dpi=300, width=7.5,height=3.5,units="in")

