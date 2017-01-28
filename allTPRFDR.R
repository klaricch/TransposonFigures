#!/usr/bin/R
## this script plots the TPR and FDR for each transposon detection method
## produces graphs for both overall transposon detection and family-aware detection
## plots each method seaparately and on the same graph
## creates histograms displaying the distance between the position of the simulated TE and the position of the detected transposon
## NOTE: set x limits separately based on histograms

## USE: (navigate to directory with BEDCOMAPRE results files) graph_TFPN_distances.R

library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)

setwd("/Users/kristen/Documents/transposon_figure_data/simulations/round19")
summarydata <- read.table("BEDCOMPARE_MEANS.txt",header=TRUE)


FAM <- summarydata[ summarydata$fam=="family_aware", ]

#FAMILY_AWARE_TPR
N2_FAM <- FAM[(grep("N2", FAM$M2)), ]
N2_FAM<-mutate(N2_FAM, method=stringr::str_split_fixed(N2_FAM$M2, "_",8)[,6])
N2_FAM<-mutate(N2_FAM, calc=stringr::str_split_fixed(N2_FAM$M2, "_",8)[,8])
N2_FAM<-gather(N2_FAM,rate_name, rate_value, TPR:FDR)

TPR<- N2_FAM[(N2_FAM$rate_name=="TPR") & (N2_FAM$calc=="F"),]
FDR<- N2_FAM[(N2_FAM$rate_name=="FDR"& (N2_FAM$calc=="F")),]
TPRE<- N2_FAM[(N2_FAM$rate_name=="TPR"& (N2_FAM$calc=="F_error")),]
FDRE<- N2_FAM[(N2_FAM$rate_name=="FDR"& (N2_FAM$calc=="F_error")),]

TPR<-merge(TPR,TPRE, by=c("method","DistanceCutoff"))
TPR<-dplyr::select(TPR, method, DistanceCutoff,rate_name.x,rate_value.x,rate_value.y)
colnames(TPR)<-c("Method","DistanceCutoff","Stat","Value","Error")
FDR<-merge(FDR,FDRE, by=c("method","DistanceCutoff"))
FDR<-dplyr::select(FDR, method, DistanceCutoff,rate_name.x,rate_value.x,rate_value.y)
colnames(FDR)<-c("Method","DistanceCutoff","Stat","Value","Error")

N2_FAM<-rbind(TPR,FDR)
N2_FAM$Method <- factor(N2_FAM$Method,
                            levels = c("retroseq", "telocate","temp"),
                            labels = c("RetroSeq", "TE-Locate", "TEMP"))

N2_FAM$Stat <- factor(N2_FAM$Stat, levels = N2_FAM$Stat[order(N2_FAM$Stat, decreasing = TRUE)])
#pdf(file = "FAMILY_AWARE_ALL_TPR.pdf")
a <- ggplot(data = N2_FAM, aes(x = DistanceCutoff, y = Value,colour=Method))+
  geom_errorbar(aes(ymin=Value-Error, ymax=Value+Error),color="gray60") +
 geom_line()+ xlim(0,400)+
  facet_grid(Stat ~ .,scale="free_y")+
theme(strip.background = element_blank(),
      strip.text.x = element_text(size = 9, colour = "black",face="bold"),
      strip.text.y = element_text(size = 9, colour = "black",face="bold",angle=90),
      panel.background = element_rect(fill = "white"),
      axis.ticks =element_line(colour = "black"),
      axis.text.y = element_text(colour = "black",size=9),
      axis.text.x = element_text(colour = "black",size=9),
      axis.line=element_line(linetype="solid"),
      axis.line.y = element_line(colour = "black"),
      axis.line.x = element_line(colour = "black"),
      axis.title=element_text(size=9),
      axis.title.x=element_text(face="bold"),
      axis.title.y=element_text(vjust=1.75,face="bold", colour="white"),
      legend.title=element_blank(),
      legend.background = element_rect(fill=FALSE),
      legend.key=element_rect(fill=NA),
      plot.title = element_text(size=9))+
  scale_colour_manual(values = c("tan1","paleturquoise1","indianred1"))
a

a <- ggplotGrob(a)
a$layout[a$layout$name == "strip-right",c("l", "r")] <- 2
plot(a)
grid.draw(a)

setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(a, filename="te_progs.tiff",dpi=300, width=7.5,height=3.5,units="in")
ggsave(a, filename="te_progs.png",dpi=300, width=7.5,height=3.5,units="in")
