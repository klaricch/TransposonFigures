#!/usr/bin/R
# this script plots:
# 1) the FD, FN, and TP count vs Depth of Coverage 
# 2) the FD and TP count vs Read Support while displaying the Read Support Threshold
# 3) Read Support vs Coverage 
# NOTE: One of each of the above 3 graphs are produced for the full, half, quarter, and tenth subsets of the original coverage amount plus one combined plot
# NOTE: This script is for reference calls 
# USE: depth_TPFD_ref.R

library(ggplot2)
library(grid)
library(cowplot)
library(dplyr)

#Coverage references
file_list=c("/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/RSV_SIM",
            "/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/RSV_SIM_half",
            "/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/RSV_SIM_quarter",
            "/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/RSV_SIM_tenth")

for (i in file_list){
  print(i)
  references<-paste(i,"ref",sep="/")
  setwd(references)
  
  fileID<-basename(i)

  
  if (fileID=="RSV_SIM"){ID<-"130x Depth of Coverage"}
  if (fileID=="RSV_SIM_half"){ID<-"65x Depth of Coverage"}
  if (fileID=="RSV_SIM_quarter"){ID<-"32.5x Depth of Coverage"}
  if (fileID=="RSV_SIM_tenth"){ID<-"13x Depth of Coverage"}

  
  
  summarydata <- read.table("summary_depth_TPFD_ref.txt",header=TRUE,row.names=NULL)
  coverage<- read.table("mean_coverage_and_sd.txt",header=TRUE,row.names=NULL)
  
  summarydata$call <- factor(summarydata$call, levels = summarydata$call[order(summarydata$call, decreasing = TRUE)])
  
  a <- ggplot(summarydata, aes(x=coverage,fill=call))
  a <- a + geom_bar(binwidth=5)+ 
    facet_grid(call ~ ., scale="free_y")+
    labs(x="Coverage (Number of Reads)", y="Count", title=ID)+
    theme(strip.background = element_blank(),
          strip.text.x = element_text(size = 8, colour = "black",face="bold"),
          strip.text.y = element_text(size = 8, colour = "black",face="bold",angle=90),
          panel.margin = unit(.6, "lines"),
          panel.border = element_rect(fill=NA,colour = "black"),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(colour = "grey97"),
          axis.title=element_text(size=8),
          axis.text.y = element_text(colour = "black",size=8),
          axis.text.x = element_text(colour = "black",size=8),
          axis.title.y=element_text(vjust=1.75),
          axis.line.y = element_line(colour = "black"),
          axis.line.x = element_line(colour = "black"),
          plot.title = element_text(colour = "black",size=8),
          # legend.text=element_text(size=9),
          legend.position=('none'))+
    geom_point(aes(x=0, y=0),alpha=0)+
    scale_x_continuous(expand = c(0,0),limits=c(0, max(summarydata$coverage)))+
    scale_fill_manual(values = c("purple3", "turquoise3", "darkorange"))
  a
  a <- ggplotGrob(a)
  a$layout[a$layout$name == "strip-right",c("l", "r")] <- 2
  plot(a)
  grid.draw(a)
  
  filename <- paste("depth_TPFD_ref_cov",fileID, sep ="_")
  filename <- paste(filename,"tiff", sep =".")
  ggsave(filename,dpi=300, width=7.5,height=3.5,units="in")
  
  if (fileID=="RSV_SIM"){first<-a}
  if (fileID=="RSV_SIM_half"){second<-a}
  if (fileID=="RSV_SIM_quarter"){third<-a}
  if (fileID=="RSV_SIM_tenth"){fourth<-a}
  
}

setwd("/Users/kristen/Documents/transposon_figure_data/figures")
a_all<-plot_grid(first,second,third,fourth)
ggsave(a_all,filename="Combined_Ref_Cov.tiff",dpi=300, width=7.5,height=5,units="in")
ggsave(a_all,filename="Combined_Ref_Cov.png",dpi=300, width=7.5,height=5,units="in")


########################################################################################################################
########################################################################################################################
########################################################################################################################
# references READ SUPPORT

file_list=c("/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/RSV_SIM",
            "/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/RSV_SIM_half",
            "/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/RSV_SIM_quarter",
            "/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/RSV_SIM_tenth")

for (i in file_list){
  print(i)
  references<-paste(i,"ref",sep="/")
  setwd(references)
  
  fileID<-basename(i)

  
  if (fileID=="RSV_SIM"){ID<-"130x Depth of Coverage"}
  if (fileID=="RSV_SIM_half"){ID<-"65x Depth of Coverage"}
  if (fileID=="RSV_SIM_quarter"){ID<-"32.5x Depth of Coverage"}
  if (fileID=="RSV_SIM_tenth"){ID<-"13x Depth of Coverage"}
  
  
  summarydata <- read.table("summary_depth_TPFD_ref.txt",header=TRUE,row.names=NULL)
  new_df<- summarydata[summarydata$call != "FN",]
  summarydata <- new_df
  summarydata$call <- factor(summarydata$call, levels = summarydata$call[order(summarydata$call, decreasing = TRUE)])
  a <- ggplot(summarydata, aes(x=N,fill=call))
  a <- a + geom_bar(binwidth=5)+
    geom_vline(xintercept=c(3,3), linetype="dashed", color="black")+
    facet_grid(call ~ ., scale="free_y")+
    labs(x="Read Support", y="Count", title=ID)+
    theme(strip.background = element_blank(),
          strip.text.x = element_text(size = 8, colour = "black",face="bold"),
          strip.text.y = element_text(size = 8, colour = "black",face="bold",angle=90),
          panel.margin = unit(.6, "lines"),
          panel.border = element_rect(fill=NA,colour = "black"),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(colour = "grey97"),
          axis.title=element_text(size=8),
          axis.text.y = element_text(colour = "black",size=8),
          axis.text.x = element_text(colour = "black",size=8),
          axis.title.y=element_text(vjust=1.75),
          axis.line.y = element_line(colour = "black"),
          axis.line.x = element_line(colour = "black"),
          plot.title = element_text(colour = "black",size=8),
          # legend.text=element_text(size=9),
          legend.position=('none'))+
    geom_point(aes(x=0, y=0),alpha=0)+
    scale_x_continuous(expand = c(0,0),limits=c(0, max(summarydata$N)))+
    scale_fill_manual(values = c("purple3", "darkorange"))
  a
  a <- ggplotGrob(a)
  a$layout[a$layout$name == "strip-right",c("l", "r")] <- 2
  plot(a)
  grid.draw(a)
  
  filename <- paste("depth_TPFD_ref_RS",fileID, sep ="_")
  filename <- paste(filename,"tiff", sep =".")
  ggsave(filename,dpi=300, width=7.5,height=3.5,units="in")
  
  if (fileID=="RSV_SIM"){first<-a}
  if (fileID=="RSV_SIM_half"){second<-a}
  if (fileID=="RSV_SIM_quarter"){third<-a}
  if (fileID=="RSV_SIM_tenth"){fourth<-a}
  
}

setwd("/Users/kristen/Documents/transposon_figure_data/figures")
a_all<-plot_grid(first,second,third,fourth)
ggsave(a_all,filename="Combined_Ref_R.tiff",dpi=300, width=7.5,height=5,units="in")
ggsave(a_all,filename="Combined_Ref_R.png",dpi=300, width=7.5,height=5,units="in")

########################################################################################################################
########################################################################################################################
########################################################################################################################
#SCATTER ref
file_list=c("/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/RSV_SIM",
            "/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/RSV_SIM_half",
            "/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/RSV_SIM_quarter",
            "/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/RSV_SIM_tenth")

for (i in file_list){
  print(i)
  references<-paste(i,"ref",sep="/")
  setwd(references)
  
  fileID<-basename(i)

  
  if (fileID=="RSV_SIM"){ID<-"130x Depth of Coverage"}
  if (fileID=="RSV_SIM_half"){ID<-"65x Depth of Coverage"}
  if (fileID=="RSV_SIM_quarter"){ID<-"32.5x Depth of Coverage"}
  if (fileID=="RSV_SIM_tenth"){ID<-"13x Depth of Coverage"}
  
  summarydata <- read.table("summary_depth_TPFD_ref.txt",header=TRUE,row.names=NULL)
  coverage<- read.table("mean_coverage_and_sd.txt",header=TRUE,row.names=NULL)
  
  new_df<- summarydata[summarydata$call != "FN",]
  summarydata <- new_df
  
  
  ###CHECK THIS OVER!!!!!!!!
  #calculate percent of False Discoveries removed
  str(summarydata)
  above_FD<-filter(summarydata,N>=3,call=="FD")
  aFD<-nrow(above_FD)
  below_FD<-filter(summarydata,N<3,call=="FD")
  bFD<-nrow(below_FD)
  per_FD=round(bFD/(bFD+aFD)*100,digits=2)
  print(per_FD)
  
  #calculate percent of True Positives removed
  above_TP<-filter(summarydata,N>=3,call=="TP")
  aTP<-nrow(above_TP)
  below_TP<-filter(summarydata,N<3,call=="TP")
  bTP<-nrow(below_TP)
  per_TP=round(bTP/(bTP+aTP)*100,digits=2)
  print(per_TP)
  
  #plot discarded percentages of TP and FD calls the same relative distance above the cutoff line
  #.97 max coverage for x coordinate
  x_ann<-.97*max(summarydata$coverage)
  # for the TP y coordinate:
  Ty_ann<-max(summarydata[summarydata$call == "TP",]$N)
  Ty_ann<-.50*(Ty_ann/3)+3
  # for the FD y coordinate:
  Fy_ann<-max(summarydata[summarydata$call == "FD",]$N)
  Fy_ann<-.50*(Fy_ann/3)+3
  
  #reorder factor levels so that TP is above FD in the facets
  summarydata$call <- factor(summarydata$call, levels = summarydata$call[order(summarydata$call, decreasing = TRUE)])
  a <- ggplot(data = summarydata, aes(x = coverage,y=N, colour=call))
  
  a <- a + geom_point(size=.75)+
    geom_hline(aes(yintercept = 3),linetype="dashed", color="black")+
    facet_grid(call ~ ., scale="free_y")+
    geom_text(data = subset(summarydata, call=="TP"),label = paste(per_TP,"%"),x=x_ann,y=Ty_ann,size=2.2)+
    geom_text(data = subset(summarydata, call=="FD"),label = paste(per_FD,"%"),x=x_ann,y=Fy_ann,size=2.2)+
    theme(strip.background = element_blank(),
          strip.text.x = element_text(size = 8, colour = "black",face="bold"),
          strip.text.y = element_text(margin = margin(t=0,r=0,l=10,b=0),size = 8, colour = "black",face="bold",angle=90),
          axis.title.y = element_text(margin = margin(r=18),face="bold"),
          axis.title.x=element_text(face="bold"),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(colour = "grey97"),
          axis.ticks =element_line(colour = "black"),
          axis.text.y = element_text(colour = "black",size=8),
          axis.text.x = element_text(colour = "black",size=8),
          axis.line=element_line(linetype="solid"),
          axis.line.y = element_line(colour = "black"),
          axis.line.x = element_line(colour = "black"),
          axis.title=element_text(size=8),
          plot.title = element_text(colour = "black",size=8),
          legend.position=('none'))+
    geom_point(aes(x=0, y=0),alpha=0)+
    labs(x="Coverage (Number of Reads)", y= "Read Support", title=ID)+
    scale_colour_manual(values = c("purple3","darkorange"))
  a
  
  if (fileID=="RSV_SIM_quarter"){a<-a+geom_point(aes(x=0,y=4),size=.75,alpha=0)}
  if (fileID=="RSV_SIM_tenth"){a<-a+geom_point(aes(x=0,y=4),size=.75,alpha=0)}
  
  a <- ggplotGrob(a)
  a$layout[a$layout$name == "strip-right",c("l", "r")] <- 2
  plot(a)
  grid.draw(a)
  
  filename <- paste("depth_TPFD_ref_scatter",fileID, sep ="_")
  filename <- paste(filename,"tiff", sep =".")
  ggsave(filename,dpi=300, width=7.5,height=3.5,units="in")
  ggsave(filename,dpi=300, width=7.5,height=3.5,units="in")
  

  if (fileID=="RSV_SIM"){first<-a}
  if (fileID=="RSV_SIM_half"){second<-a}
  if (fileID=="RSV_SIM_quarter"){third<-a}
  if (fileID=="RSV_SIM_tenth"){fourth<-a}
  
}

setwd("/Users/kristen/Documents/transposon_figure_data/figures")
a_all<-plot_grid(first,second,third,fourth,labels=c("A","B","C","D"))
ggsave(a_all,filename="Combined_Ref_Scatter.tiff",dpi=300, width=7.5,height=5,units="in")
ggsave(a_all,filename="Combined_Ref_Scatter.png",dpi=300, width=7.5,height=5,units="in")
