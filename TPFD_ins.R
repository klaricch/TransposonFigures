#!/usr/bin/R
# this script plots:
# 1) the FD, FN, and TP count vs Depth of Coverage 
# 2) the FD and TP count vs Read Support while displaying the Read Support Threshold
# 3) Read Support vs Coverage 
# NOTE: One of each of the above 3 graphs are produced for the full, half, quarter, and tenth subsets of the original coverage amount plus one combined plot
# NOTE: This script is for insertion calls 
# USE: depth_TPFD_ins.R

library(ggplot2)
library(grid)
library(cowplot)
library(dplyr)



#Coverage absences
file_list=c("/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/feb/round20",
            "/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/feb/round21",
            "/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/feb/round22",
            "/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/feb/round23")

for (i in file_list){
  print(i)
  setwd(i)
  
  fileID<-basename(i)
  
  if (fileID=="round20"){ID<-"130x Depth of Coverage"}
  if (fileID=="round23"){ID<-"65x Depth of Coverage"}
  if (fileID=="round22"){ID<-"32.5x Depth of Coverage"}
  if (fileID=="round21"){ID<-"13x Depth of Coverage"}
  
  
  
  summarydata <- read.table("summary_depth_TPFD_ins.txt",header=TRUE,row.names=NULL)
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
  
  filename <- paste("depth_TPFD_ins_cov",fileID, sep ="_")
  filename <- paste(filename,"tiff", sep =".")
  ggsave(filename,dpi=300, width=7.5,height=3.5,units="in")
  
  if (fileID=="round20"){first<-a}
  if (fileID=="round23"){second<-a}
  if (fileID=="round22"){third<-a}
  if (fileID=="round21"){fourth<-a}
  
}

setwd("/Users/kristen/Documents/transposon_figure_data/figures")
a_all<-plot_grid(first,second,third,fourth)
ggsave(a_all,filename="Combined_Ins_Cov.tiff",dpi=300, width=7.5,height=5,units="in")
ggsave(a_all,filename="Combined_Ins_Cov.png",dpi=300, width=7.5,height=5,units="in")


########################################################################################################################
########################################################################################################################
########################################################################################################################
# absences READ SUPPORT

file_list=c("/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/feb/round20",
            "/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/feb/round23",
            "/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/feb/round22",
            "/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/feb/round21")

for (i in file_list){
  print(i)
  setwd(i)
  
  fileID<-basename(i)

  
  if (fileID=="round20"){ID<-"130x Depth of Coverage"}
  if (fileID=="round23"){ID<-"65x Depth of Coverage"}
  if (fileID=="round22"){ID<-"32.5x Depth of Coverage"}
  if (fileID=="round21"){ID<-"13x Depth of Coverage"}
  
  
  summarydata <- read.table("summary_depth_TPFD_ins.txt",header=TRUE,row.names=NULL)
  new_df<- summarydata[summarydata$call != "FN",]
  summarydata <- new_df
  summarydata$call <- factor(summarydata$call, levels = summarydata$call[order(summarydata$call, decreasing = TRUE)])
  a <- ggplot(summarydata, aes(x=N,fill=call))
  a <- a + geom_bar(binwidth=5)+
    geom_vline(xintercept=c(8,8), linetype="dashed", color="black")+
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
  
  filename <- paste("depth_TPFD_ins_RS",fileID, sep ="_")
  filename <- paste(filename,"tiff", sep =".")
  ggsave(filename,dpi=300, width=7.5,height=3.5,units="in")
  
  if (fileID=="round20"){first<-a}
  if (fileID=="round23"){second<-a}
  if (fileID=="round22"){third<-a}
  if (fileID=="round21"){fourth<-a}
  
}

setwd("/Users/kristen/Documents/transposon_figure_data/figures")
a_all<-plot_grid(first,second,third,fourth)
ggsave(a_all,filename="Combined_Ins_R.tiff",dpi=300, width=7.5,height=5,units="in")
ggsave(a_all,filename="Combined_Ins_R.png",dpi=300, width=7.5,height=5,units="in")

########################################################################################################################
########################################################################################################################
########################################################################################################################

#SCATTER abs
file_list=c("/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/feb/round20",
            "/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/feb/round23",
            "/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/feb/round22",
            "/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/feb/round21")

for (i in file_list){
  print(i)
  setwd(i)
  
  fileID<-basename(i)
  
  if (fileID=="round20"){ID<-"130x Depth of Coverage"}
  if (fileID=="round23"){ID<-"65x Depth of Coverage"}
  if (fileID=="round22"){ID<-"32.5x Depth of Coverage"}
  if (fileID=="round21"){ID<-"13x Depth of Coverage"}
  
  summarydata <- read.table("summary_depth_TPFD_ins.txt",header=TRUE,row.names=NULL)
  coverage<- read.table("mean_coverage_and_sd.txt",header=TRUE,row.names=NULL)
  
  new_df<- summarydata[summarydata$call != "FN",]
  summarydata <- new_df
  
  
  
  
  # EDIT HERE
  nr_dir<-paste(i,"nr_files",sep="/")
  setwd(nr_dir)
  all_runs<-read.table("all_runs.txt")
  colnames(all_runs)<-c("row.names","CHROM","Start", "End", "TE", "Read_Support","Orient","Variant_Support")
  all_runs<-dplyr::select(all_runs,row.names,TE,Read_Support,Variant_Support)
  summarydata<-left_join(summarydata,all_runs,by=c("row.names","TE"))
  
  
  
  
  ###CHECK THIS OVER!!!!!!!!
  #calculate percent of False Discoveries removed
  str(summarydata)
  above_FD<-filter(summarydata,N>=8,Variant_Support>.25,call=="FD")
  aFD<-nrow(above_FD)
  below_FD<-filter(summarydata,N<8|Variant_Support<=.25,call=="FD")
  bFD<-nrow(below_FD)
  per_FD=round(bFD/(bFD+aFD)*100,digits=2)

  
  #calculate percent of True Positives removed
  above_TP<-filter(summarydata,N>=8,Variant_Support>.25,call=="TP")
  aTP<-nrow(above_TP)
  below_TP<-filter(summarydata,N<8|Variant_Support<=.25,call=="TP")
  bTP<-nrow(below_TP)
  per_TP=round(bTP/(bTP+aTP)*100,digits=2)
  

  #per_TP
  #per_FD
  #plot discarded percentages of TP and FD calls the same relative distance above the cutoff line
  #.97 max coverage for x coordinate
  x_ann<-.97*max(summarydata$coverage)
  # for the TP y coordinate:
  Ty_ann<-max(summarydata[summarydata$call == "TP",]$N)
  Ty_ann<-.60*(Ty_ann/8)+8
  # for the FD y coordinate:
  Fy_ann<-max(summarydata[summarydata$call == "FD",]$N)
  Fy_ann<-.60*(Fy_ann/8)+8
  
  #reorder factor levels so that TP is above FD in the facets
  
  #set colors
  summarydata <- mutate(summarydata,color_code=ifelse(Variant_Support>.25 & call=="TP" , 'purple3', ifelse(Variant_Support>.25 & call=="FD","darkorange",ifelse(Variant_Support<=.25 & call=="TP", "orchid3","goldenrod1"))))
  
  #unique(summarydata$color_code)
  #unique(summarydata$call)
  #test1<-filter(summarydata,Variant_Support>.25 & call=="TP")
  #test2<-filter(summarydata,Variant_Support<=.25 & call=="TP")
  #test3<-filter(summarydata,Variant_Support>.25 & call=="FD")
  #test4<-filter(summarydata,Variant_Support<=.25 & call=="FD")
  
  #summarydata<-arrange(summarydata,desc(call))
  summarydata$call <- factor(summarydata$call, levels = summarydata$call[order(summarydata$call, decreasing = TRUE)])
  a <- ggplot(data = summarydata, aes(x = coverage,y=N, colour=color_code))
  
  a <- a + geom_point(size=.75)+
  #a<-a + geom_point(aes( color=ifelse(Variant_Support>=.25 , 'red', 'black')),size=.75)+
    
    
    
    
    geom_hline(aes(yintercept = 8),linetype="dashed", color="black")+
    facet_grid(call ~ ., scale="free_y")+
    geom_text(data = subset(summarydata, call=="TP"),label = paste(per_TP,"%"),x=x_ann,y=Ty_ann,size=2.2,color="purple3")+
    geom_text(data = subset(summarydata, call=="FD"),label = paste(per_FD,"%"),x=x_ann,y=Fy_ann,size=2.2,color="darkorange")+
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
    labs(x="Coverage (Number of Reads)", y= "Read Support", title=ID,face="bold")+
    scale_colour_manual(values = c('purple3'="purple3",'darkorange'="darkorange",'orchid3'="orchid3",'goldenrod1'="goldenrod1"))
  a

  
  a <- ggplotGrob(a)
  a$layout[a$layout$name == "strip-right",c("l", "r")] <- 2
  plot(a)
  grid.draw(a)
  
  filename <- paste("depth_TPFD_ins_scatter",fileID, sep ="_")
  filename <- paste(filename,"tiff", sep =".")
  ggsave(filename,dpi=300, width=7.5,height=3.5,units="in")
  ggsave(filename,dpi=300, width=7.5,height=3.5,units="in")
  
  if (fileID=="round20"){first<-a}
  if (fileID=="round23"){second<-a}
  if (fileID=="round22"){third<-a}
  if (fileID=="round21"){fourth<-a}
  
}

setwd("/Users/kristen/Documents/transposon_figure_data/figures")
a_all<-plot_grid(first,second,third,fourth,labels=c("A","B","C","D"))

a_all
ggsave(a_all,filename="Combined_Ins_Scatter.tiff",dpi=300, width=7.5,height=5,units="in")
ggsave(a_all,filename="Combined_Ins_Scatter.png",dpi=300, width=7.5,height=5,units="in")
