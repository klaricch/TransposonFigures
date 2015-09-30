#!/usr/bin/R
# this script plots:
# 1) the FD, FN, and TP count vs Depth of Coverage
# 2) the FD and TP count vs Read Support
# 3) Read Support vs Coverage
# NOTE: One of each of the above 3 graphs are produced for the full, half, quarter, and tenth subsets of the original coverage amount plus one combined plot
# NOTE: This script is for reference calls
# USE: depth_TPFD_ref.R

library(ggplot2)
library(grid)
library(cowplot)
library(dplyr)

# references full
setwd("/Users/kristen/Desktop/Fig_p/depth_TPFD_files/RSV_SIM/ref")
summarydata <- read.table("summary_depth_TPFD_ref.txt",header=TRUE,row.names=NULL)
coverage<- read.table("mean_coverage_and_sd.txt",header=TRUE,row.names=NULL)

m1 <- ggplot(summarydata, aes(x=coverage,fill=call))
m1 <- m1 + geom_bar(binwidth=5)+ 
  facet_grid(call ~ ., scale="free_y")+
  labs(x="Coverage (Number of Reads)", y="Count", title="130x Depth of Coverage")+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 8, colour = "black",face="bold"),
        panel.margin = unit(.6, "lines"),
        panel.border = element_rect(fill=NA,colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey97"),
        axis.title=element_text(size=8),
        axis.text.y = element_text(colour = "black",size=8),
        axis.text.x = element_text(colour = "black",size=8),
        plot.title = element_text(colour = "black",size=8),
        # legend.text=element_text(size=9),
        legend.position=('none'))+
  geom_point(aes(x=0, y=0),alpha=0)+
  scale_x_continuous(expand = c(0,0),limits=c(0, max(summarydata$coverage)))+
  scale_fill_manual(values = c("darkorange", "turquoise3", "slateblue1"))
m1
ggsave(filename="depth_TPFD_ref_full.tiff",dpi=300, width=7.5,height=3.5,units="in")

# references half
setwd("/Users/kristen/Desktop/Fig_p/depth_TPFD_files/RSV_SIM_half/ref")
summarydata <- read.table("summary_depth_TPFD_ref.txt",header=TRUE,row.names=NULL)
coverage<- read.table("mean_coverage_and_sd.txt",header=TRUE,row.names=NULL)

m2 <- ggplot(summarydata, aes(x=coverage,fill=call))
m2 <- m2 + geom_bar(binwidth=5)+ 
  facet_grid(call ~ ., scale="free_y")+
  labs(x="Coverage (Number of Reads)", y="Count", title="65x Depth of Coverage")+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 8, colour = "black",face="bold"),
        panel.margin = unit(.6, "lines"),
        panel.border = element_rect(fill=NA,colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey97"),
        axis.title=element_text(size=8),
        axis.text.y = element_text(colour = "black",size=8),
        axis.text.x = element_text(colour = "black",size=8),
        plot.title = element_text(colour = "black",size=8),
        # legend.text=element_text(size=9),
        legend.position=('none'))+
  geom_point(aes(x=0, y=0),alpha=0)+
  scale_x_continuous(expand = c(0,0),limits=c(0, max(summarydata$coverage)))+
  scale_fill_manual(values = c("darkorange", "turquoise3", "slateblue1"))
m2
ggsave(filename="depth_TPFD_ref_half.tiff",dpi=300, width=7.5,height=3.5,units="in")

# references quarter
setwd("/Users/kristen/Desktop/Fig_p/depth_TPFD_files/RSV_SIM_quarter/ref")
summarydata <- read.table("summary_depth_TPFD_ref.txt",header=TRUE,row.names=NULL)
coverage<- read.table("mean_coverage_and_sd.txt",header=TRUE,row.names=NULL)

m3 <- ggplot(summarydata, aes(x=coverage,fill=call))
m3 <- m3 + geom_bar(binwidth=5)+ 
  facet_grid(call ~ ., scale="free_y")+
  labs(x="Coverage (Number of Reads)", y="Count", title="32.5x Depth of Coverage")+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 8, colour = "black",face="bold"),
        panel.margin = unit(.6, "lines"),
        panel.border = element_rect(fill=NA,colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey97"),
        axis.title=element_text(size=8),
        axis.text.y = element_text(colour = "black",size=8),
        axis.text.x = element_text(colour = "black",size=8),
        plot.title = element_text(colour = "black",size=8),
        # legend.text=element_text(size=9),
        legend.position=('none'))+
  geom_point(aes(x=0, y=0),alpha=0)+
  scale_x_continuous(expand = c(0,0),limits=c(0, max(summarydata$coverage)))+
  scale_fill_manual(values = c("darkorange", "turquoise3", "slateblue1"))
m3
ggsave(filename="depth_TPFD_ref_quarter.tiff",dpi=300, width=7.5,height=3.5,units="in")

# references tenth
setwd("/Users/kristen/Desktop/Fig_p/depth_TPFD_files/RSV_SIM_tenth/ref")
summarydata <- read.table("summary_depth_TPFD_ref.txt",header=TRUE,row.names=NULL)
coverage<- read.table("mean_coverage_and_sd.txt",header=TRUE,row.names=NULL)

m4 <- ggplot(summarydata, aes(x=coverage,fill=call))
m4 <- m4 + geom_bar(binwidth=5)+ 
  facet_grid(call ~ ., scale="free_y")+
  labs(x="Coverage (Number of Reads)", y="Count", title="13x Depth of Coverage")+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 8, colour = "black",face="bold"),
        panel.margin = unit(.6, "lines"),
        panel.border = element_rect(fill=NA,colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey97"),
        axis.title=element_text(size=8),
        axis.text.y = element_text(colour = "black",size=8),
        axis.text.x = element_text(colour = "black",size=8),
        plot.title = element_text(colour = "black",size=8),
        # legend.text=element_text(size=9),
        legend.position=('none'))+
  geom_point(aes(x=0, y=0),alpha=0)+
  scale_x_continuous(expand = c(0,0),limits=c(0, max(summarydata$coverage)))+
  scale_fill_manual(values = c("turquoise3", "slateblue1"))

m4
ggsave(filename="depth_TPFD_ref_tenth.tiff",dpi=300, width=7.5,height=3.5,units="in")

m_all<-plot_grid(m1,m2,m3,m4)
setwd("/Users/kristen/Desktop/Fig_p")
ggsave(filename="combined_ref_F.tiff",dpi=300, width=7.5,height=5,units="in")

########################################################################################################################
########################################################################################################################
########################################################################################################################

# references READ SUPPORT full
setwd("/Users/kristen/Desktop/Fig_p/depth_TPFD_files/RSV_SIM/ref")
summarydata <- read.table("summary_depth_TPFD_ref.txt",header=TRUE,row.names=NULL)
new_df<- summarydata[summarydata$call != "FN",]
summarydata <- new_df
m <- ggplot(summarydata, aes(x=N,fill=call))
m <- m + geom_bar(binwidth=50)+
 # geom_vline(xintercept=c(8,8), linetype="dashed", color="black")+
  facet_grid(call ~ ., scale="free_y")+
  labs(x="Read Support", y="Count")+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        panel.margin = unit(.6, "lines"),
        panel.border = element_rect(fill=NA,colour = "black"),
        panel.background = element_rect(fill = "white"),
        axis.title=element_text(size=9),
        axis.text.y = element_text(colour = "black",size=9),
        axis.text.x = element_text(colour = "black",size=9),
        legend.text=element_text(size=9))
m
ggsave(filename="depth_TPFD_ref_RS_full.tiff",dpi=300, width=7.5,height=3.5,units="in")

# references READ SUPPORT full Diff Bin
setwd("/Users/kristen/Desktop/Fig_p/depth_TPFD_files/RSV_SIM/ref")
summarydata <- read.table("summary_depth_TPFD_ref.txt",header=TRUE,row.names=NULL)
new_df<- summarydata[summarydata$call != "FN",]
summarydata <- new_df
m1 <- ggplot(summarydata, aes(x=N,fill=call))
m1 <- m1 + geom_bar(binwidth=5)+
  geom_vline(xintercept=c(3,3), linetype="dashed", color="black")+
  facet_grid(call ~ ., scale="free_y")+
  labs(x="Read Support", y="Count", title="130x Depth of Coverage")+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 8, colour = "black",face="bold"),
        panel.margin = unit(.6, "lines"),
        panel.border = element_rect(fill=NA,colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey97"),
        axis.title=element_text(size=8),
        axis.text.y = element_text(colour = "black",size=8),
        axis.text.x = element_text(colour = "black",size=8),
        plot.title = element_text(colour = "black",size=8),
        # legend.text=element_text(size=9),
        legend.position=('none'))+
  geom_point(aes(x=0, y=0),alpha=0)+
  scale_x_continuous(expand = c(0,0),limits=c(0, max(summarydata$N)))+
  scale_fill_manual(values = c("darkorange", "slateblue1"))
m1
ggsave(filename="depth_TPFD_ref_RS_full_BIN.tiff",dpi=300, width=7.5,height=3.5,units="in")

# references READ SUPPORT half
setwd("/Users/kristen/Desktop/Fig_p/depth_TPFD_files/RSV_SIM_half/ref")
summarydata <- read.table("summary_depth_TPFD_ref.txt",header=TRUE,row.names=NULL)
new_df<- summarydata[summarydata$call != "FN",]
summarydata <- new_df
m2 <- ggplot(summarydata, aes(x=N,fill=call))
m2 <- m2 + geom_bar(binwidth=5)+
  geom_vline(xintercept=c(3,3), linetype="dashed", color="black")+
  facet_grid(call ~ ., scale="free_y")+
  labs(x="Read Support", y="Count", title="65x Depth of Coverage")+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 8, colour = "black",face="bold"),
        panel.margin = unit(.6, "lines"),
        panel.border = element_rect(fill=NA,colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey97"),
        axis.title=element_text(size=8),
        axis.text.y = element_text(colour = "black",size=8),
        axis.text.x = element_text(colour = "black",size=8),
        plot.title = element_text(colour = "black",size=8),
        # legend.text=element_text(size=9),
        legend.position=('none'))+
  geom_point(aes(x=0, y=0),alpha=0)+
  scale_x_continuous(expand = c(0,0),limits=c(0, max(summarydata$N)))+
  scale_fill_manual(values = c("darkorange", "slateblue1"))
m2
ggsave(filename="depth_TPFD_ref_RS_half.tiff",dpi=300, width=7.5,height=3.5,units="in")

# references READ SUPPORT quarter
setwd("/Users/kristen/Desktop/Fig_p/depth_TPFD_files/RSV_SIM_quarter/ref")
summarydata <- read.table("summary_depth_TPFD_ref.txt",header=TRUE,row.names=NULL)
new_df<- summarydata[summarydata$call != "FN",]
summarydata <- new_df
m3 <- ggplot(summarydata, aes(x=N,fill=call))
m3 <- m3 + geom_bar(binwidth=5)+
  geom_vline(xintercept=c(3,3), linetype="dashed", color="black")+
  facet_grid(call ~ ., scale="free_y")+
  labs(x="Read Support", y="Count", title="32.5x Depth of Coverage")+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 8, colour = "black",face="bold"),
        panel.margin = unit(.6, "lines"),
        panel.border = element_rect(fill=NA,colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey97"),
        axis.title=element_text(size=8),
        axis.text.y = element_text(colour = "black",size=8),
        axis.text.x = element_text(colour = "black",size=8),
        plot.title = element_text(colour = "black",size=8),
        # legend.text=element_text(size=9),
        legend.position=('none'))+
  geom_point(aes(x=0, y=0),alpha=0)+
  scale_x_continuous(expand = c(0,0),limits=c(0, max(summarydata$N)))+
  scale_fill_manual(values = c("darkorange", "slateblue1"))
m3
ggsave(filename="depth_TPFD_ref_RS_quarter.tiff",dpi=300, width=7.5,height=3.5,units="in")

# references READ SUPPORT tenth
setwd("/Users/kristen/Desktop/Fig_p/depth_TPFD_files/RSV_SIM_tenth/ref")
summarydata <- read.table("summary_depth_TPFD_ref.txt",header=TRUE,row.names=NULL)
new_df<- summarydata[summarydata$call != "FN",]
summarydata <- new_df
m4 <- ggplot(summarydata, aes(x=N,fill=call))
m4 <- m4 + geom_bar(binwidth=5)+
  geom_vline(xintercept=c(3,3), linetype="dashed", color="black")+
  facet_grid(call ~ ., scale="free_y")+
  labs(x="Read Support", y="Count", title="13x Depth of Coverage")+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 8, colour = "black",face="bold"),
        panel.margin = unit(.6, "lines"),
        panel.border = element_rect(fill=NA,colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey97"),
        axis.title=element_text(size=8),
        axis.text.y = element_text(colour = "black",size=8),
        axis.text.x = element_text(colour = "black",size=8),
        plot.title = element_text(colour = "black",size=8),
        # legend.text=element_text(size=9),
        legend.position=('none'))+
  geom_point(aes(x=0, y=0),alpha=0)+
  scale_x_continuous(expand = c(0,0),limits=c(0, max(summarydata$N)))+
  scale_fill_manual(values = c("slateblue1"))

m4
ggsave(filename="depth_TPFD_ref_RS_tenth.tiff",dpi=300, width=7.5,height=3.5,units="in")
m_all<-plot_grid(m1,m2,m3,m4)
setwd("/Users/kristen/Desktop/Fig_p")
ggsave(filename="combined_ref_RS_F.tiff",dpi=300, width=7.5,height=5,units="in")

########################################################################################################################
########################################################################################################################
########################################################################################################################

#SCATTER ref full
setwd("/Users/kristen/Desktop/Fig_p/depth_TPFD_files/RSV_SIM/ref")
summarydata <- read.table("summary_depth_TPFD_ref.txt",header=TRUE,row.names=NULL)
coverage<- read.table("mean_coverage_and_sd.txt",header=TRUE,row.names=NULL)


new_df<- summarydata[summarydata$call != "FN",]
summarydata <- new_df

#calculate percent of False Discoveries removed
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

#.97 max coverage
x_ann<-.97*max(summarydata$coverage)
# .90 max N of TP calls
Ty_ann<-.90*max(summarydata[summarydata$call == "TP",]$N)
# .90 max N of FD calls
Fy_ann<-.90*max(summarydata[summarydata$call == "FD",]$N)


a1 <- ggplot(data = summarydata, aes(x = coverage,y=N, colour=call))
a1 <- a1 + geom_point(size=.75)+
  #geom_hline(xintercept=c(8,8), linetype="dashed", color="black")+
  #if onyl adding to one facet:
  #geom_hline(data = summarydata[summarydata$call=="FD",], aes(yintercept = 8))+
  geom_hline(aes(yintercept = 3),linetype="dashed", color="black")+
  facet_grid(call ~ ., scale="free_y")+
  geom_text(data = subset(summarydata, call=="TP"),label = paste(per_TP,"%"),x=x_ann,y=Ty_ann,size=2)+
  geom_text(data = subset(summarydata, call=="FD"),label = paste(per_FD,"%"),x=x_ann,y=Fy_ann,size=2)+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 8, colour = "black",face="bold"),
        strip.text.y = element_text(size = 8, colour = "black",face="bold"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey97"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=8),
        axis.text.x = element_text(colour = "black",size=8),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=8),
        plot.title = element_text(colour = "black",size=8),
        #legend.position=c(.90,0.85),
        #legend.background = element_rect(fill=FALSE),
        #legend.key=element_rect(fill=NA),
        # legend.text=element_text(size=9),
        legend.position=('none'))+
  geom_point(aes(x=0, y=0),alpha=0)+
  labs(x="Coverage (Number of Reads)", y= "Read Support", title="130x Depth of Coverage")+
  #scale_x_continuous(expand = c(0,0),limits=c(0, max(summarydata$coverage)))+
  scale_colour_manual(values = c("darkorange", "slateblue1"))
a1
ggsave(filename="depth_TPFD_ref_scatter_full.tiff",dpi=300, width=7.5,height=3.5,units="in")

#SCATTER ref half
setwd("/Users/kristen/Desktop/Fig_p/depth_TPFD_files/RSV_SIM_half/ref")
summarydata <- read.table("summary_depth_TPFD_ref.txt",header=TRUE,row.names=NULL)
coverage<- read.table("mean_coverage_and_sd.txt",header=TRUE,row.names=NULL)

###either merge a new column or add in separate data...
new_df<- summarydata[summarydata$call != "FN",]
summarydata <- new_df

#calculate percent of False Discoveries removed
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

#.97 max coverage
x_ann<-.97*max(summarydata$coverage)
# .90 max N of TP calls
Ty_ann<-.90*max(summarydata[summarydata$call == "TP",]$N)
# .90 max N of FD calls
Fy_ann<-.90*max(summarydata[summarydata$call == "FD",]$N)


a2 <- ggplot(data = summarydata, aes(x = coverage,y=N, colour=call))
a2 <- a2 + geom_point(size=.75)+
  geom_hline(aes(yintercept = 3),linetype="dashed", color="black")+
  facet_grid(call ~ ., scale="free_y")+
  geom_text(data = subset(summarydata, call=="TP"),label = paste(per_TP,"%"),x=x_ann,y=Ty_ann,size=2)+
  geom_text(data = subset(summarydata, call=="FD"),label = paste(per_FD,"%"),x=x_ann,y=Fy_ann,size=2)+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 8, colour = "black",face="bold"),
        strip.text.y = element_text(size = 8, colour = "black",face="bold"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey97"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=8),
        axis.text.x = element_text(colour = "black",size=8),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=8),
        plot.title = element_text(colour = "black",size=8),
        #legend.position=c(.90,0.85),
        #legend.background = element_rect(fill=FALSE),
        #legend.key=element_rect(fill=NA),
        # legend.text=element_text(size=9),
        legend.position=('none'))+
  geom_point(aes(x=0, y=0),alpha=0)+
  labs(x="Coverage (Number of Reads)", y= "Read Support", title="65x Depth of Coverage")+
  #scale_x_continuous(expand = c(0,0),limits=c(0, max(summarydata$coverage)))+
  scale_colour_manual(values = c("darkorange", "slateblue1"))
a2
ggsave(filename="depth_TPFD_ref_scatter_half.tiff",dpi=300, width=7.5,height=3.5,units="in")

#SCATTER ref quarter
setwd("/Users/kristen/Desktop/Fig_p/depth_TPFD_files/RSV_SIM_quarter/ref")
summarydata <- read.table("summary_depth_TPFD_ref.txt",header=TRUE,row.names=NULL)
coverage<- read.table("mean_coverage_and_sd.txt",header=TRUE,row.names=NULL)

new_df<- summarydata[summarydata$call != "FN",]
summarydata <- new_df

#calculate percent of False Discoveries removed
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

#.97 max coverage
x_ann<-.97*max(summarydata$coverage)
# .90 max N of TP calls
Ty_ann<-.90*max(summarydata[summarydata$call == "TP",]$N)
# .90 max N of FD calls
Fy_ann<-.90*max(summarydata[summarydata$call == "FD",]$N)


a3 <- ggplot(data = summarydata, aes(x = coverage,y=N, colour=call))
a3 <- a3 + geom_point(size=.75)+
  geom_hline(aes(yintercept = 3),linetype="dashed", color="black")+
  facet_grid(call ~ ., scale="free_y")+
  geom_text(data = subset(summarydata, call=="TP"),label = paste(per_TP,"%"),x=x_ann,y=Ty_ann,size=2)+
  geom_text(data = subset(summarydata, call=="FD"),label = paste(per_FD,"%"),x=x_ann,y=Fy_ann,size=2)+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 8, colour = "black",face="bold"),
        strip.text.y = element_text(size = 8, colour = "black",face="bold"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey97"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=8),
        axis.text.x = element_text(colour = "black",size=8),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=8),
        plot.title = element_text(colour = "black",size=8),
        #legend.position=c(.90,0.85),
        #legend.background = element_rect(fill=FALSE),
        #legend.key=element_rect(fill=NA),
        # legend.text=element_text(size=9),
        legend.position=('none'))+
  geom_point(aes(x=0, y=0),alpha=0)+
  labs(x="Coverage (Number of Reads)", y= "Read Support", title="32.5x Depth of Coverage")+
  #scale_x_continuous(expand = c(0,0),limits=c(0, max(summarydata$coverage)))+
  scale_colour_manual(values = c("darkorange", "slateblue1"))
a3
ggsave(filename="depth_TPFD_ref_scatter_quarter.tiff",dpi=300, width=7.5,height=3.5,units="in")

#SCATTER ref tenth
setwd("/Users/kristen/Desktop/Fig_p/depth_TPFD_files/RSV_SIM_tenth/ref")
summarydata <- read.table("summary_depth_TPFD_ref.txt",header=TRUE,row.names=NULL)
coverage<- read.table("mean_coverage_and_sd.txt",header=TRUE,row.names=NULL)

new_df<- summarydata[summarydata$call != "FN",]
summarydata <- new_df

#calculate percent of False Discoveries removed
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

#.97 max coverage
x_ann<-.97*max(summarydata$coverage)
# .90 max N of TP calls
Ty_ann<-.90*max(summarydata[summarydata$call == "TP",]$N)
# .90 max N of FD calls
Fy_ann<-.90*max(summarydata[summarydata$call == "FD",]$N)


a4 <- ggplot(data = summarydata, aes(x = coverage,y=N, colour=call))
a4 <- a4 + geom_point(size=.75)+
  geom_hline(aes(yintercept = 3),linetype="dashed", color="black")+
  facet_grid(call ~ ., scale="free_y")+
  geom_text(data = subset(summarydata, call=="TP"),label = paste(per_TP,"%"),x=x_ann,y=Ty_ann,size=2)+
  geom_text(data = subset(summarydata, call=="FD"),label = paste(per_FD,"%"),x=x_ann,y=Fy_ann,size=2)+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 8, colour = "black",face="bold"),
        strip.text.y = element_text(size = 8, colour = "black",face="bold"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey97"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=8),
        axis.text.x = element_text(colour = "black",size=8),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=8),
        plot.title = element_text(colour = "black",size=8),
        #legend.position=c(.90,0.85),
        #legend.background = element_rect(fill=FALSE),
        #legend.key=element_rect(fill=NA),
        # legend.text=element_text(size=9),
        legend.position=('none'))+
  geom_point(aes(x=0, y=0),alpha=0)+
  labs(x="Coverage (Number of Reads)", y= "Read Support", title="13x Depth of Coverage")+
  #scale_x_continuous(expand = c(0,0),limits=c(0, max(summarydata$coverage)))+
  scale_colour_manual(values = c("slateblue1"))
a4
ggsave(filename="depth_TPFD_ref_scatter_tenth.tiff",dpi=300, width=7.5,height=3.5,units="in")

m_all<-plot_grid(a1,a2,a3,a4)
setwd("/Users/kristen/Desktop/Fig_p")
ggsave(filename="combined_ref_scatter_F.tiff",dpi=300, width=7.5,height=5,units="in")
