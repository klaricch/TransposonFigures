#!/usr/bin/R
# this script plots:
# 1) the TPR and FDR vs the read support threshold for the TEMP absence caller
# 2) the TPR and FDR vs the read support threshold for the TELOCATE reference caller
# USE: graph_BEDCOMPARE_MEAN_RSV_SIMs.R

library(tidyr)
library(dplyr)
library(ggplot2)
library(gtable)
library(cowplot)

#Read Support Cutoff for TEMP
file_list=c("/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/RSV_SIM",
            "/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/RSV_SIM_half",
            "/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/RSV_SIM_quarter",
            "/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/RSV_SIM_tenth")

#loop through the coverage level simulations
for (i in file_list){
  print(i)
  
setwd(i)
    
fileID<-basename(i)
print(file)

if (fileID=="RSV_SIM"){ID<-"130x Depth of Coverage"}
if (fileID=="RSV_SIM_half"){ID<-"65x Depth of Coverage"}
if (fileID=="RSV_SIM_quarter"){ID<-"32.5x Depth of Coverage"}
if (fileID=="RSV_SIM_tenth"){ID<-"13x Depth of Coverage"}


summarydata <- read.table("BEDCOMPARE_MEANS.txt",header=TRUE)

TEMP <- summarydata[(summarydata$M2=="new_CT_run_3_N2_temp_nonredundant.bed_F")& summarydata$fam=="family_aware", ]
SD_TEMP <- summarydata[(summarydata$M2=="new_CT_run_3_N2_temp_nonredundant.bed_F_error")& summarydata$fam=="family_aware", ]
merged_TEMP<-merge(TEMP,SD_TEMP, by="DistanceCutoff")

merged_ALL<-merged_TEMP

#collapse 2 columns to prep for faceting
names(merged_ALL)
merged<-merged_ALL %>%
  gather(rate_name, rate_value, TPR.x:FDR.x)
#gather(rate_name, rate_value, TPR.x:FDR.x, FNR.x)
#separate TPR and FDR and associate with the proper error number
names(merged)
TPR<- merged[(merged$rate_name=="TPR.x"),c(1:7,8,10:11)]
FDR<- merged[(merged$rate_name=="FDR.x"),c(1:7,9,10:11)]
#rename TPR.y and FDR.y to error
#rename TPR.y and FDR.y to error
colnames(TPR)[8] <- "error"
colnames(FDR)[8] <- "error"

names(TPR)
names(FDR)
final <-rbind(TPR,FDR)

names(final)
levels(final$rate_name)

final$rate_name <- factor(final$rate_name,
                          levels = c("TPR.x", "FDR.x"),
                          labels = c("TPR", "FDR"))

a <- ggplot(data = final, aes(x = DistanceCutoff, y = rate_value))
a <- a + geom_line(color="black")+ xlim(0,31)+
  geom_errorbar(aes(ymin=rate_value-error, ymax=rate_value+error),color="black") +
  geom_point(data = subset(final, rate_name=="TPR"),aes(x=0, y=100),alpha=0)+
  geom_point(data = subset(final, rate_name=="FDR"),aes(x=0, y=0),alpha=0)+
  geom_vline(xintercept=c(3,3), linetype="dashed", color="gray52")+
  facet_grid(rate_name ~ .,scale="free_y")+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        strip.text.y = element_text(size = 9, colour = "black",face="bold",angle=90),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=8),
        axis.text.x = element_text(colour = "black",size=8),
        axis.line=element_line(linetype="solid"),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.title=element_text(size=8),
        axis.title.x=element_text(face="bold"),
        plot.title = element_text(size=9),
        legend.position=('none'))+
  labs(x = "Minimum Read Support Threshold", y="",title=ID)+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0),limits=c(0,31))
a

#flip y-axis facet labels
a<- ggplotGrob(a)
a$layout[a$layout$name == "strip-right",c("l", "r")] <- 2
plot(a)
grid.draw(a)

filename <- paste("rs_temp",fileID, sep ="_")
filename <- paste(filename,"tiff", sep =".")
ggsave(filename,dpi=300, width=7.5,height=3.5,units="in")

#Read Support Cutoff for TELOCATE
summarydata <- read.table("BEDCOMPARE_MEANS.txt",header=TRUE)

TELOCATE <- summarydata[(summarydata$M2=="new_CT_run_3_N2_telocate_nonredundant.bed_F")& summarydata$fam=="family_aware", ]
SD_TELOCATE <- summarydata[(summarydata$M2=="new_CT_run_3_N2_telocate_nonredundant.bed_F_error")& summarydata$fam=="family_aware", ]
merged_TELOCATE<-merge(TELOCATE,SD_TELOCATE, by="DistanceCutoff")

merged_ALL<-merged_TELOCATE

#collapse 2 columns to prep for faceting
names(merged_ALL)
merged<-merged_ALL %>%
  gather(rate_name, rate_value, TPR.x:FDR.x)
#gather(rate_name, rate_value, TPR.x:FDR.x, FNR.x)
#separate TPR and FDR and associate with the proper error number
names(merged)
TPR<- merged[(merged$rate_name=="TPR.x"),c(1:7,8,10:11)]
FDR<- merged[(merged$rate_name=="FDR.x"),c(1:7,9,10:11)]
#FNR<- merged[(merged$rate_name=="FNR.x"),c(1:7,10:12)]

#rename TPR.y and FDR.y to error
#rename TPR.y and FDR.y to error
colnames(TPR)[8] <- "error"
colnames(FDR)[8] <- "error"

names(TPR)
names(FDR)
final <-rbind(TPR,FDR)
#final <-rbind(TPR,FDR,FNR)
names(final)
levels(final$rate_name)

final$rate_name <- factor(final$rate_name,
                          levels = c("TPR.x", "FDR.x"),
                          labels = c("TPR", "FDR"))

m <- ggplot(data = final, aes(x = DistanceCutoff, y = rate_value))
m <- m + geom_line(color="black")+ xlim(0,31)+
  geom_errorbar(aes(ymin=rate_value-error, ymax=rate_value+error),color="black") +
  geom_point(data = subset(final, rate_name=="TPR"),aes(x=0, y=100),alpha=0)+
  geom_point(data = subset(final, rate_name=="FDR"),aes(x=0, y=0),alpha=0)+
  geom_vline(xintercept=c(3,3), linetype="dashed", color="gray52")+
  facet_grid(rate_name ~ .,scale="free_y")+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        strip.text.y = element_text(size = 9, colour = "black",face="bold",angle=90),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=8),
        axis.text.x = element_text(colour = "black",size=8),
        axis.line=element_line(linetype="solid"),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.title=element_text(size=8),
        axis.title.x=element_text(face="bold"),
        plot.title = element_text(size=9),
        legend.position=('none'))+
  labs(x = "Minimum Read Support Threshold", y="",title=ID)+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0),limits=c(0,31))
m

#flip y-axis facet labels
m <- ggplotGrob(m)
m$layout[m$layout$name == "strip-right",c("l", "r")] <- 2
plot(m)
grid.draw(m)

filename <- paste("rs_telocate",fileID, sep ="_")
filename <- paste(filename,"tiff", sep =".")
ggsave(filename,dpi=300, width=7.5,height=3.5,units="in")


if (fileID=="RSV_SIM"){first<-a;m_first<-m}
if (fileID=="RSV_SIM_half"){second<-a;m_second<-m}
if (fileID=="RSV_SIM_quarter"){third<-a;m_third<-m}
if (fileID=="RSV_SIM_tenth"){fourth<-a;m_fourth<-m}
}

setwd("/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/")
a_all<-plot_grid(first,second,third,fourth,labels=c("A","B","C","D"))
m_all<-plot_grid(m_first,m_second,m_third,m_fourth,labels=c("A","B","C","D"))


setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(a_all,filename="Combined_Temp.tiff",dpi=300, width=7.5,height=5,units="in")
ggsave(m_all,filename="Combined_Telocate.tiff",dpi=300, width=7.5,height=5,units="in")


ggsave(a_all,filename="Combined_Temp.png",dpi=300, width=7.5,height=5,units="in")
ggsave(m_all,filename="Combined_Telocate.png",dpi=300, width=7.5,height=5,units="in")
