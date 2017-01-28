#!/usr/bin/R
# this script plots:
# 1) the TPR and FDR vs the distance cutoff
# 2) the TPR and FDR vs the read support threshold
# 3) the TPR and FDR vs the popFreq support threshold
# USE: graph_BEDCOMPARE_MEAN_rounds20-23.R
# NOTE: BEDCOMPARE_MEANS_ <reg|rs|vs> were named manually after running the filtering loop on the simulations 


library(tidyr)
library(dplyr)
library(ggplot2)
library(gtable)
library(cowplot)


####DISTANCE

file_list=c("/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/feb/round20",
            "/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/feb/round21",
            "/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/feb/round22",
            "/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/feb/round23")
i="/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/feb/round20"
#loop through levels of coverage simulations
for (i in file_list){
  print(i)
  setwd(i)
  fileID<-basename(i)
  print(file)
  
 
  if (fileID=="round20"){ID<-"130x Depth of Coverage"}
  if (fileID=="round23"){ID<-"65x Depth of Coverage"}
  if (fileID=="round22"){ID<-"32.5x Depth of Coverage"}
  if (fileID=="round21"){ID<-"13x Depth of Coverage"}
  

summarydata <- read.table("BEDCOMPARE_MEANS_reg.txt",header=TRUE)
names(summarydata)

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
#FNR<- merged[(merged$rate_name=="FNR.x"),c(1:7,10:12)]

#rename TPR.y and FDR.y to error
#rename TPR.y and FDR.y to error
colnames(TPR)[8] <- "error"
colnames(FDR)[8] <- "error"
#colnames(FNR)[8] <- "error"
names(TPR)
names(FDR)
final <-rbind(TPR,FDR)
#final <-rbind(TPR,FDR,FNR)
names(final)
levels(final$rate_name)

final$rate_name <- factor(final$rate_name,
                         levels = c("TPR.x", "FDR.x"),
                         labels = c("TPR", "FDR"))
#Distance Cuttoff
a <- ggplot(data = final, aes(x = DistanceCutoff, y = rate_value))
a <- a + geom_line(color="black")+ xlim(0,51)+
  geom_errorbar(aes(ymin=rate_value-error, ymax=rate_value+error),color="black") +
  geom_point(data = subset(final, rate_name=="TPR"),aes(x=0, y=100),alpha=0)+
  geom_point(data = subset(final, rate_name=="FDR"),aes(x=0, y=0),alpha=0)+
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
        axis.title=element_text(size=9),
        axis.title.x=element_text(face="bold"),
        axis.title.y=element_text(vjust=1.75),
        plot.title = element_text(size=9),
        legend.position=('none'))+
  labs(x = "Distance Cutoff (No.of Reads)", y="",title=ID)+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0),limits=c(0,50.5))
a

a <- ggplotGrob(a)
a$layout[a$layout$name == "strip-right",c("l", "r")] <- 2
plot(a)
grid.draw(a)


filename <- paste("ds",fileID, sep ="_")
filename <- paste(filename,"tiff", sep =".")
ggsave(filename,dpi=300, width=7.5,height=3.5,units="in")

##Read Support Cutoff
summarydata <- read.table("BEDCOMPARE_MEANS_rs.txt",header=TRUE)

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
#FNR<- merged[(merged$rate_name=="FNR.x"),c(1:7,10:12)]

#rename TPR.y and FDR.y to error
#rename TPR.y and FDR.y to error
colnames(TPR)[8] <- "error"
colnames(FDR)[8] <- "error"
#colnames(FNR)[8] <- "error"
names(TPR)
names(FDR)
final <-rbind(TPR,FDR)
#final <-rbind(TPR,FDR,FNR)
names(final)
levels(final$rate_name)

final$rate_name <- factor(final$rate_name,
                          levels = c("TPR.x", "FDR.x"),
                          labels = c("TPR", "FDR"))

b <- ggplot(data = final, aes(x = DistanceCutoff, y = rate_value))
b <- b + geom_line(color="black")+ xlim(0,51)+
  geom_errorbar(aes(ymin=rate_value-error, ymax=rate_value+error),color="black") +
  geom_point(data = subset(final, rate_name=="TPR"),aes(x=0, y=100),alpha=0)+
  geom_point(data = subset(final, rate_name=="FDR"),aes(x=0, y=0),alpha=0)+
  geom_vline(xintercept=c(8,8), linetype="dashed", color="gray52")+
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
        axis.title=element_text(size=9),
        axis.title.y=element_text(vjust=1.75),
        axis.title.x=element_text(face="bold"),
        plot.title = element_text(size=9),
        legend.position=('none'))+
  labs(x = "Minimum Read Support Threshold", y="",title=ID)+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0),limits=c(0,31))
b

b <- ggplotGrob(b)
b$layout[b$layout$name == "strip-right",c("l", "r")] <- 2
plot(b)
grid.draw(b)

filename <- paste("rs",fileID, sep ="_")
filename <- paste(filename,"tiff", sep =".")
ggsave(filename,dpi=300, width=7.5,height=3.5,units="in")

#PopFreq Cutoff
summarydata <- read.table("BEDCOMPARE_MEANS_vs.txt",header=TRUE)

TEMP <- summarydata[(summarydata$M2=="new_CT_run_3_N2_temp_nonredundant.bed_F")& summarydata$fam=="family_aware", ]
SD_TEMP <- summarydata[(summarydata$M2=="new_CT_run_3_N2_temp_nonredundant.bed_F_error")& summarydata$fam=="family_aware", ]
merged_TEMP<-merge(TEMP,SD_TEMP, by="DistanceCutoff")
names(TEMP)

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
#FNR<- merged[(merged$rate_name=="FNR.x"),c(1:7,10:12)]

#rename TPR.y and FDR.y to error
#rename TPR.y and FDR.y to error
colnames(TPR)[8] <- "error"
colnames(FDR)[8] <- "error"
#colnames(FNR)[8] <- "error"
names(TPR)
names(FDR)
final <-rbind(TPR,FDR)
#final <-rbind(TPR,FDR,FNR)
names(final)
levels(final$rate_name)

final$rate_name <- factor(final$rate_name,
                          levels = c("TPR.x", "FDR.x"),
                          labels = c("TPR", "FDR"))

c <- ggplot(data = final, aes(x = DistanceCutoff, y = rate_value))
c <- c + geom_line(color="black")+ xlim(0,51)+
  geom_errorbar(aes(ymin=rate_value-error, ymax=rate_value+error),color="black") +
  geom_point(data = subset(final, rate_name=="TPR"),aes(x=0, y=100),alpha=0)+
  geom_point(data = subset(final, rate_name=="FDR"),aes(x=0, y=0),alpha=0)+
  geom_vline(xintercept=c(.25,.25), linetype="dashed", color="gray52")+
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
        axis.title=element_text(size=9),
        axis.title.y=element_text(vjust=1.75),
        axis.title.x=element_text(face="bold"),
        plot.title = element_text(size=9),
        legend.position=('none'))+
  labs(x = "Minimum Population Frequency\nSupport Threshold", y="",title=ID)+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0),limits=c(0,1))
c

c <- ggplotGrob(c)
c$layout[c$layout$name == "strip-right",c("l", "r")] <- 2
plot(c)
grid.draw(c)

filename <- paste("vs",fileID, sep ="_")
filename <- paste(filename,"tiff", sep =".")
ggsave(filename,dpi=300, width=7.5,height=3.5,units="in",title=ID)


if (fileID=="round20"){first<-a;b_first<-b;c_first<-c}
if (fileID=="round23"){second<-a;b_second<-b;c_second<-c}
if (fileID=="round22"){third<-a;b_third<-b;c_third<-c}
if (fileID=="round21"){fourth<-a;b_fourth<-b;c_fourth<-c}
}


setwd("/Users/kristen/Documents/transposon_figure_data/simulations/depth_TPFD_files/feb/")
a_all<-plot_grid(first,second,third,fourth, labels=c("A","B","C","D"))
b_all<-plot_grid(b_first,b_second,b_third,b_fourth,labels=c("A","B","C","D"))
c_all<-plot_grid(c_first,c_second,c_third,c_fourth,labels=c("A","B","C","D"))


setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(a_all,filename="Combined_Insertion_DS.tiff",dpi=300, width=7.5,height=5,units="in")
ggsave(b_all,filename="Combined_Insertion_RS.tiff",dpi=300, width=7.5,height=5,units="in")
ggsave(c_all,filename="Combined_Insertion_PF.tiff",dpi=300, width=7.5,height=5,units="in")

ggsave(a_all,filename="Combined_Insertion_DS.png",dpi=300, width=7.5,height=5,units="in")
ggsave(b_all,filename="Combined_Insertion_RS.png",dpi=300, width=7.5,height=5,units="in")
ggsave(c_all,filename="Combined_Insertion_PF.png",dpi=300, width=7.5,height=5,units="in")