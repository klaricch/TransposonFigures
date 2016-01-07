#!/usr/bin/R
## this script plots the TPR and FDR for each of the telocate  and temp results of RSVSIM
## USE: (navigate to directory with BEDCOMAPRE results files) graph_TFPN_distances_RSVSIM.R
# NOTE: script has been spread out to other scripts, kept for notes

#install.packages("tidyr")
#install.packages("cowplot")
#library(cowplot)
#install.packages("scales")
#library(scales)
#install.packages("ggplot2")
library(tidyr)
library(dplyr)
library(ggplot2)
library(gtable)

#library(devtools)
#install_github("hadley/scales")
#library(scales)
#install_github("hadley/ggplot2movies")
#library(ggplot2movies)
#install_github("hadley/ggplot2")
#install.packages("scales")
#library(scales)
#install.packages("ggplot2movies")
#library(ggplot2movies)
#install.packages("ggplot2")
#library(ggplot2)


#######################################################################################

#                                 TEMP, RETROSEQ, TELOCATE

#######################################################################################

directory = getwd()
setwd(directory)
##REMOVE BELOW LATER
setwd("/Users/kristen/Desktop/Fig_p/round19")
summarydata <- read.table("BEDCOMPARE_MEANS.txt",header=TRUE)
names(summarydata)
#print(summarydata)


#ALL <- summarydata[(summarydata$M2=="new_CT_run_3_N2_temp_nonredundant.bed_F"|summarydata$M2=="new_CT_run_3_N2_telocate_nonredundant.bed_F"|summarydata$M2=="new_CT_run_3_N2_retroseq_nonredundant.bed_F")& summarydata$fam=="family_aware", ]
#SD <- summarydata[(summarydata$M2=="new_CT_run_3_N2_temp_nonredundant.bed_F_error"|summarydata$M2=="new_CT_run_3_N2_telocate_nonredundant.bed_F_error"|summarydata$M2=="new_CT_run_3_N2_retroseq_nonredundant.bed_F_error")& summarydata$fam=="family_aware", ]

TEMP <- summarydata[(summarydata$M2=="new_CT_run_3_N2_temp_nonredundant.bed_F")& summarydata$fam=="family_aware", ]
SD_TEMP <- summarydata[(summarydata$M2=="new_CT_run_3_N2_temp_nonredundant.bed_F_error")& summarydata$fam=="family_aware", ]
merged_TEMP<-merge(TEMP,SD_TEMP, by="DistanceCutoff")
names(TEMP)
print(merged_TEMP)

RETROSEQ <- summarydata[(summarydata$M2=="new_CT_run_3_N2_retroseq_nonredundant.bed_F")& summarydata$fam=="family_aware", ]
SD_RETROSEQ <- summarydata[(summarydata$M2=="new_CT_run_3_N2_retroseq_nonredundant.bed_F_error")& summarydata$fam=="family_aware", ]
merged_RETROSEQ<-merge(RETROSEQ,SD_RETROSEQ, by="DistanceCutoff")
print(merged_RETROSEQ)

TELOCATE <- summarydata[(summarydata$M2=="new_CT_run_3_N2_telocate_nonredundant.bed_F")& summarydata$fam=="family_aware", ]
SD_TELOCATE <- summarydata[(summarydata$M2=="new_CT_run_3_N2_telocate_nonredundant.bed_F_error")& summarydata$fam=="family_aware", ]
merged_TELOCATE<-merge(TELOCATE,SD_TELOCATE, by="DistanceCutoff")

merged_ALL<-rbind(merged_TEMP,merged_TELOCATE,merged_RETROSEQ)

#collapse 2 columns to prep for faceting
names(merged_ALL)
names(TEMP)
names(merged)
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
names(FNR)
final <-rbind(TPR,FDR)
#final <-rbind(TPR,FDR,FNR)
names(final)
levels(final$rate_name)

method_names <- list(
  'TPR.x'="TPR",
  'FDR.x'="FDR"
)

method_labeller <- function(variable,value){
  if (variable=='rate_name') {
    return(method_names[value])
  }else {
    return(as.character(value))
  }
}
labeller()
method_labeller("rate_name", "FDR.x")
names(TEMP)
TELOCATE <- final[(final$M2.x=="new_CT_run_3_N2_telocate_nonredundant.bed_F"), ]
TEMP <- final[(final$M2.x=="new_CT_run_3_N2_temp_nonredundant.bed_F"), ]
RETROSEQ <- final[(final$M2.x=="new_CT_run_3_N2_retroseq_nonredundant.bed_F"), ]

#a <- ggplot(data = TELOCATE, aes(x = DistanceCutoff, y = rate_value))
#a <- a + geom_line()+
#  facet_grid(rate_name ~ .,scale="free_y")+
#  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0),limits=c(0,100))
#a
names(TEMP)
levels(TEMP$rate_name)
#TEMP
# TEST 
labels<-data.frame(rate_name = "FDR")
names(labels)

labels <- lapply(labels, function(values) {
  if (is.logical(values)) {
    as.integer(values) + 1
  }
  if (is.factor(values)) {
    as.character(values)
  }
  else if (is.numeric(values) && .as_character) {
    as.character(values)
  }
  else {
    values
  }
})
labels

a <- ggplot(data = TEMP, aes(x = DistanceCutoff, y = rate_value))
a <- a + geom_line(color="black")+ xlim(0,51)+
  geom_errorbar(aes(ymin=rate_value-error, ymax=rate_value+error),color="black") +
  facet_grid(rate_name ~ .,scale="free_y",labeller=labels)+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        strip.text.y = element_text(size = 9, colour = "black",face="bold",angle=90),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9),
        legend.position=('none'))+
  labs(x = "Distance Cutoff", y="")+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0),limits=c(0,100))
a

# TEST
#http://www.statmethods.net/input/valuelabels.html
TEMP$rate_name <- factor(TEMP$rate_name,
                    levels = c("TPR.x", "FDR.x"),
                    labels = c("TPR", "FDR"))

a <- ggplot(data = TEMP, aes(x = DistanceCutoff, y = rate_value))
  a <- a + geom_line(color="black")+ xlim(0,51)+
  geom_errorbar(aes(ymin=rate_value-error, ymax=rate_value+error),color="black") +
  facet_grid(rate_name ~ .,scale="free_y")+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        strip.text.y = element_text(size = 9, colour = "black",face="bold",angle=90),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9),
        legend.position=('none'))+
  labs(x = "Distance Cutoff", y="")+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0),limits=c(0,100))
a
g <- ggplotGrob(a)
g$layout[g$layout$name == "strip-right",c("l", "r")] <- 2
plot(g)
#grid.newpage()
grid.draw(g)
#g
ggsave(g,filename="FAMILY_AWARE_TEMP_TPR_and_FDR.tiff",dpi=300, width=7.5,height=3.5,units="in")


names(TEMP)
#TEMP ROC
a <- ggplot(data = merged_TEMP, aes(x = FDR.x, y = TPR.x))
a <- a + geom_line(color="black")+ xlim(0,51)+
  #geom_errorbar(aes(ymin=rate_value-error, ymax=rate_value+error),color="black") +
  #facet_grid(rate_name ~ .,scale="free_y",labeller=method_labeller)+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        strip.text.y = element_text(size = 9, colour = "black",face="bold",angle=90),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9),
        legend.position=('none'))+
  labs(x = "Distance Cutoff", y="")+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0),limits=c(0,100))
a

#RETROSEQ

a <- ggplot(data = RETROSEQ, aes(x = DistanceCutoff, y = rate_value))
a <- a + geom_line(color="black")+ xlim(0,51)+
  geom_errorbar(aes(ymin=rate_value-error, ymax=rate_value+error),color="black") +
  facet_grid(rate_name ~ .,scale="free_y",labeller=method_labeller)+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        strip.text.y = element_text(size = 9, colour = "black",face="bold",angle=90),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9),
        legend.position=('none'))+
  labs(x = "Distance Cutoff", y="")+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0),limits=c(0,100))
a

#RETROSEQ ROC

a <- ggplot(data = merged_RETROSEQ, aes(x = FDR.x, y = TPR.x))
a <- a + geom_point(color="black")+
  #geom_errorbar(aes(ymin=rate_value-error, ymax=rate_value+error),color="black") +
  #facet_grid(rate_name ~ .,scale="free_y",labeller=method_labeller)+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        strip.text.y = element_text(size = 9, colour = "black",face="bold",angle=90),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9),
        legend.position=('none'))+
  labs(x = "Distance Cutoff", y="")+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0),limits=c(0,100))
a



#TELOCATE
a <- ggplot(data = TELOCATE, aes(x = DistanceCutoff, y = rate_value))
a <- a + geom_line(color="black")+ xlim(0,51)+
  geom_errorbar(aes(ymin=rate_value-error, ymax=rate_value+error),color="black") +
  facet_grid(rate_name ~ .,scale="free_y",labeller=method_labeller)+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        strip.text.y = element_text(size = 9, colour = "black",face="bold",angle=90),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9),
        legend.position=('none'))+
  labs(x = "Distance Cutoff", y="")+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0),limits=c(0,300))
a

a
view(final)
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
directory = getwd()
setwd(directory)
##REMOVE BELOW LATER
setwd("/Users/kristen/Desktop/Fig_p")
summarydata <- read.table("BEDCOMPARE_MEANS.txt",header=TRUE)
#print(summarydata)

## CHECK NAME ORDER!!!!!!!
colnames(summarydata)[3] <- "Read_Support"

#######################################################################################

#                                 TEMP

#######################################################################################

FAM <- summarydata[ summarydata$fam=="family_aware", ]
TEMP <- summarydata[(summarydata$M2=="new_CT_run_3_N2_temp_nonredundant.bed_F")& summarydata$fam=="family_aware", ]
SD<- summarydata[(summarydata$M2=="new_CT_run_3_N2_temp_nonredundant.bed_F_error")& summarydata$fam=="family_aware",]

names(FAM)
merged<-merge(TEMP,SD, by="Read_Support")
#collapse 2 columns to prep for faceting
merged<-merged %>%
  gather(rate_name, rate_value, TPR.x:FDR.x)
#separate TPR and FDR and associate with the proper error number
TPR<- merged[(merged$rate_name=="TPR.x"),c(1:7,8,10:11)]
FDR<- merged[(merged$rate_name=="FDR.x"),c(1:7,9:11)]
#rename TPR.y and FDR.y to error
colnames(TPR)[8] <- "error"
colnames(FDR)[8] <- "error"
final <-merge(TPR, FDR, all=TRUE)
names(final)
method_names <- list(
  'TPR.x'="TPR",
  'FDR.x'="FDR"
)

method_labeller <- function(variable,value){
  if (variable=='method') {
    return(method_names[value])
  }else {
    return(as.character(value))
  }
}




#final$label.new <- factor(df$label, levels=sort(c(""," ",levels(df$label))))
#p <- ggplot(df, aes(x=x, y=y)) + geom_point() + 
#  facet_wrap(~ label.new, ncol=3,drop=FALSE)

#ggsave <- ggplot2::ggsave; body(ggsave) <- body(ggplot2::ggsave)[-2]
names(summarydata)
#
labeller()
a <- ggplot(data = final, aes(x = Read_Support, y=rate_value))
a <- a + geom_line(color="darkorange")+ xlim(0,51)+
  geom_errorbar(aes(ymin=rate_value-error, ymax=rate_value+error),color="black") +
  facet_grid(rate_name ~ .,scale="free_y", labeller=method_labeller)+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        strip.text.y = element_text(size = 9, colour = "black",face="bold",angle=90),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9),
        legend.position=('none'))+
        labs(x = "Minimum Read Support Threshold", y="")+
  geom_vline(xintercept=c(4,4), linetype="dashed", color="gray52")
a
?png
ggsave(filename="FAMILY_AWARE_TEMP_TPR_and_FDR.tiff",dpi=300, width=7.5,height=3.5,units="in")
        
  #legend=FALSE +

#a
#flip y-axis facet labels
g <- ggplotGrob(a)
 g$layout[g$layout$name == "strip-right",c("l", "r")] <- 2
plot(g)
#grid.newpage()
grid.draw(g)
#g
ggsave(g,filename="FAMILY_AWARE_TEMP_TPR_and_FDR.tiff",dpi=300, width=7.5,height=3.5,units="in")

names(final)
print(final$rate_value)
#CONTINUOUS Y AXIS:
a <- ggplot(data = final, aes(x = Read_Support, y=rate_value, colour=rate_name))
a <- a + geom_line()+ xlim(0,51)+
geom_errorbar(aes(ymin=rate_value-error, ymax=rate_value+error),color="black") +
  geom_vline(xintercept=c(4,4), linetype="dashed", color="gray52")+
  #facet_grid(rate_name ~ ., labeller=method_labeller) +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        strip.text.y = element_text(size = 9, colour = "black",face="bold",angle=90),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9),
        legend.position=('none'))+
  scale_color_manual(values = c("darkorange", "darkorange"))+
  labs(x = "Minimum Read Support Threshold", y="FDR                                              TPR")
a
ggsave(filename="FAMILY_AWARE_TEMP_TPR_and_FDR_continuous_y.tiff",dpi=300, width=7.5,height=3.5,units="in")


#######################################################################################

#                                 TELOCATE

#######################################################################################

TELOCATE <- summarydata[(summarydata$M2=="new_CT_run_3_N2_telocate_nonredundant.bed_F")& summarydata$fam=="family_aware", ]
SD<- summarydata[(summarydata$M2=="new_CT_run_3_N2_telocate_nonredundant.bed_F_error")& summarydata$fam=="family_aware",]

names(SD)
merged<-merge(TELOCATE,SD, by="Read_Support")
#collapse 2 columns to prep for faceting
merged<-merged %>%
  gather(rate_name, rate_value, TPR.x:FDR.x)
#separate TPR and FDR and associate with the proper error number
TPR<- merged[(merged$rate_name=="TPR.x"),c(1:7,8,10:11)]
FDR<- merged[(merged$rate_name=="FDR.x"),c(1:7,9:11)]
#rename TPR.y and FDR.y to error
colnames(TPR)[8] <- "error"
colnames(FDR)[8] <- "error"
final <-merge(TPR, FDR, all=TRUE)

method_names <- list(
  'TPR.x'="TPR",
  'FDR.x'="FDR"
)

method_labeller <- function(variable,value){
  if (variable=='rate_name') {
    return(method_names[value])
  }else {
    return(as.character(value))
  }
}
a <- ggplot(data = TPR, aes(x = Read_Support, y=rate_value))
a <- a + geom_line(color="slateblue1")+ xlim(0,51)+
  geom_errorbar(aes(ymin=rate_value-error, ymax=rate_value+error),color="black") +
  geom_vline(xintercept=c(3,3), linetype="dashed", color="gray52")+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9),
        legend.position=('none'))+
        #background_grid() +
        
  #legend=FALSE +

  labs(x = "Minimum Read Support Threshold", y="TPR") 

a
#ggdraw(switch_axis_position(a, axis = 'y')) 
ggsave(filename="FAMILY_AWARE_TELOCATE_TPR.tiff",dpi=300, width=7.5,height=3.5,units="in")


stop("Stopping....comment this line to run alternative code")




a <- ggplot(data = final, aes(x = Read_Support, y=rate_value))
a <- a + geom_line(color="slateblue1")+ xlim(0,50)+
  geom_errorbar(aes(ymin=rate_value-error, ymax=rate_value+error),color="black") +
  facet_grid(rate_name ~ .,scale="free_y",labeller=method_labeller) +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9),
        legend.position=('none'))+
  #legend=FALSE +
  labs(x = "Read Support", y="") 
a
ggsave(filename="FAMILY_AWARE_TELOCATE_TPR.tiff",dpi=300, width=7.5,height=3.5,units="in")



#TEMP TPR
TEMP <- summarydata[(summarydata$M2=="new_CT_run_3_N2_temp_nonredundant.bed_F"| summarydata$M2=="new_CT_run_3_N2_temp_nonredundant.bed_NF" )& summarydata$fam=="family_aware", ]
pdf(file = "FAMILY_AWARE_TEMP_TPR.pdf")
a <- ggplot(data = TEMP, aes(x = Read_Support, y = TPR,colour=M2))
a <- a + geom_line()+ xlim(0,50)+theme(axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))+scale_color_manual(labels=c("TEMP_filtered","TEMP_unfiltered"), values = c("red", "darkturquoise"))+labs(x = "Read Support")
a
dev.off()

#TELOCATE FDR
TELOCATE <- summarydata[(summarydata$M2=="new_CT_run_3_N2_telocate_nonredundant.bed_F"| summarydata$M2=="new_CT_run_3_N2_telocate_nonredundant.bed_NF" )& summarydata$fam=="family_aware", ]
pdf(file = "FAMILY_AWARE_TELOCATE_FDR.pdf")
a <- ggplot(data = TELOCATE, aes(x = Read_Support, y = FDR,colour=M2))
a <- a + geom_line()+ xlim(0,50)+theme(axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))+scale_color_manual(labels=c("TELOCATE_filtered","TELOCATE_unfiltered"), values = c("red", "darkturquoise"))+labs(x = "Read Support")
a
dev.off()
#TELOCATE TPR
TELOCATE <- summarydata[(summarydata$M2=="new_CT_run_3_N2_telocate_nonredundant.bed_F"| summarydata$M2=="new_CT_run_3_N2_telocate_nonredundant.bed_NF" )& summarydata$fam=="family_aware", ]
pdf(file = "FAMILY_AWARE_TELOCATE_TPR.pdf")
a <- ggplot(data = TELOCATE, aes(x = Read_Support, y = TPR,colour=M2))
a <- a + geom_line()+ xlim(0,50)+theme(axis.ticks =element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black"))+scale_color_manual(labels=c("TELOCATE_filtered","TELOCATE_unfiltered"), values = c("red", "darkturquoise"))+labs(x = "Read Support")
a
dev.off()


