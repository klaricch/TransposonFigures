#!/usr/bin/R
# this script plots the distrubition of transposons (totals across all samples) according to their chromosome positions
# USE: TE_density.R

library(ggplot2)
library(grid)
library(dplyr)

setwd("/Users/kristen/Documents/transposon_figure_data/data")
summarydata <- read.table("CtCp_all_nonredundant.txt",header=TRUE)
names(summarydata)
names(summarydata)<-c("chr","start","end","TE","orientation","method","strain","class")

#3X-BIN .25MB
method_names <- list(
  'absent'="Absence",
  'new'="Insertion",
  'reference'="Reference"
)

method_labeller <- function(variable,value){
  if (variable=='method') {
    return(method_names[value])
  }else {
    return(as.character(value))
  }
}


summarydata <- distinct(summarydata, chr,start,method, orientation,class)

# Add y coordinates for "phantom" points
names(summarydata)
summarydata$top <- NA
summarydata$top[summarydata$method=="absent"] <- 6
summarydata$top[summarydata$method=="new"] <- 30
summarydata$top[summarydata$method=="reference"] <- 8
levels(summarydata$class)

#revalue classes
summarydata$class <- factor(summarydata$class,
                          levels = c("dnatransposon", "retrotransposon","unknown"),
                          labels = c("DNA Transposon", "Retrotransposon", "Unknown"))




m <- ggplot(summarydata, aes(x=start/1e6,fill=class))
m <-m + geom_bar(binwidth=.25)+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))+
  facet_grid(method ~ chr,scale="free",space = "free_x",labeller=method_labeller)+
  geom_point(data = subset(summarydata, method=="absent"),aes(y=top),alpha=0) +
  geom_point(data = subset(summarydata, method=="new"),aes(y=top),alpha=0) +
  geom_point(data = subset(summarydata, method=="reference"),aes(y=top),alpha=0) +
  
  labs(x="Chromosome Position (Mb)", y="Number of Transposition Events")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 9, colour = "black",face="bold"),
        panel.margin = unit(.25, "lines"),
        panel.border = element_rect(fill=NA,colour = "black"),
        panel.background = element_blank(),
        axis.title=element_text(size=9,face="bold"),
        axis.text.y = element_text(colour = "black",size=9),
        axis.text.x=element_blank(),
        #axis.text.x = element_text(colour = "black",size=9),
        axis.ticks =element_line(colour = "black"),
        legend.title=element_blank(),
        legend.position="bottom",
        legend.key.size=unit(.1,"cm"),
        legend.text=element_text(size=8))+
  scale_fill_manual(values = c("navy", "brown3", "darkgoldenrod2"))
m
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="Chromosome_Distribution.tiff",
       dpi=300,
       width=7.5,
       height=3.5,
       units="in")







#ARMS  AND  CENTERS (in bp)

I_A<-c(1:3858000,11040001:15072434)
I_C<-(3858001:11040000)

II_A<-c(1:4879000,12020001:15279421)
II_C<-(4879001:12020000)

III_A<-c(1:3722000,10340001:13783801)
III_C<-(3722001:10340000)

IV_A<-c(1:3896000,12970001:17493829)
IV_C<-(3896001:12970000)

V_A<-c(1:5897000,16550001:20924180)
V_C<-(5897001:16550000)

X_A<-c(1:6137000,12480001:17718942)
X_C<-(6137001:12480000)


#total the number of bases on arms,centers, and overall
arm_bases<-length(c(I_A,II_A,III_A,IV_A,V_A))
center_bases<-length(c(I_C,II_C,III_C,IV_C,V_C))
total_bases<-arm_bases+center_bases

#calculate probability of base on arms or centers
prob_arms <- arm_bases/total_bases
prob_center <- center_bases/total_bases

summarydata<-mutate(summarydata, arm=paste(chr,"A",sep="_"))
summarydata<-mutate(summarydata, center=paste(chr,"C",sep="_"))

one<-filter(summarydata, chr=="I")
two<-filter(summarydata, chr=="II")
three<-filter(summarydata, chr=="III")
four<-filter(summarydata, chr=="IV")
five <-filter(summarydata, chr=="V")
six<-filter(summarydata, chr=="X")

one<-mutate(one, region=ifelse(start %in% I_A,"ARM",ifelse(start %in% I_C,"CENTER","tip")))
two<-mutate(two, region=ifelse(start %in% II_A,"ARM",ifelse(start %in% II_C,"CENTER","tip")))
three<-mutate(three, region=ifelse(start %in% III_A,"ARM",ifelse(start %in% III_C,"CENTER","tip")))
four<-mutate(four, region=ifelse(start %in% IV_A,"ARM",ifelse(start %in% IV_C,"CENTER","tip")))
five<-mutate(five, region=ifelse(start %in% V_A,"ARM",ifelse(start %in% V_C,"CENTER","tip")))
six<-mutate(six, region=ifelse(start %in% X_A,"ARM",ifelse(start %in% X_C,"CENTER","tip")))

summarydata<-rbind(one,two,three,four,five,six)


#summarydata<-mutate(summarydata, region=ifelse(start %in% get(arm),"ARM",ifelse(start %in% get(center),"CENTER","ERROR")))
errors<-filter(summarydata,region=="ERROR")

get(I_A)

#remove TEs on chromosome tips
AC<-filter(summarydata,region!="tip")
#remove X chromosome
AC<-filter(AC,chr!="X")

#generate chi square contingency tables
chi_total <- table(AC$region)
chi_class <- table(AC$class,AC$region)
chi_method <- table(AC$method,AC$region)

chi_DNA <- chi_class[1,]
chi_Retro <- chi_class[2,]
chi_Unknown <- chi_class[3,]
chi_absent <- chi_method[1,]
chi_new <- chi_method[2,]
chi_reference <- chi_method[3,]

# chi square test with the hypothesized probabilities
chisq.test(chi_total,p=c(prob_arms,prob_center))
chisq.test(chi_DNA,p=c(prob_arms,prob_center))
chisq.test(chi_Retro,p=c(prob_arms,prob_center))
chisq.test(chi_Unknown,p=c(prob_arms,prob_center))
chisq.test(chi_absent,p=c(prob_arms,prob_center))
chisq.test(chi_new,p=c(prob_arms,prob_center))
chisq.test(chi_reference,p=c(prob_arms,prob_center))




