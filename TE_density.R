#!/usr/bin/R
# this script plots the distrubition of transposons (totals across all samples) according to their chromosome positions
# USE: TE_density.R

library(ggplot2)
library(grid)
library(dplyr)
library(tidyr)

setwd("/Users/kristen/Documents/transposon_figure_data/data")
summarydata <- read.table("CtCp_all_nonredundant.txt")
singletons <- read.table("singletons.txt")


names(summarydata)
names(singletons)<-"te"
names(summarydata)<-c("chr","start","end","TE","orientation","method","strain","class")

summarydata<-mutate(summarydata,frequency=ifelse(TE %in% singletons$te, "singleton","non-singleton"))
#3X-BIN .25MB
summarydata <- distinct(summarydata, chr,start,method, orientation,class,.keep_all=TRUE)
test<-filter(summarydata,method=="reference")


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

#revalue methods
summarydata$method <- factor(summarydata$method,
                            levels = c("new","reference","absent"),
                            labels = c("Insertion", "Reference","Absence"))
summarydata<-filter(summarydata,method!="Absence")

#summarydata<-filter(summarydata,method=="Insertion")
m <- ggplot(summarydata, aes(x=start/1e6,fill=class))
m <-m + geom_histogram(binwidth=.25)+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))+
  facet_grid(. ~ chr,scale="free",space = "free_x")+
  labs(x="Chromosome Position (Mb)", y="Number of Sites")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 11, colour = "black",face="bold"),
        panel.border = element_rect(fill=NA, colour="black",size=1, linetype="solid"),
        panel.background = element_blank(),
        plot.margin=unit(c(.1,.1,0,.1), "cm"),
        panel.spacing = unit(.5, "cm"),
        axis.title=element_text(size=11,face="bold"),
        axis.ticks =element_line(colour = "black"),
        axis.text.x=element_text(colour="black",size=11),
        axis.text.y=element_text(colour="black",size=11),
        legend.title=element_blank(),
        legend.position="none",
        legend.key.size=unit(1,"cm"),
        legend.text=element_text(size=11))+
  scale_fill_manual(values = c("DNA Transposon" = "navy", "Retrotransposon"="brown3","Unknown"="darkgoldenrod2"))

m
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="Chromosome_Distribution.tiff", dpi=350, width=6.75, height=3.375, units="in")
ggsave(filename="Chromosome_Distribution.png", dpi=350, width=6.75, height=3.375, units="in")


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
arm_bases<-length(c(I_A,II_A,III_A,IV_A,V_A,X_A))
center_bases<-length(c(I_C,II_C,III_C,IV_C,V_C,X_C))
total_bases<-arm_bases+center_bases

#calculate probability of base on arms or centers
prob_arms <- arm_bases/total_bases
prob_center <- center_bases/total_bases



#for autosome and X
auto_arm_bases<-length(c(I_A,II_A,III_A,IV_A,V_A))
auto_center_bases<-length(c(I_C,II_C,III_C,IV_C,V_C))
auto_total_bases<-auto_arm_bases+auto_center_bases
auto_prob_arms <- auto_arm_bases/auto_total_bases
auto_prob_center <- auto_center_bases/auto_total_bases

arm_bases
center_bases
auto_center_bases
x_arm_bases<-length(c(X_A))
x_center_bases<-length(c(X_C))
x_total_bases<-x_arm_bases+x_center_bases
x_prob_arms <- x_arm_bases/x_total_bases
x_prob_center <- x_center_bases/x_total_bases



#cds_arms <- 15939156
#cds_centers <- 18463383
#cds_total <- 34402539
#promoter_arms <- 12593834
#promoter_centers <- 10786591
#promoter_total <- 23380425
#intron_arms <- 36751268
#intron_centers <- 26031590
#intron_total <- 62782858
#intergenic_arms <- 21497784
#intergenic_centers <- 15774816
#intergenic_total <- 37272600


cds_arms <- 11402987
cds_centers <- 11960466
cds_total <- 23363453
promoter_arms <- 11039036
promoter_centers <- 9261981
promoter_total <- 20301017
intron_arms <- 19885426
intron_centers <- 12307061
intron_total <- 32192487
utr_arms <- 1603044
utr_centers <- 1877748
utr_total <- 3480792
intergenic_arms <- 22261883
intergenic_centers <- 16766696
intergenic_total <- 39028579


prob_cds_arms <- cds_arms/cds_total
prob_cds_centers <- cds_centers/cds_total
prob_promoter_arms <- promoter_arms/promoter_total
prob_promoter_centers <- promoter_centers/promoter_total
prob_intron_arms <- intron_arms/intron_total
prob_intron_centers <- intron_centers/intron_total
prob_intergenic_arms <- intergenic_arms/intergenic_total
prob_intergenic_centers <- intergenic_centers/intergenic_total
prob_utr_arms <- utr_arms/utr_total
prob_utr_centers <- utr_centers/utr_total


prob_cds<-cds_total/total_bases
prob_promoter<-promoter_total/total_bases
prob_intron<-intron_total/total_bases
prob_intergenic<-intergenic_total/total_bases
prob_utr<-utr_total/total_bases
#prob_genic
prob_utr

summarydata<-mutate(summarydata, arm=paste(chr,"A",sep="_"))
summarydata<-mutate(summarydata, center=paste(chr,"C",sep="_"))

head(summarydata)
test1<-filter(summarydata, method=="Insertion")
test2<-filter(summarydata, method=="Reference")
str(summarydata)


summarydata<-filter(summarydata, method=="Insertion")

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



#remove TEs on chromosome tips
AC<-filter(summarydata,region!="tip") # nothign should be filtered here
#remove X chromosome
#AC<-filter(AC,chr!="X")

#generate chi square contingency tables
chi_total <- table(AC$region)
chi_class <- table(AC$class,AC$region)
chi_method <- table(AC$method,AC$region)
chi_class
chi_method

table(summarydata$frequency)
#autosome and X
autosome<-filter(AC,chr!="X")
x_chrom<-filter(AC,chr=="X")

#single and non-single
autosome_single<-filter(AC,chr!="X",frequency=="singleton")
autosome_non_single<-filter(AC,chr!="X",frequency=="non-singleton")
x_chrom_single<-filter(AC,chr=="X",frequency=="singleton")
x_chrom_non_single<-filter(AC,chr=="X",frequency=="non-singleton")


chi_auto_total <- table(autosome$region)
chi_x_total <- table(x_chrom$region)

chi_auto_total_single <- table(autosome_single$region)
chi_auto_total_non_single <- table(autosome_non_single$region)
chi_x_total_single <- table(x_chrom_single$region)
chi_x_total_non_single <- table(x_chrom_non_single$region)
table(summarydata$method)
chi_auto_total
chi_x_total

chi_DNA <- chi_class[1,]
chi_Retro <- chi_class[2,]
chi_Unknown <- chi_class[3,]
chi_new <- chi_method[1,]
chi_reference <- chi_method[2,]
chi_absent<- chi_method[3,]

all <- data.frame(Factor=character(),
                  Arms=integer(),
                  Centers=integer(),
                    Chi2=numeric(),
                    PValue=integer(),
                    stringsAsFactors=FALSE)

all4 <- data.frame(Factor=character(),
                  Observed_Arms=integer(),
                  Observed_Centers=integer(),
                  Chi2=numeric(),
                  PValue=integer(),
                  Expected_Centers=integer(),
                  Expected_Arms=integer(),
                  Stat=character(),
                  stringsAsFactors=FALSE)

all2 <- data.frame(Factor=character(),
                  Observed=integer(),
                  Expected=integer(),
                  Chi2=numeric(),
                  PValue=integer(),
                  stringsAsFactors=FALSE)

all5 <- data.frame(Factor=character(),
                   Observed_Arms=integer(),
                   Observed_Centers=integer(),
                   Chi2=numeric(),
                   PValue=integer(),
                   Expected_Centers=integer(),
                   Expected_Arms=integer(),
                   Stat=character(),
                   stringsAsFactors=FALSE)
# chi square test with the hypothesized probabilities
ctest<-chisq.test(chi_total,p=c(prob_arms,prob_center))
all[1,]<-c("All Transposons",chi_total[1],chi_total[2],ctest$statistic,ctest$p.value)
#xx<-rbind(blank,ttest)
ctest<-chisq.test(chi_DNA,p=c(prob_arms,prob_center))
all[2,]<-c("DNA Transposons",chi_DNA[1],chi_DNA[2],ctest$statistic,ctest$p.value)
ctest<-chisq.test(chi_Retro,p=c(prob_arms,prob_center))
all[3,]<-c("Retrotransposons",chi_Retro[1],chi_Retro[2],ctest$statistic,ctest$p.value)
ctest<-chisq.test(chi_Unknown,p=c(prob_arms,prob_center))
all[4,]<-c("Unknown Transposons",chi_Unknown[1],chi_Unknown[2],ctest$statistic,ctest$p.value)
#chisq.test(chi_absent,p=c(prob_arms,prob_center))
chisq.test(chi_new,p=c(prob_arms,prob_center))
chisq.test(chi_reference,p=c(prob_arms,prob_center))


#Chi-Square tests per genomic feature
setwd("/Users/kristen/Documents/transposon_figure_data/data")
data <- read.table("essentiality_nonredundant_GO.txt",sep="\t",header=TRUE)
singleton_pos <- read.table("singleton_positions.txt",sep="\t",header=TRUE)
colnames(singleton_pos)<-"te"

#data<-mutate(data,position=paste(Chromosome,TE_start,sep="_"))
data<-filter(data, Method=="new")
data <- mutate(data, ID=paste(Chromosome, TE_start,sep="_"))
data<-mutate(data,frequency=ifelse(ID %in% singleton_pos$te, "singleton","non-singleton"))
data<-filter(data,frequency=="singleton")
#data<-filter(data,Chromosome!="X")
AC <- mutate(AC, ID=paste(chr,start,sep="_"))
ACS <- dplyr::select(AC, ID, region) #take subset of columns in AC
ACS<- distinct(ACS,.keep_all=TRUE)
data<-left_join(data, ACS, by = "ID")
nrow(data)
nrow(summarydata)

test<-distinct(data,Chromosome,TE_start)
##############################
# INCREASED vs DECREASED
##############################
#TEs on arms/# base pairs CDS on arms
#TEs on centers/ #base pairs CDS on centers
data$Region
levels(data$Region)[levels(data$Region)=="three_prime_UTR"] <- "UTR"
levels(data$Region)[levels(data$Region)=="five_prime_UTR"] <- "UTR"

data<-mutate(data,Region2=ifelse(Region=="intergenic","inter","genic"))
#generate chi square contingency tables
chi_feature <- table(data$Region,data$region)
chi_feature
chi_cds <- chi_feature[1,]
chi_utr <- chi_feature[3,]
chi_intergenic <- chi_feature[5,]
chi_intron <- chi_feature[6,]
chi_promoter <- chi_feature[7,]
chi_cds
chi_utr
chi_intron
chi_feature
chi_promoter
library(tibble)

data2<-distinct(data,Chromosome,TE_start)
num_tes<-nrow(data2)

chi_feats<-rbind(chi_cds,chi_promoter,chi_intergenic,chi_intron,chi_utr)
feats_df<-data.frame(chi_feats)
class(feats_df)
feats_df$total<-feats_df$ARM+feats_df$CENTER
feats_df<-dplyr::select(feats_df,total)
feats_df <- rownames_to_column(feats_df, "type")
#feats_df$non<-tes-feats_df$total

feats_df<-mutate(feats_df,non=num_tes-total)
feats_df$total<-as.numeric(feats_df$total)
feats_df$non<-as.numeric(feats_df$non)

varnames<-c("total","non")
Fcds<-feats_df[1,c(2,3)]
Fpromoter<-feats_df[2,c(2,3)]
Fintergenic<-feats_df[3,c(2,3)]
Fintron<-feats_df[4,c(2,3)]
Futr<-feats_df[5,c(2,3)]

ctable<-data.frame(rbind(Fcds,Fpromoter,Fintergenic,Fintron,Futr))
colnames(ctable)<-varnames
tes<-as.table(as.matrix(ctable))
tes


tes
chi2_cds <- tes[1,]
data_cds<-filter(data,Region=="CDS")
num_cds=nrow(data_cds)
ctest<-chisq.test(chi2_cds,p=c(prob_cds,1-prob_cds))
all2[1,]<-c("CDS",num_cds,prob_cds*num_tes,ctest$statistic,ctest$p.value)

chi2_promoter <- tes[2,]
data_promoter<-filter(data,Region=="promoter")
num_promoter=nrow(data_promoter)
ctest<-chisq.test(chi2_promoter,p=c(prob_promoter,1-prob_promoter))
all2[2,]<-c("Promoter",num_promoter,prob_promoter*num_tes,ctest$statistic,ctest$p.value)

chi2_intergenic <- tes[3,]
data_intergenic<-filter(data,Region=="intergenic")
num_intergenic=nrow(data_intergenic)
ctest<-chisq.test(chi2_intergenic,p=c(prob_intergenic,1-prob_intergenic))
all2[3,]<-c("Intergenic",num_intergenic,prob_intergenic*num_tes,ctest$statistic,ctest$p.value)

chi2_intron <- tes[4,]
data_intron<-filter(data,Region=="intron")
num_intron=nrow(data_intron)
ctest<-chisq.test(chi2_intron,p=c(prob_intron,1-prob_intron))
all2[4,]<-c("Intron",num_intron,prob_intron*num_tes,ctest$statistic,ctest$p.value)

chi2_utr <- tes[5,]
data_utr<-filter(data,Region=="UTR")
num_utr=nrow(data_utr)
ctest<-chisq.test(chi2_utr,p=c(prob_utr,1-prob_utr))
all2[5,]<-c("UTR",num_utr,prob_utr*num_tes,ctest$statistic,ctest$p.value)



data_auto=distinct(autosome,chr,start)
num_auto<-nrow(data_auto)
data_X=distinct(x_chrom,chr,start)
num_X<-nrow(data_X)


distinct_autosome_single<-distinct(autosome_single,chr,start)
num_auto_single<-nrow(distinct_autosome_single)

distinct_autosome_non_single<-distinct(autosome_non_single,chr,start)
num_auto_non_single<-nrow(distinct_autosome_non_single)


distinct_x_chrom_single<-distinct(x_chrom_single,chr,start)
num_x_chrom_single<-nrow(distinct_x_chrom_single)

distinct_x_chrom_non_single<-distinct(x_chrom_non_single,chr,start)
num_x_chrom_non_single<-nrow(distinct_x_chrom_non_single)

num_auto
auto_prob_arms
auto_prob_center
#autosome and X
ctest<-chisq.test(chi_auto_total,p=c(auto_prob_arms,auto_prob_center))
all2[6,]<-c("Autosome Arms",chi_auto_total[1],auto_prob_arms*num_auto,ctest$statistic,ctest$p.value)


ctest<-chisq.test(chi_auto_total,p=c(auto_prob_arms,auto_prob_center))
all2[7,]<-c("Autosome Centers",chi_auto_total[2],auto_prob_center*num_auto,ctest$statistic,ctest$p.value)
chi_auto_total

ctest<-chisq.test(chi_x_total,p=c(x_prob_arms,x_prob_center))
all2[8,]<-c("X Chromosome Arms",chi_x_total[1],x_prob_arms*num_X,ctest$statistic,ctest$p.value)


ctest<-chisq.test(chi_x_total,p=c(x_prob_arms,x_prob_center))
all2[9,]<-c("X Chromosome Center",chi_x_total[2],x_prob_center*num_X,ctest$statistic,ctest$p.value)




#SPLIT MORE HERE:

all4[1,]<-c("CDS",chi_cds[1],chi_cds[2],ctest$statistic,ctest$p.value,prob_cds_centers*num_cds,prob_cds_arms*num_cds)

chi_auto_total_single
ctest<-chisq.test(chi_auto_total_single,p=c(auto_prob_arms,auto_prob_center))
all5[1,]<-c("Sin.",chi_auto_total_single[1],chi_auto_total_single[2],ctest$statistic,ctest$p.value,auto_prob_center*num_auto_single,auto_prob_arms*num_auto_single,"Autosomes")

ctest<-chisq.test(chi_auto_total_non_single,p=c(auto_prob_arms,auto_prob_center))
all5[2,]<-c("Non-Sin.",chi_auto_total_non_single[1],chi_auto_total_non_single[2],ctest$statistic,ctest$p.value,auto_prob_center*num_auto_non_single,auto_prob_arms*num_auto_non_single,"Autosomes")

ctest<-chisq.test(chi_x_total_single,p=c(x_prob_arms,x_prob_center))
all5[3,]<-c("Sin.",chi_x_total_single[1],chi_x_total_single[2],ctest$statistic,ctest$p.value, x_prob_center*num_x_chrom_single,x_prob_arms*num_x_chrom_single,"X Chromosome")

ctest<-chisq.test(chi_x_total_non_single,p=c(x_prob_arms,x_prob_center))
all5[4,]<-c("Non-Sin.",chi_x_total_non_single[1],chi_x_total_non_single[2],ctest$statistic,ctest$p.value, x_prob_center*num_x_chrom_non_single,x_prob_arms*num_x_chrom_non_single, "X Chromosome")








chi_x_total
.libPaths()

all2$Expected<-as.numeric(all2$Expected)
all2$Observed<-as.numeric(all2$Observed)
all2<-mutate(all2, Change=ifelse(Observed>Expected,"Increased","Decreased"))
all2$Chi2<-as.numeric(all2$Chi2)
all2$PValue<-as.numeric(all2$PValue)
all2$Expected<-signif(all2$Expected,4)
all2$Chi2<-signif(all2$Chi2,4)
all2$PValue<-signif(all2$PValue,4)


colnames(all2)[1] <- "Region"
##############################
# CHI FEATURE
##############################

chi_feature2 <- table(data$Region2,data$region)
chi_genic<-chi_feature2[1,]
chi_genic
colnames(data)
# chi square test with the hypothesized probabilities
#chi_cds
#chi_x_total[2]
#num_X
#num_cds
#cds_prob_arms
ctest<-chisq.test(chi_cds,p=c(prob_cds_arms,prob_cds_centers))
all[5,]<-c("CDS Sin.",chi_cds[1],chi_cds[2],ctest$statistic,ctest$p.value)
all4[1,]<-c("Sin.",chi_cds[1],chi_cds[2],ctest$statistic,ctest$p.value,prob_cds_centers*num_cds,prob_cds_arms*num_cds,"CDS")


#ctest<-chisq.test(chi_cds_total_single,p=c(chi_prob_arms,chi_prob_center))
#all4[1,]<-c("CDS Singleton",chi_cds_total_single[1],chi_cds_total_single[2],ctest$statistic,ctest$p.value, cds_prob_center*num_cds_single,cds_prob_arms*num_single)

#all2[9,]<-c("X Chromosome Centers",chi_x_total[2],x_prob_center*num_X,ctest$statistic,ctest$p.value)

ctest<-chisq.test(chi_promoter,p=c(prob_promoter_arms,prob_promoter_centers))
all[6,]<-c("Promoter Sin.",chi_promoter[1],chi_promoter[2],ctest$statistic,ctest$p.value)
all4[2,]<-c("Sin.",chi_promoter[1],chi_promoter[2],ctest$statistic,ctest$p.value,prob_promoter_centers*num_promoter,prob_promoter_arms*num_promoter,"Promoter")

ctest<-chisq.test(chi_intergenic,p=c(prob_intergenic_arms,prob_intergenic_centers))
all[7,]<-c("Intergenic Sin.",chi_intergenic[1],chi_intergenic[2],ctest$statistic,ctest$p.value)
all4[3,]<-c("Sin.",chi_intergenic[1],chi_intergenic[2],ctest$statistic,ctest$p.value,prob_intergenic_centers*num_intergenic,prob_intergenic_arms*num_intergenic,"Intergenic")

ctest<-chisq.test(chi_intron ,p=c(prob_intron_arms,prob_intron_centers))
all[8,]<-c("Intron Sin.",chi_intron[1],chi_intron[2],ctest$statistic,ctest$p.value)
all4[4,]<-c("Sin.",chi_intron[1],chi_intron[2],ctest$statistic,ctest$p.value,prob_intron_centers*num_intron,prob_intron_arms*num_intron,"Intron")

ctest<-chisq.test(chi_utr ,p=c(prob_utr_arms,prob_utr_centers))
all[9,]<-c("UTR Sin.",chi_utr[1],chi_utr[2],ctest$statistic,ctest$p.value)
all4[5,]<-c("Sin.",chi_utr[1],chi_utr[2],ctest$statistic,ctest$p.value,prob_utr_centers*num_utr,prob_utr_arms*num_utr,"UTR")






########## REDO ABOVE WITH NON-SINGELTONS (lazy shortcut) #####################
########## REDO ABOVE WITH NON-SINGELTONS (lazy shortcut) #####################
########## REDO ABOVE WITH NON-SINGELTONS (lazy shortcut) #####################


#Chi-Square tests per genomic feature
setwd("/Users/kristen/Documents/transposon_figure_data/data")
data <- read.table("essentiality_nonredundant_GO.txt",sep="\t",header=TRUE)
singleton_pos <- read.table("singleton_positions.txt",sep="\t",header=TRUE)
colnames(singleton_pos)<-"te"

#data<-mutate(data,position=paste(Chromosome,TE_start,sep="_"))
data<-filter(data, Method=="new")
data <- mutate(data, ID=paste(Chromosome, TE_start,sep="_"))
data<-mutate(data,frequency=ifelse(ID %in% singleton_pos$te, "singleton","non-singleton"))
data<-filter(data,frequency=="non-singleton")
#data<-filter(data,Chromosome!="X")
AC <- mutate(AC, ID=paste(chr,start,sep="_"))
ACS <- dplyr::select(AC, ID, region) #take subset of columns in AC
ACS<- distinct(ACS,.keep_all=TRUE)
data<-left_join(data, ACS, by = "ID")
nrow(data)
nrow(summarydata)

test<-distinct(data,Chromosome,TE_start)
##############################
# INCREASED vs DECREASED
##############################
#TEs on arms/# base pairs CDS on arms
#TEs on centers/ #base pairs CDS on centers
data$Region
levels(data$Region)[levels(data$Region)=="three_prime_UTR"] <- "UTR"
levels(data$Region)[levels(data$Region)=="five_prime_UTR"] <- "UTR"

data<-mutate(data,Region2=ifelse(Region=="intergenic","inter","genic"))
#generate chi square contingency tables
chi_feature <- table(data$Region,data$region)
chi_feature
chi_cds <- chi_feature[1,]
chi_utr <- chi_feature[3,]
chi_intergenic <- chi_feature[5,]
chi_intron <- chi_feature[6,]
chi_promoter <- chi_feature[7,]
chi_cds
chi_utr
chi_intron
chi_feature
chi_promoter
library(tibble)

data2<-distinct(data,Chromosome,TE_start)
num_tes<-nrow(data2)

chi_feats<-rbind(chi_cds,chi_promoter,chi_intergenic,chi_intron,chi_utr)
feats_df<-data.frame(chi_feats)
class(feats_df)
feats_df$total<-feats_df$ARM+feats_df$CENTER
feats_df<-dplyr::select(feats_df,total)
feats_df <- rownames_to_column(feats_df, "type")
#feats_df$non<-tes-feats_df$total

feats_df<-mutate(feats_df,non=num_tes-total)
feats_df$total<-as.numeric(feats_df$total)
feats_df$non<-as.numeric(feats_df$non)

varnames<-c("total","non")
Fcds<-feats_df[1,c(2,3)]
Fpromoter<-feats_df[2,c(2,3)]
Fintergenic<-feats_df[3,c(2,3)]
Fintron<-feats_df[4,c(2,3)]
Futr<-feats_df[5,c(2,3)]

ctable<-data.frame(rbind(Fcds,Fpromoter,Fintergenic,Fintron,Futr))
colnames(ctable)<-varnames
tes<-as.table(as.matrix(ctable))
tes


tes
chi2_cds <- tes[1,]
data_cds<-filter(data,Region=="CDS")
num_cds=nrow(data_cds)
ctest<-chisq.test(chi2_cds,p=c(prob_cds,1-prob_cds))
#all2[1,]<-c("CDS",num_cds,prob_cds*num_tes,ctest$statistic,ctest$p.value)

chi2_promoter <- tes[2,]
data_promoter<-filter(data,Region=="promoter")
num_promoter=nrow(data_promoter)
ctest<-chisq.test(chi2_promoter,p=c(prob_promoter,1-prob_promoter))
#all2[2,]<-c("Promoter",num_promoter,prob_promoter*num_tes,ctest$statistic,ctest$p.value)

chi2_intergenic <- tes[3,]
data_intergenic<-filter(data,Region=="intergenic")
num_intergenic=nrow(data_intergenic)
ctest<-chisq.test(chi2_intergenic,p=c(prob_intergenic,1-prob_intergenic))
#all2[3,]<-c("Intergenic",num_intergenic,prob_intergenic*num_tes,ctest$statistic,ctest$p.value)

chi2_intron <- tes[4,]
data_intron<-filter(data,Region=="intron")
num_intron=nrow(data_intron)
ctest<-chisq.test(chi2_intron,p=c(prob_intron,1-prob_intron))
#all2[4,]<-c("Intron",num_intron,prob_intron*num_tes,ctest$statistic,ctest$p.value)

chi2_utr <- tes[5,]
data_utr<-filter(data,Region=="UTR")
num_utr=nrow(data_utr)
ctest<-chisq.test(chi2_utr,p=c(prob_utr,1-prob_utr))
#all2[5,]<-c("UTR",num_utr,prob_utr*num_tes,ctest$statistic,ctest$p.value)



data_auto=distinct(autosome,chr,start)
num_auto<-nrow(data_auto)
data_X=distinct(x_chrom,chr,start)
num_X<-nrow(data_X)

num_auto
auto_prob_arms
auto_prob_center
#autosome and X
ctest<-chisq.test(chi_auto_total,p=c(auto_prob_arms,auto_prob_center))
#all2[6,]<-c("Autosome Arms",chi_auto_total[1],auto_prob_arms*num_auto,ctest$statistic,ctest$p.value)


ctest<-chisq.test(chi_auto_total,p=c(auto_prob_arms,auto_prob_center))
#all2[7,]<-c("Autosome Centers",chi_auto_total[2],auto_prob_center*num_auto,ctest$statistic,ctest$p.value)
chi_auto_total

ctest<-chisq.test(chi_x_total,p=c(x_prob_arms,x_prob_center))
#all2[8,]<-c("X Chromosome Arms",chi_x_total[1],x_prob_arms*num_X,ctest$statistic,ctest$p.value)


ctest<-chisq.test(chi_x_total,p=c(x_prob_arms,x_prob_center))
#all2[9,]<-c("X Chromosome Centers",chi_x_total[2],x_prob_center*num_X,ctest$statistic,ctest$p.value)




#SPLIT MORE HERE:






all2$Expected<-as.numeric(all2$Expected)
all2$Observed<-as.numeric(all2$Observed)
all2<-mutate(all2, Change=ifelse(Observed>Expected,"Increased","Decreased"))
all2$Chi2<-as.numeric(all2$Chi2)
all2$PValue<-as.numeric(all2$PValue)
all2$Expected<-signif(all2$Expected,4)
all2$Chi2<-signif(all2$Chi2,4)
all2$PValue<-signif(all2$PValue,4)


colnames(all2)[1] <- "Region"
##############################
# CHI FEATURE
##############################

chi_feature2 <- table(data$Region2,data$region)
chi_genic<-chi_feature2[1,]
chi_genic
colnames(data)
# chi square test with the hypothesized probabilities
#chi_cds
#chi_x_total[2]
#num_X
#num_cds
#cds_prob_arms
ctest<-chisq.test(chi_cds,p=c(prob_cds_arms,prob_cds_centers))
#all[5,]<-c("CDS Singletons",chi_cds[1],chi_cds[2],ctest$statistic,ctest$p.value)
all4[6,]<-c("Non-Sin.",chi_cds[1],chi_cds[2],ctest$statistic,ctest$p.value,prob_cds_centers*num_cds,prob_cds_arms*num_cds,"CDS")


#ctest<-chisq.test(chi_cds_total_single,p=c(chi_prob_arms,chi_prob_center))
#all4[1,]<-c("CDS Singleton",chi_cds_total_single[1],chi_cds_total_single[2],ctest$statistic,ctest$p.value, cds_prob_center*num_cds_single,cds_prob_arms*num_single)

#all2[9,]<-c("X Chromosome Centers",chi_x_total[2],x_prob_center*num_X,ctest$statistic,ctest$p.value)

ctest<-chisq.test(chi_promoter,p=c(prob_promoter_arms,prob_promoter_centers))
#all[6,]<-c("Promoter Singletons",chi_promoter[1],chi_promoter[2],ctest$statistic,ctest$p.value)
all4[7,]<-c("Non-Sin.",chi_promoter[1],chi_promoter[2],ctest$statistic,ctest$p.value,prob_promoter_centers*num_promoter,prob_promoter_arms*num_promoter,"Promoter")

ctest<-chisq.test(chi_intergenic,p=c(prob_intergenic_arms,prob_intergenic_centers))
#all[7,]<-c("Intergenic Singletons",chi_intergenic[1],chi_intergenic[2],ctest$statistic,ctest$p.value)
all4[8,]<-c("Non-Sin.",chi_intergenic[1],chi_intergenic[2],ctest$statistic,ctest$p.value,prob_intergenic_centers*num_intergenic,prob_intergenic_arms*num_intergenic,"Intergenic")

ctest<-chisq.test(chi_intron ,p=c(prob_intron_arms,prob_intron_centers))
#all[8,]<-c("Intron Singletons",chi_intron[1],chi_intron[2],ctest$statistic,ctest$p.value)
all4[9,]<-c("Non-Sin.",chi_intron[1],chi_intron[2],ctest$statistic,ctest$p.value,prob_intron_centers*num_intron,prob_intron_arms*num_intron,"Intron")

ctest<-chisq.test(chi_utr ,p=c(prob_utr_arms,prob_utr_centers))
#all[9,]<-c("UTR Singletons",chi_utr[1],chi_utr[2],ctest$statistic,ctest$p.value)
all4[10,]<-c("Non-Sin.",chi_utr[1],chi_utr[2],ctest$statistic,ctest$p.value,prob_utr_centers*num_utr,prob_utr_arms*num_utr,"UTR")





###############################################################################
###############################################################################
###############################################################################



#ctest<-chisq.test(chi_genic,p=c(prob_genic_arms,prob_genic_center))
#all[10,]<-c("Genic",chi_genic[1],chi_genic[2],ctest$statistic,ctest$p.value)
#ctest<-chisq.test(chi_auto_total,p=c(auto_prob_arms,auto_prob_center))
#all[11,]<-c("Autosome",chi_auto_total[1],chi_auto_total[2],ctest$statistic,ctest$p.value)
#ctest<-chisq.test(chi_x_total,p=c(x_prob_arms,x_prob_center))
#all[12,]<-c("X Chromosome",chi_x_total[1],chi_x_total[2],ctest$statistic,ctest$p.value)

ctest
chi_x_total
auto_prob_arms
auto_prob_center

chi_promoter
str(all)
all$Chi2<-as.numeric(all$Chi2)
all$PValue<-as.numeric(all$PValue)
all$Chi2<-signif(all$Chi2,4)
all$PValue<-signif(all$PValue,4)

colnames(all)<-c("Factor","Number on Arms","Number on Centers","Chi-Squared Statistic", "P-Value")
setwd("/Users/kristen/Documents/transposon_figure_data/figures")

all
all2


all4$Expected_Arms<-as.numeric(all4$Expected_Arms)
all4$Expected_Centers<-as.numeric(all4$Expected_Centers)
all4$Observed_Arms<-as.numeric(all4$Observed_Arms)
all4$Observed_Centers<-as.numeric(all4$Observed_Centers)
#all4<-mutate(all4, Change=ifelse(Observed>Expected,"Increased","Decreased"))
all4$Chi2<-as.numeric(all4$Chi2)
all4$PValue<-as.numeric(all4$PValue)
#all4$Expected<-signif(all4$Expected,4)
all4$Chi2<-signif(all4$Chi2,4)
all4$PValue<-signif(all4$PValue,4)


m<-ggplot(all2, aes(Expected,Observed, label=rownames(all2)))  +
  geom_point()+
  geom_text(aes(x=Expected,y=Observed,label=Region), size=3,vjust=-.75)+
  scale_x_continuous(limits=c(-20,1350))+scale_y_continuous(limits=c(0,1800))+
  xlab("Expected Number of Transposons") +
  ylab("Observed Number of Transposons") 
m

all3<-gather(all2,cat,value,Expected,Observed)


all3$Region <- factor(all3$Region,
                          levels = c("Autosome Arms", "Autosome Centers","CDS","Intergenic","Intron","Promoter","UTR", "X Chromosome Arms","X Chromosome Centers"),
                          labels = c("Autosome\nArms", "Autosome\nCenters", "CDS","Intergenic","Intron","Promoter","UTR", "X Chromosome\nArms","X Chromosome\nCenters"))
b<-ggplot(all3, aes(x=Region,y=value,fill=cat))  +
  theme(axis.text.x=element_text(angle=35,hjust=1,size=11),
        axis.text.y=element_text(size=11),
        axis.title=element_text(size=11,face="bold"))+
  #scale_fill_manual 
  geom_bar(stat="identity",position="dodge")+
  scale_y_continuous(expand = c(0,0),limits=c(0,1700))+
  scale_fill_manual(values = c('Expected'="gray60",'Observed' = "gray17")) +	
  guides(fill=FALSE) +
  ylab("Number of Transposons")+xlab("")
b
ggsave(b,filename="Obs_v_Ex_TEsb.tiff",dpi=350, width=6.75,height=3.375,units="in")

all3
#all4$Factor <- factor(all4$Factor,
                     # levels = c("Intergenic Singletons", "Intergenic Non-Singletons",
                                # "Promoter Singletons","Promoter Non-Singletons",
                                 #"CDS Singletons", "CDS Non-Singletons",
                                 #"Intron Singletons", "Intron Non-Singletons",
                                 #"UTR Singletons","UTR Non-Singletons"),
                     # labels = c("Intergenic Singletons", "Intergenic Non-Singletons",
                                 #"Promoter Singletons","Promoter Non-Singletons",
                                 #"CDS Singletons", "CDS Non-Singletons",
                                 #"Intron Singletons", "Intron Non-Singletons",
                                 #"UTR Singletons","UTR Non-Singletons"))

all4$Factor <- factor(all4$Factor,
                      levels = c("Sin.", "Non-Sin."),
                      labels = c("Sin.", "Non-Sin."))

all5$Factor <- factor(all5$Factor,
                      levels = c("Sin.", "Non-Sin."),
                      labels = c("Sin.", "Non-Sin."))

all4<-mutate(all4,Significance=ifelse(PValue<0.05,"sig","non_sig"))
all5<-mutate(all5,Significance=ifelse(PValue<0.05,"sig","non_sig"))
all4_centers<-select(all4,-Observed_Arms,-Expected_Arms)
all4_arms<-select(all4,-Observed_Centers,-Expected_Centers)

all4_centers<-mutate(all4_centers,text_pos=Expected_Centers+.10*Expected_Centers)
all4_arms<-mutate(all4_arms,text_pos=Observed_Arms+.10*Observed_Arms)


all4_centers<-gather(all4_centers,cat,value,Expected_Centers,Observed_Centers)
all4_arms<-gather(all4_arms,cat,value,Expected_Arms,Observed_Arms)

#labels = c("Autosome\nArms", "Autosome\nCenters", "CDS","Intergenic","Intron","Promoter","UTR", "X Chromosome\nArms","X Chromosome\nCenters"))
c<-ggplot(all4_centers, aes(x=Factor,y=value,fill=cat))  +
  facet_grid(.~Stat) +
  theme(axis.text.x=element_text(angle=35,hjust=1,size=11),
        axis.text.y=element_text(size=11),
        axis.title=element_text(size=11,face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(angle=55, size = 9, colour = "black",face="bold"))+
  #scale_fill_manual 
  geom_bar(stat="identity",position="dodge") +
  scale_y_continuous(expand = c(0,0),limits=c(0,325))+
  scale_fill_manual(values = c('Expected_Centers'="gray60",'Observed_Centers' = "gray17")) +	
  guides(fill=FALSE) +
  ylab("Number of Transposons\non Centers")+xlab("")+
  geom_text(aes(x=Factor,y=text_pos,label=ifelse(Significance=="sig","*","")))
c


d<-ggplot(all4_arms, aes(x=Factor,y=value,fill=cat))  +
  facet_grid(.~Stat) +
  theme(axis.text.x=element_text(angle=35,hjust=1,size=11),
        axis.text.y=element_text(size=11),
        axis.title=element_text(size=11,face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(angle=55, size = 9, colour = "black",face="bold"))+
  #scale_fill_manual 
  geom_bar(stat="identity",position="dodge") +
 scale_y_continuous(expand = c(0,0),limits=c(0,550))+
  scale_fill_manual(values = c('Expected_Arms'="gray60",'Observed_Arms' = "gray17")) +	
  guides(fill=FALSE) +
  ylab("Number of Transposons\non Arms")+xlab("")+
  geom_text(aes(x=Factor,y=text_pos,label=ifelse(Significance=="sig","*","")))
d



all5$Expected_Arms<-as.numeric(all5$Expected_Arms)
all5$Expected_Centers<-as.numeric(all5$Expected_Centers)
all5$Observed_Arms<-as.numeric(all5$Observed_Arms)
all5$Observed_Centers<-as.numeric(all5$Observed_Centers)
#all4<-mutate(all4, Change=ifelse(Observed>Expected,"Increased","Decreased"))
all5$Chi2<-as.numeric(all5$Chi2)
all5$PValue<-as.numeric(all5$PValue)
#all4$Expected<-signif(all4$Expected,4)
#all5$Chi2<-signif(all5$Chi2,4)
#all5$PValue<-signif(all5$PValue,4)


all5_centers<-select(all5,-Observed_Arms,-Expected_Arms)
all5_arms<-select(all5,-Observed_Centers,-Expected_Centers)

all5_centers<-mutate(all5_centers,text_pos=Expected_Centers+.10*Expected_Centers)
all5_arms<-mutate(all5_arms,text_pos=Observed_Arms+.10*Observed_Arms)

all5_centers<-gather(all5_centers,cat,value,Expected_Centers,Observed_Centers)
all5_arms<-gather(all5_arms,cat,value,Expected_Arms,Observed_Arms)
#all5<-gather(all5,cat,value,Expected_Arms,Observed_Arms,Expected_Centers,Observed_Centers)


e<-ggplot(all5_centers, aes(x=Factor,y=value,fill=cat))  +
  facet_grid(.~Stat) +
  theme(axis.text.x=element_text(angle=35,hjust=1,size=11),
        axis.text.y=element_text(size=11),
        axis.title=element_text(size=11,face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(angle=55, size = 9, colour = "black",face="bold"))+
  #scale_fill_manual 
  geom_bar(stat="identity",position="dodge") +
  scale_y_continuous(expand = c(0,0),limits=c(0,850))+
  scale_fill_manual(values = c('Expected_Centers'="gray60",'Observed_Centers' = "gray17")) +	
  guides(fill=FALSE) +
  ylab("Number of Transposons\non Centers")+xlab("")+
  geom_text(aes(x=Factor,y=text_pos,label=ifelse(Significance=="sig","*","")))

e


f<-ggplot(all5_arms, aes(x=Factor,y=value,fill=cat))  +
  facet_grid(.~Stat) +
  theme(axis.text.x=element_text(angle=35,hjust=1,size=11),
        axis.text.y=element_text(size=11),
        axis.title=element_text(size=11,face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(angle=55, size = 9, colour = "black",face="bold"))+
  #scale_fill_manual 
  geom_bar(stat="identity",position="dodge") +
  scale_y_continuous(expand = c(0,0),limits=c(0,1050))+
  scale_fill_manual(values = c('Expected_Arms'="gray60",'Observed_Arms' = "gray17")) +	
  guides(fill=FALSE) +
  ylab("Number of Transposons\non Arms")+xlab("")+
  geom_text(aes(x=Factor,y=text_pos,label=ifelse(Significance=="sig","*","")))

f
library(cowplot)
ef<-plot_grid(e,f,c,d,ncol=2,nrow=2,labels=c('A', 'B','C','D'))
ef
ggsave(filename="ARM_CENTER_SINGLETON.tiff", dpi=350, width=6.75, height=7.5, units="in")
ggsave(filename="ARM_CENTER_SINGLETON.png", dpi=350, width=6.75, height=7.5, units="in")

all6<-rbind(all4,all5)

names(all6)[names(all6) == 'Factor'] <- 'Frequency'
names(all6)[names(all6) == 'Stat'] <- 'Region'
names(all6)[names(all6) == 'PValue'] <- 'P-Value'
names(all6)[names(all6) == 'Chi2'] <- 'Chi-Squared Statistic'

all6$Significance <- factor(all6$Significance,
                      levels = c("sig", "non_sig"),
                      labels = c("Sig.", "Not Sig."))


all6<-all6[,c(8,1,6,3,7,2,4,5,9)]
all6$Expected_Centers<-signif(all6$Expected_Centers,4)
all6$Expected_Arms<-signif(all6$Expected_Arms,4)


all<-filter(all, Factor %in% c("All Transposons", "DNA Transposons", "Retrotransposons","Unknown Transposons"))
write.table(all, file="Chi_Table.txt",sep="\t",quote=FALSE,row.names=FALSE)
write.table(all2, file="Chi_IncDec.txt",sep="\t",quote=FALSE,row.names=FALSE)
write.table(all6, file="Chi_AC_Sin.txt",sep="\t",quote=FALSE,row.names=FALSE)

###autosome arm c
###X chromosome arm c
#genic
####intergenic

###CDS
###intron
#utr
####promoter

#TAJIMA D
m <- ggplot(summarydata, aes(x=start/1e4))
m <-m + geom_histogram(binwidth=1)+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))+
  facet_grid(. ~ chr,scale="free",space = "free_x")+

  labs(x="Chromosome Position (Mb)", y="Number of Sites")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 9, colour = "black",face="bold"),
        panel.border = element_rect(fill=NA, colour="black",size=1, linetype="solid"),
        panel.background = element_blank(),
        plot.margin=unit(c(.1,.1,0,.1), "cm"),
        panel.spacing = unit(.5, "cm"),
        axis.title=element_text(size=9,face="bold"),
        axis.ticks =element_line(colour = "black"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        legend.title=element_blank(),
        legend.position="none",
        legend.key.size=unit(1,"cm"),
        legend.text=element_text(size=9))+
  #scale_fill_manual(values = c("navy", "brown3", "darkgoldenrod2"))
  scale_fill_manual(values = c("DNA Transposon" = "navy", "Retrotransposon"="brown3","Unknown"="darkgoldenrod2"))

m
a<-ggplot_build(m)
counts<-a$data[[1]]
head(counts)
counts<-dplyr::select(counts, x, count,PANEL)

colnames(counts)<-c("Bin","Count","CHROM")
counts$CHROM<-gsub("6","X",counts$CHROM)
counts$CHROM<-gsub("5","V",counts$CHROM)
counts$CHROM<-gsub("4","IV",counts$CHROM)
counts$CHROM<-gsub("3","III",counts$CHROM)
counts$CHROM<-gsub("2","II",counts$CHROM)
counts$CHROM<-gsub("1","I",counts$CHROM)

setwd("/Users/kristen/Documents/transposon_figure_data/data")
tj <- read.table("TD10.txt",header=TRUE)
tj<-dplyr::select(tj,CHROM,BIN_START,TajimaD)
counts<-mutate(counts,BIN_START=Bin*10000)

counts<-left_join(counts, tj, by = c("CHROM","BIN_START"))
counts<-mutate(counts,MEAN=mean(counts$Count))
counts<-mutate(counts,SD=sd(counts$Count))



#Find top %5 of values
n <- 5
top_five<-counts[counts$Count > quantile(counts$Count,prob=1-n/100),]

counts_table<-dplyr::select(top_five,CHROM, BIN_START,Count,TajimaD) 
test<-filter(counts_table,Count>=2)
threes<-filter(counts_table,Count>=3)


counts_table<-filter(counts_table,Count>=3,TajimaD>2.0 |TajimaD<(-2))
bed<-counts_table
counts_table<-mutate(counts_table,TajimaD=signif(TajimaD,4))

colnames(counts_table)<-c("Chromosome","Bin Start (bp)","Number of Transposons", "Tajima's D")
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
write.table(counts_table, file="TajimaD_Table.txt",sep="\t",quote=FALSE,row.names=FALSE)


###T test
#ROI<-filter(top_five,TajimaD>2.0 |TajimaD<(-2))

top_five<-mutate(top_five, score="TOPFIVE")
top_five<-dplyr::select(top_five,CHROM, BIN_START,score)
counts<-left_join(counts,top_five, by = c("CHROM","BIN_START"))
counts<-mutate(counts, SCORE=ifelse(!is.na(score),"TOP","BOTTOM"))
#T test on Tajima's D between top 5% and bottom 95%
counts<-mutate(counts,abs_taj=abs(TajimaD))
wilcox.test(counts$abs_taj~counts$SCORE)
#t.test(counts$abs_taj~counts$SCORE)aaa


top_five$Taji

max(counts$Count)

bb<- ggplot(counts, aes(x=TajimaD))
bb<-bb + geom_histogram(binwidth=.2)
bb


#arm vs center for tajima D




one<-filter(threes, CHROM=="I")
two<-filter(threes, CHROM=="II")
three<-filter(threes, CHROM=="III")
four<-filter(threes, CHROM=="IV")
five <-filter(threes, CHROM=="V")
six<-filter(threes, CHROM=="X")

one<-mutate(one, region=ifelse(BIN_START %in% I_A,"ARM",ifelse(BIN_START %in% I_C,"CENTER","tip")))
two<-mutate(two, region=ifelse(BIN_START %in% II_A,"ARM",ifelse(BIN_START %in% II_C,"CENTER","tip")))
three<-mutate(three, region=ifelse(BIN_START %in% III_A,"ARM",ifelse(BIN_START %in% III_C,"CENTER","tip")))
four<-mutate(four, region=ifelse(BIN_START %in% IV_A,"ARM",ifelse(BIN_START %in% IV_C,"CENTER","tip")))
five<-mutate(five, region=ifelse(BIN_START %in% V_A,"ARM",ifelse(BIN_START %in% V_C,"CENTER","tip")))
six<-mutate(six, region=ifelse(BIN_START %in% X_A,"ARM",ifelse(BIN_START %in% X_C,"CENTER","tip")))

threes<-rbind(one,two,three,four,five,six)

# basing arm/center on the start position of the bin
chi_threes <- table(threes$region)
chi_threes
ctest<-chisq.test(chi_threes,p=c(prob_arms,prob_center))
ctest


# make bedfile
bed<-mutate(bed,BIN_END=BIN_START+10001)
bed<-dplyr::select(bed,CHROM,BIN_START,BIN_END)
bed$blank1<-"."
bed$blank2<-"."
bed$blank3<-"."
write.table(bed, file="TajimaD_interest.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

