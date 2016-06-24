#!/usr/bin/R

#1. Correlations of absence with insertion and reference with insertion (square 4x4 inch, 16 point fonts, no title)
#2. Numbers of TE families, numbers of TEs, etc. I can probably get them from the results section, but I don't know if it is done enough for me to copy
#3. Manhattan plots for our good traits (9 inches wide, 2.5 inches high, 14 point font, no title, family on the side) - just like your paper figure, maybe just your paper figure
#4. Plot of TE locations across genome by type (9 inches wide, 4.5 inches high, 16 point font, no title)
#5. Plot of TE insertions into genomic features (9 inches wide, 4.5 inches high, 16 point font, no title)

#############################################################################
#############################################################################
#############################################################################
#5

library(ggplot2)
library(grid)
library(dplyr)
library(cowplot)

setwd("/Users/kristen/Documents/transposon_figure_data/data")
data <- read.table("essentiality_nonredundant_GO.txt",sep="\t",header=TRUE,stringsAsFactors = F)
data<-filter(data, Method=="new")

# simplify UTRs
data<-mutate(data, Region=ifelse(Region=="three_prime_UTR"|Region=="five_prime_UTR","UTR",Region))
data<-distinct(data,Chromosome, TE_start)

# simplify Biotypes
data<-mutate(data,final_bio=ifelse(Region=="intergenic","Intergenic",ifelse(Biotype=="pseudogene"|Biotype=="transposon_pseudogene","Pseudogene","Genic")))
#-split plot: A) intergenic, genic, pseudogene, B) CDS, promoter, intron
#-potential table with pseudogenes for loss of function caused by TE

a <- ggplot(data,aes(x=TE_start/1e6,fill=final_bio))
a <- a + geom_histogram(binwidth=.25)+
  facet_grid(.~Chromosome, scale="free", space="free_x")+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))+
  geom_point(aes(y=30), alpha=0)+
  labs(x="", y= "Count")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 16, colour = "black", face = "bold"),
        panel.margin = unit(.25, "lines"),
        panel.border = element_rect(fill=NA, colour="black",size=1, linetype="solid"),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.text=element_text(size=16),
        legend.text.align = 0,
        plot.margin=unit(c(.1,.1,-.5,.1), "cm"),
        axis.title = element_text(size=16,face="bold"),
        axis.text.y = element_text(colour="black", size=16,face="bold"),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(colour="black", size=11,face="bold"),
        axis.ticks = element_line(colour="black"),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"))+
  
  scale_fill_manual(values = c('Genic'="gray17",'Intergenic' = "gray60", "Pseudogene"="tan3"))
a


max_y<-ggplot_build(a)$panel$ranges[[1]]$y.range
max_y<-max_y[2]
a<- a + scale_y_continuous(expand = c(0,0),limits=c(0,max_y*1.075))
a


protein_coding<-filter(data,final_bio=="Genic", Biotype=="protein_coding")
protein_coding<-filter(protein_coding,Region!="exon")
protein_coding<-filter(protein_coding,Region!="gene")
protein_coding$Region <- factor(protein_coding$Region,
                                levels = c("promoter", "CDS","intron","UTR"),
                                labels = c("Promoter", "CDS","Intron","UTR"))

b <- ggplot(protein_coding,aes(x=TE_start/1e6,fill=Region))
b <- b + geom_histogram(binwidth=.25)+
  facet_grid(.~Chromosome, scale="free", space="free_x")+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))+
  geom_point(aes(y=25), alpha=0)+
  theme(strip.background = element_blank(), 
        strip.text = element_blank(),
        #strip.text = element_text(size = 11, colour = "black", face = "bold"),
        panel.margin = unit(.25, "lines"),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA, colour="black",size=1, linetype="solid"),
        legend.title = element_blank(),
        legend.text=element_text(size=16),
        legend.text.align = 0,
        #plot.margin=unit(c(-.5,.1,.1,.1), "cm"),
        axis.title = element_text(size=16,face="bold"),
        axis.text.y = element_text(colour="black", size=16,face="bold"),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(colour="black", size=11,face="bold"),
        axis.ticks = element_line(colour="black"),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"))+
  labs(x="Chromosome Position (Mb)", y= "Count")+
  scale_fill_manual(values = c('CDS'="orange", 'Intron' = "plum2", 'Promoter' = "cornflowerblue","UTR"="olivedrab3")) 
b
max_y<-ggplot_build(b)$panel$ranges[[1]]$y.range
max_y<-max_y[2]
b<- b + scale_y_continuous(expand = c(0,0),limits=c(0,max_y*1.075))
b

all<-plot_grid(a,b,ncol=1,align="v" )+ background_grid(major = "xy", minor = "none")
all
setwd("/Users/kristen/Documents/transposon_figure_data/fig_for_e2")
ggsave(filename="Genic_Features.tiff",
       dpi=300,
       width=9,
       height=4.5,
       units="in")


b<-b+theme(strip.background = element_blank(),
             strip.text = element_text(size = 16, colour = "black", face = "bold"))
ggsave(b,filename="Genic_Features_b.tiff",
       dpi=300,
       width=10,
       height=4,
       units="in")



#############################################################################
#############################################################################
#############################################################################
#4
library(ggplot2)
library(grid)
library(dplyr)

setwd("/Users/kristen/Documents/transposon_figure_data/data")
summarydata <- read.table("CtCp_all_nonredundant.txt")
names(summarydata)
names(summarydata)<-c("chr","start","end","TE","orientation","method","strain","class")

#3X-BIN .25MB
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

#revalue methods
summarydata$method <- factor(summarydata$method,
                             levels = c("new","reference","absent"),
                             labels = c("Insertion", "Reference","Absence"))



m <- ggplot(summarydata, aes(x=start/1e6,fill=class))
m <-m + geom_histogram(binwidth=.25)+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))+
  facet_grid(method ~ chr,scale="free",space = "free_x")+
  geom_point(data = subset(summarydata, method=="Absence"),aes(y=top),alpha=0) +
  geom_point(data = subset(summarydata, method=="Insertion"),aes(y=top),alpha=0) +
  geom_point(data = subset(summarydata, method=="Reference"),aes(y=top),alpha=0) +
  
  labs(x="Chromosome Position (Mb)", y="Number of Sites")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 16, colour = "black",face="bold"),
        #panel.margin = unit(.25, "lines"),
        panel.border = element_rect(fill=NA, colour="black",size=1, linetype="solid"),
        panel.background = element_blank(),
        panel.margin.y=unit(.75,"cm"),
        plot.margin=unit(c(.1,.1,0,.1), "cm"),
        #panel.margin = unit(.75, "cm"),
        #panel.margin = unit(c(.5,.5,.5,.5), "cm"),
        #panel.margin = unit(c(.5,.5,.5,.5), "cm"),
        axis.title=element_text(size=16,face="bold"),
        axis.text.y = element_text(colour = "black",size=16),
        axis.text.x=element_blank(),
        #axis.text.x = element_text(colour = "black",size=9),
        axis.ticks =element_line(colour = "black"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        legend.title=element_blank(),
        legend.position="none",
        legend.key.size=unit(1,"cm"),
        legend.text=element_text(size=16))+
  #scale_fill_manual(values = c("navy", "brown3", "darkgoldenrod2"))
  scale_fill_manual(values = c("DNA Transposon" = "navy", "Retrotransposon"="brown3","Unknown"="darkgoldenrod2"))

m <- m
m
setwd("/Users/kristen/Documents/transposon_figure_data/fig_for_e2")
ggsave(filename="Chromosome_Distribution.tiff",
       dpi=300,
       width=9,
       height=4.5,
       units="in")

#############################################################################
#############################################################################
#############################################################################
#3


library(dplyr)
library(ggplot2)
library(data.table)
library(grid)
library(stringr)
library(gridExtra)
library(tidyr)
library(scales)
library(gtable)
library(cowplot)

setwd("/Users/kristen/Documents/transposon_figure_data/data")
load("Processed_Transposon_Mappings_SUBSET2.Rda")


unique(processed_mapping_df$trait)
load("count_QTL.Rda")

# pull unique combos, remove strain column(don't need specific strain info at this point)
processed_mapping_df<- processed_mapping_df %>% distinct(trait,marker,strain)

#create family and method columns
processed_mapping_df$family <- stringr::str_split_fixed(processed_mapping_df$trait, "_TRANS_",2)[,2]
processed_mapping_df$method <- stringr::str_split_fixed(processed_mapping_df$trait, "_TRANS_",2)[,1]

#read in position data and create family column
positions <- read.table("CtCp_all_nonredundant.txt")
names(positions)<-c("CHROM","start","end","TE","orientation","method","strain","class")
positions$family<- stringr::str_split_fixed(positions$TE, regex("_(non-)?reference"),2)[,1]
positions$family<- paste(stringr::str_split_fixed(positions$family, "_",4)[,3],stringr::str_split_fixed(positions$family, "_",4)[,4],sep="_")
positions$family <- gsub("_$" ,"",positions$family)
positions$family <- gsub("_non-reference(.*)$" ,"",positions$family)

#select traits above BF.....this step not needed, double checking everything is above BF
selection<-filter(processed_mapping_df, log10p > BF)

#extract the count base traits
base_traits <-selection[(selection$method=="absent"| selection$method=="new" |selection$method=="reference"|selection$method=="ZERO_new"|selection$method=="ONE_new"), ]
counts<-subset(base_traits, grepl("_C$", base_traits$family))
counts$family <- gsub("_C$" ,"",counts$family)

processed_mapping_df <- distinct(select(processed_mapping_df, -strain,-allele,-value))
processed_mapping_df<- processed_mapping_df %>% distinct(trait,marker)

#pull out only position traits from mappings dataframe
position_traits<-subset(selection,
                        grepl('^I', selection$trait) |
                          grepl('^V', selection$trait) |
                          grepl('^X', selection$trait))

#create family column
position_traits$family  <- paste(stringr::str_split_fixed(position_traits$trait, "_",4)[,3],stringr::str_split_fixed(position_traits$trait, "_",4)[,4],sep="_")
position_traits$family <- gsub("_$" ,"",position_traits$family)
position_traits$family <- gsub("_non-reference(.*)$" ,"",position_traits$family)

# add position TRAIT_col family info to processed_mapping_df
processed_mapping_df<-processed_mapping_df %>%mutate(family = ifelse(processed_mapping_df$trait %in% position_traits$trait, (paste(stringr::str_split_fixed(processed_mapping_df$trait, "_",4)[,3],stringr::str_split_fixed(processed_mapping_df$trait, "_",4)[,4],sep="_")), processed_mapping_df$family))

selection<-counts

#strip count marker and remnant marks from dataframes
selection$trait <- gsub("_C$" ,"",selection$trait)
processed_mapping_df$trait <- gsub("_C$" ,"",processed_mapping_df$trait)
processed_mapping_df$family <- gsub("_C$" ,"",processed_mapping_df$family)
processed_mapping_df$family <- gsub("_$" ,"",processed_mapping_df$family)
processed_mapping_df$family <- gsub("_non-reference(.*)$" ,"",processed_mapping_df$family)

processed_mapping_df<-mutate(processed_mapping_df,ID=paste(trait,peak_id,sep="_"))
copy<-processed_mapping_df
processed_mapping_df<-mutate(processed_mapping_df,SNP_col=ifelse(ID %in% count_QTL$trait,"PASS","FAIL" ))

count_QTL<-mutate(count_QTL, trait2=gsub("_\\d+$","",trait)) 
selection <- filter(selection, (trait %in% count_QTL$trait2))

processed_mapping_df<-filter(processed_mapping_df,CHROM != "MtDNA")
class_subset<- positions %>% distinct(class,family) %>% select(class,family)
selection <-merge(selection, class_subset, by="family")
selection<-arrange(selection,class,family,method)

unique(selection$trait)
count<-0
for (i in unique(selection$trait)){
  specific_trait<- processed_mapping_df[processed_mapping_df$trait == i, ]
  empty <-specific_trait[specific_trait$method==NA,]
  #specific_trait_mx <- max(specific_trait$log10p)
  class_TE<-unique(filter(selection,trait==i)$class)
  pvalues<-filter(specific_trait,log10p !="Inf") #
  specific_trait_mx <- max(pvalues$log10p) #
  TE<-specific_trait$family[1]
  rect_data<-filter(specific_trait,SNP_col=="PASS")
  plot_title<-gsub(".*_TRANS_","",i)
  plot_title<-gsub("_CE$","",plot_title)
  plot_title<-gsub("WBTransposon","WBT",plot_title)
  
  A<- processed_mapping_df %>%
    filter(trait == i)%>%
    .[order(.$peak_id,na.last=FALSE),]%>% 
    ggplot(.)+
    aes(x=POS/1e6,y=log10p,fill=BF)+ #fill to get legend
    geom_rect(data=rect_data,mapping=aes(xmin=startPOS/1e6, xmax=endPOS/1e6, ymin=0, ymax= Inf),fill="thistle1", alpha=1) +
    geom_point(aes( color=ifelse(log10p> BF & SNP_col=="PASS", 'red', 'black')),size=1)+
    facet_grid(.~CHROM,scale="free_x",space = "free_x") + #scale_color_identity() +
    geom_hline(aes(yintercept=BF),color="grey60",linetype="dashed")+
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(color="black", size=0.5, linetype="solid", fill=NA),
          panel.margin = unit(.6, "lines"),
          panel.background = element_blank(),
          axis.ticks =element_line(colour = "black"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(colour = "black",size=14),
          axis.title.y = element_text(size=14,colour=ifelse(class_TE=="dnatransposon","navy",ifelse(class_TE=="retrotransposon","brown3","darkgoldenrod2"))),
          axis.title=element_text(size=14),
          plot.margin=unit(c(.05,.30,-.5,.30), "cm"),
          legend.title = element_text(size = 14, colour = ifelse(class_TE=="dnatransposon","navy",ifelse(class_TE=="retrotransposon","brown3","darkgoldenrod2")), angle = 270),
          legend.text = element_blank(),
          legend.key.size=unit(0,"cm"),
          legend.key = element_rect(colour = "pink"),
          legend.position=('right'))+
    labs(x="",y="",colour="black",size=16,face="bold")+
    scale_color_identity()+
    scale_fill_continuous(name=plot_title)+
    scale_y_continuous(breaks= pretty_breaks(),expand=c(0,0),limits=c(0,specific_trait_mx+.075*specific_trait_mx),labels = function(x) format(x,width = 4))
  A
  
  if (count==0){B<-A+theme(strip.background = element_rect(fill = "white"),
                           strip.text.x = element_text(size = 14, colour = "black",face="bold"));
                first<-B}
  if (count==1){second<-A}
  if (count==2){third<-A}
  if (count==3){fourth<-A}
  if (count==4){fifth<-A}
  if (count==5){sixth<-A}
  count<-count+1
}

a_all<-plot_grid(first,second,third,fourth,fifth,ncol=1) #ZER)_new_TRANS_NeSL-1_C no  longer in here so don't need fifth 
label<-expression(bold(-log["10"](p)))

a_all<- a_all + draw_label(label, x = .04, y = 0.5, hjust = .5, vjust = .5,
                           fontfamily = "", fontface = "bold", colour = "black", size = 14,
                           angle = 90, lineheight = 0.9, alpha = 1)

df <- data.frame(1,2)
blank_plot<-ggplot(df,aes(x=1,y=1)) + geom_point(color="white") + theme(axis.line=element_blank(),axis.text =element_blank(),axis.ticks =element_blank(),axis.title =element_blank(),panel.background = element_blank(),panel.grid = element_blank())
a_all<-plot_grid(a_all,blank_plot,ncol=1,rel_heights = c(1, .03))


a_all<- a_all + draw_label("Chromosome Position (Mb)", x = .5, y = 0.020, hjust = .5, vjust = .5,
                           fontfamily = "", fontface = "bold", colour = "black", size = 14,
                           angle = 0, lineheight = 0.9, alpha = 1)

setwd("/Users/kristen/Documents/transposon_figure_data/fig_for_e2")
ggsave(filename="five_trait_QTL.tiff",
       dpi=300,
       width=7.5,
       height=12.5,
       units="in")

fourth<-fourth+theme(strip.background = element_rect(fill = "white"),
                     plot.margin=unit(c(.05,.30,.1,.30), "cm"),
                     axis.title=element_text(size=16,colour="black",face="bold"),
                     strip.text.x = element_text(size = 14, colour = "black",face="bold"))+
                     labs(x="Chromosome Position (Mb)",y="",colour="black",size=16)

  fourth<-plot_grid(fourth) + draw_label(label, x = .04, y = 0.5, hjust = .5, vjust = .5,
             fontfamily = "", fontface = "bold", colour = "black", size = 16,
             angle = 90, lineheight = 0.9, alpha = 1)
  

fourth
ggsave(fourth,filename="fourth.tiff",
       dpi=300,
       width=7.5,
       height=2.5,
       units="in")


#############################################################################
#############################################################################
#############################################################################
#1

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(cowplot)
library(grid)
setwd("/Users/kristen/Documents/transposon_figure_data/data")
summarydata <- read.table("T_kin_C_matrix_full.txt",header=TRUE)
#remove ZERO_new traits
summarydata<-subset(summarydata,!grepl('^ZERO_new', summarydata$trait))
summarydata<-subset(summarydata,!grepl('^coverage', summarydata$trait))
#clean trait names
summarydata$trait <- gsub("_C$" ,"",summarydata$trait)
summarydata$trait <- gsub("^ONE_new" ,"new",summarydata$trait)

#new column that specifies what caller was used
summarydata$method<- stringr::str_split_fixed(summarydata$trait, "_TRANS_",2)[,1]
#new column that specifies TE family
summarydata$transposon<- stringr::str_split_fixed(summarydata$trait, "_TRANS_",2)[,2]
summarydata<-filter(summarydata,transposon=="total") # this will get total ins,ref,abs calls, NOT total DNA, Retro, Unknonwn
unique(summarydata$transposon)

#names(summarydata)
summarydata<-gather(summarydata, "sample","value",2:(ncol(summarydata)-2))
summarydata<-rename(summarydata,total_tes=value)

#reformat the data
total_absence<-filter(summarydata,method=="absent")
total_reference<-filter(summarydata,method=="reference")
total_insertion<-filter(summarydata,method=="new")

#SCATTER
final_merge<- Reduce(function(x, y) merge(x, y, all=TRUE,by="sample"), list(total_absence, total_reference, total_insertion))
names(final_merge)<-c("sample", "trait.x",  "method.x",	"transposon.x",	"total_absences",	"trait.y",	"method.y",	"transposon.y",	"total_references",	"trait",	"method",	"transposon",	"total_insertions")


#1 ABSENCE vs INSERTION
#spearman correlation
correlation<-cor.test(final_merge$total_absences, final_merge$total_insertions,method="spearman",exact=FALSE)
rho<-round(correlation$estimate,3)

max_insertions<-max(final_merge$total_insertions)
max_absences<-max(final_merge$total_absences)

la <- paste("italic(rho) == ", rho)
m1 <- ggplot(final_merge, aes(x=total_insertions, y=total_absences))
m1 <- m1 + geom_point(size=1.25) + xlim(0,max_insertions)+ ylim(0,max_insertions)+
  geom_smooth(method="lm",se=FALSE,col="red")+
  geom_abline(slope=1,linetype="dashed",colour="gray52")+
  annotate("text", x=.2*max_insertions, y=.9*max_insertions,label=la,parse=TRUE, colour="red",size=4.5)+
  theme(strip.text.x = element_text(size = 6, colour = "black"),
        strip.background = element_blank(),
        legend.position=c(.90,0.75),
        legend.background = element_rect(fill=FALSE),
        legend.text=element_text(size=16),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=16),
        axis.text.x = element_text(colour = "black",size=16),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.title=element_text(size=16,face="bold"))+
  guides(fill=FALSE) +
  labs(x = "Insertion Sites", y = "Absence Sites")
m1
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="Absence_vs_Insertion.tiff",
       dpi=300,
       width=4,
       height=4,
       units="in")

#3 INSERTION vs REFERENCE
#spearman correlation
correlation<-cor.test(final_merge$total_insertions, final_merge$total_references,method="spearman",exact=FALSE)
rho<-round(correlation$estimate,3)

max_references<-max(final_merge$total_references)
max_insertions<-max(final_merge$total_insertions)
la <- paste("italic(rho) == ", rho)

max_references<-max(final_merge$total_references)
m3 <- ggplot(final_merge, aes(x=total_references, y=total_insertions))
m3 <- m3 + geom_point(size=1.25) + xlim(0,max_references)+ ylim(0,max_references)+
  geom_smooth(method="lm",se=FALSE,col="red")+
  geom_abline(slope=1,linetype="dashed",colour="gray52")+
  annotate("text", x=.2*max_references, y=.9*max_references,label=la,parse=TRUE, colour="red",size=4.5)+
  theme(strip.text.x = element_text(size = 16, colour = "black"),
        strip.background = element_blank(),
        legend.position=c(.90,0.75),
        legend.background = element_rect(fill=FALSE),
        legend.text=element_text(size=16),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=16),
        axis.text.x = element_text(colour = "black",size=16),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.title=element_text(size=16,face="bold"))+
  guides(fill=FALSE) +
  labs(x = "Reference Sites", y = "Insertion Sites")

ggsave(filename="Insertion_vs_Reference.tiff",
       dpi=300,
       width=4,
       height=4,
       units="in")
setwd("/Users/kristen/Documents/transposon_figure_data/fig_for_e2")
plot_grid(m1, m3,ncol=2,labels=c('A', 'B'))+ background_grid(major = "xy", minor = "none")
ggsave(filename="All_vs_All.tiff",
       dpi=300,
       width=8,
       height=4,
       units="in")



#TRANSPOSONS vs STRAINS
names(summarydata)
#INSERTIONS
insertions<-summarydata[summarydata$method=="new",]
insertions<-(insertions[ order(insertions$total_tes), ])
#plot(insertions$total_tes~insertions$sample)
#pdf(file = "insertions_per_strain.pdf")
m1 <- ggplot(insertions, aes(x=reorder(insertions$sample,insertions$total_tes), y=insertions$total_tes)) 
m1<- m1 + geom_point(size=.75) +aes(group=1)+
  theme(axis.text.x = element_text(color="black",size=8,angle=90,hjust=1),
        axis.text.y = element_text(color="black",size=16,face="bold"),
        axis.title = element_text(color="black",size=16,face="bold"),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.ticks =element_line(colour = "black"))+
  labs(x="", y="Number of Insertion Sites")
m1
ggsave(filename="Insertions_per_Strain.tiff",
       dpi=300,
       width=7.5,
       height=10,
       units="in")

#ABSENCES
absences<-summarydata[summarydata$method=="absent",]
absences<-(absences[ order(absences$total_tes), ])
#plot(absences$total_tes~absences$sample)
#pdf(file = "absences_per_strain.pdf")
m2 <- ggplot(absences, aes(x=reorder(absences$sample,absences$total_tes), y=absences$total_tes)) 
m2<- m2 + geom_point(size=.75) +aes(group=1)+
  theme(axis.text.x = element_text(color="black",size=8,angle=90,hjust=1),
        axis.text.y = element_text(color="black",size=16,face="bold"),
        axis.title = element_text(color="black",size=16,face="bold"),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.ticks =element_line(colour = "black"))+
  labs(x="", y="Number of Absence Sites")
m2
ggsave(filename="Absences_per_Strain.tiff",
       dpi=300,
       width=7.5,
       height=10,
       units="in")


plot_grid(m1, m2,ncol=1)
ggsave(filename="All_per_Strain.tiff",
       dpi=300,
       width=10.5,
       height=7,
       units="in")

