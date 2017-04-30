#!/usr/bin/R
#this script generates the GWAS plots for the traits to be used as a main figure
#USE: main_GWAS_plots.R

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
#load("Processed_Transposon_Mappings_SUBSET2.Rda")
load("Processed_Transposon_Mappings_2.Rda")

unique(processed_mapping_df$trait)
load("count_QTL.Rda")

# pull unique combos, remove strain column(don't need specific strain info at this point)
processed_mapping_df<- processed_mapping_df %>% distinct(trait,marker,strain,.keep_all=TRUE)

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
unique(processed_mapping_df$trait)
#extract the count base traits
base_traits <-selection[(selection$method=="cumulative"|selection$method=="absent"| selection$method=="new" |selection$method=="reference"|selection$method=="ZERO_new"|selection$method=="ONE_new"), ]
counts<-subset(base_traits, grepl("_C$", base_traits$family))
counts$family <- gsub("_C$" ,"",counts$family)
unique(selection$trait)
processed_mapping_df <- distinct(dplyr::select(processed_mapping_df, -strain,-allele,-value),.keep_all=TRUE)
processed_mapping_df<- processed_mapping_df %>% distinct(trait,marker,.keep_all=TRUE)

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
unique(selection$trait)
#strip count marker and remnant marks from dataframes
selection$trait <- gsub("_C$" ,"",selection$trait)
processed_mapping_df$trait <- gsub("_C$" ,"",processed_mapping_df$trait)
processed_mapping_df$family <- gsub("_C$" ,"",processed_mapping_df$family)
processed_mapping_df$family <- gsub("_$" ,"",processed_mapping_df$family)
processed_mapping_df$family <- gsub("_non-reference(.*)$" ,"",processed_mapping_df$family)

processed_mapping_df<-mutate(processed_mapping_df,ID=paste(trait,peak_id,sep="_"))
copy<-processed_mapping_df
processed_mapping_df<-mutate(processed_mapping_df,SNP_col=ifelse(ID %in% count_QTL$trait,"PASS","PASS" )) #R

count_QTL<-mutate(count_QTL, trait2=gsub("_\\d+$","",trait)) 
#selection <- filter(selection, (trait %in% count_QTL$trait2))

processed_mapping_df<-filter(processed_mapping_df,CHROM != "MtDNA")
class_subset<- positions %>% distinct(class,family,.keep_all=TRUE) %>% dplyr::select(class,family)

total_plot<-filter(selection,trait=="ONE_new_TRANS_total")
selection <-merge(selection, class_subset, by="family")
selection<-arrange(selection,class,family,method)

unique(selection$trait)
selection<-filter(selection, trait=="ONE_new_TRANS_total"|trait=="ONE_new_TRANS_PAL3A"|trait=="cumulative_TRANS_MIRAGE1"|trait== "absent_TRANS_Tc1"|trait=="absent_TRANS_CER4-I") #R
unique(selection$trait)
count<-0





unique(processed_mapping_df$trait)

for (i in unique(total_plot$trait)){
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
  plot_title<-gsub("total","Total",plot_title)
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
          panel.border = element_rect(color="black", size=0.5, linetype="solid", fill=NA),
          panel.spacing = unit(.6, "lines"),
          panel.background = element_blank(),
          axis.ticks =element_line(colour = "black"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(colour = "black"),
        #  axis.title.y = element_text(colour=ifelse(class_TE=="dnatransposon","navy",ifelse(class_TE=="retrotransposon","brown3","darkgoldenrod2"))),
          axis.title=element_text(size=9),
          plot.margin=unit(c(.05,.30,-.5,.30), "cm"),
         legend.title = element_text(size = 10, angle = 270),
          legend.text = element_blank(),
          legend.key.size=unit(0,"cm"),
          legend.key = element_rect(colour = "pink"),
          legend.position=('right'))+
    labs(x="",y="",colour="red")+
    scale_color_identity()+
    scale_fill_continuous(name=plot_title)+
    scale_y_continuous(breaks= pretty_breaks(),expand=c(0,0),limits=c(0,specific_trait_mx+.075*specific_trait_mx),labels = function(x) format(x,width = 4))
  A
  
  if (count==0){B<-A+theme(strip.background = element_rect(fill = "white"),
                           strip.text.x = element_text(size = 9, colour = "black",face="bold"));
                first<-B}
  if (count==1){second<-A}
  if (count==2){third<-A}
  if (count==3){fourth<-A}
  if (count==4){fifth<-A}
  if (count==5){sixth<-A}
  count<-count+1
}


A

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
          panel.border = element_rect(color="black", size=0.5, linetype="solid", fill=NA),
          panel.spacing = unit(.6, "lines"),
          panel.background = element_blank(),
          axis.ticks =element_line(colour = "black"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(colour = "black"),
          axis.title.y = element_text(colour=ifelse(class_TE=="dnatransposon","navy",ifelse(class_TE=="retrotransposon","brown3","darkgoldenrod2"))),
          axis.title=element_text(size=11),
          plot.margin=unit(c(.05,.30,-.5,.30), "cm"),
          legend.title = element_text(size = 11, colour = ifelse(class_TE=="dnatransposon","navy",ifelse(class_TE=="retrotransposon","brown3","darkgoldenrod2")), angle = 270),
          legend.text = element_blank(),
          legend.key.size=unit(0,"cm"),
          legend.key = element_rect(colour = "pink"),
          legend.position=('right'))+
    labs(x="",y="",colour="red")+
    scale_color_identity()+
    scale_fill_continuous(name=plot_title)+
    scale_y_continuous(breaks= pretty_breaks(),expand=c(0,0),limits=c(0,specific_trait_mx+.075*specific_trait_mx),labels = function(x) format(x,width = 4))
  A
  
  if (count==0){B<-A+theme(strip.background = element_rect(fill = "white"),
    strip.text.x = element_text(size = 11, colour = "black",face="bold"));
                  first<-B}
  if (count==1){second<-A}
  if (count==2){third<-A}
  if (count==3){fourth<-A}
  if (count==4){fifth<-A}
  if (count==5){sixth<-A}
  count<-count+1
}
#R
#Tc1, Mirage1, PAL3A, and then CER4-1.
a_all<-plot_grid(first,fourth,second,third,fifth,ncol=1) #ZER)_new_TRANS_NeSL-1_C no  longer in here so don't need fifth 
label<-expression(bold(-log["10"](p)))

a_all<- a_all + draw_label(label, x = .04, y = 0.5, hjust = .5, vjust = .5,
                           fontfamily = "", fontface = "bold", colour = "black", size = 11,
                           angle = 90, lineheight = 0.9, alpha = 1)

df <- data.frame(1,2)
blank_plot<-ggplot(df,aes(x=1,y=1)) + geom_point(color="white") + theme(axis.line=element_blank(),axis.text =element_blank(),axis.ticks =element_blank(),axis.title =element_blank(),panel.background = element_blank(),panel.grid = element_blank())
a_all<-plot_grid(a_all,blank_plot,ncol=1,rel_heights = c(1, .03))


a_all<- a_all + draw_label("Chromosome Position (Mb)", x = .5, y = 0.020, hjust = .5, vjust = .5,
                           fontfamily = "", fontface = "bold", colour = "black", size = 11,
                           angle = 0, lineheight = 0.9, alpha = 1)
a_all
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="five_trait_QTL.tiff",
       dpi=350,
       width=6.75,
       height=9,
       units="in")

ggsave(filename="five_trait_QTL.png",
       dpi=350,
       width=7.5,
       height=6,
       units="in")

