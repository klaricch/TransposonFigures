


library(pander)
library(dplyr)
library(ggplot2)
library(data.table)
library(grid)
library(stringr)
library(gridExtra)
library(knitr)
library(tidyr)
library(scales)
library(gtable)



setwd("/Users/kristen/Documents/transposon_figure_data/data")
#load("Processed_Transposon_Mappings.Rda")


#t<-filter(processed_mapping_df,trait=="ONE_new_TRANS_LINE2C_C")
load("20160321_complete_mapping_df.Rda")
load("20160321_processed_transposons.Rda")
#load("20160321_processed_transposons.Rda")

map_df <- list()
for(i in 1:length(mapping_df)){
  map_df[[i]] <- mapping_df[[i]][[2]]
}
map_df<- rbind_all(map_df)
library(cegwas)
#processed_mapping_df<-process_mappings(map_df,transposon_phenotypes, BF=5)





IDS1<-read.table("key_T_kin_C_matrix_full_id_reduced.txt")
colnames(IDS1)<-c("id", "trait")
IDS<-IDS1
processed_mapping_df<-merge(map_df,IDS, by="trait")
copy<-processed_mapping_df
processed_mapping_df<-select(processed_mapping_df, -trait)
names(processed_mapping_df)[names(processed_mapping_df)=="id"] <- "trait"


processed_mapping_df$trait<-gsub("$","_C",processed_mapping_df$trait)
processed_mapping_df$trait<-gsub("_C_C","_C",processed_mapping_df$trait)
processed_mapping_df$trait<-gsub("coverage_C","coverage",processed_mapping_df$trait)
processed_mapping_df<-subset(processed_mapping_df, !grepl('^no_', processed_mapping_df$trait))
#load("Processed_Transposon_Mappings_SUBSET.Rda")
#load("Processed_Transposon_Mappings_SUBSET2.Rda")
#load("position_QTL.Rda")
#load("count_QTL.Rda")
#pg<-read.table("paragraphs.txt",sep="\t",header=TRUE)


# kexpand=function(){
#   cat(knit(
#   text=knit_expand(text=
#                      "<<yfig-{{cap}}-,fig.cap='{{cap}}',results='markup',echo=FALSE,fig.height={{figheight}},out.height={{outheight}}>>=\n
#                    .q\n
#                    @"
#   )
# ))}
# .q=qplot(1:10);cap="first caption";figheight=9;outheight=90
# kexpand()






# pull unique combos, remove strain column(don't need specific strain info at this point)
#processed_mapping_df <- distinct(select(processed_mapping_df, -strain,-allele,-value))
processed_mapping_df<- processed_mapping_df %>% distinct(trait,marker)


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
selection<-processed_mapping_df

#extract the count base traits
base_traits <-selection[(selection$method=="absent"| selection$method=="new" |selection$method=="reference"|selection$method=="ZERO_new"|selection$method=="ONE_new"), ]
counts<-subset(base_traits, grepl("_C$", base_traits$family))
counts$family <- gsub("_C$" ,"",counts$family)

#
#
###REMOVE LATER!!!!!



#
#
processed_mapping_df <- distinct(select(processed_mapping_df))
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


#bind count and position traits option...choose oen of below three
#selection<-rbind(counts,position_traits)
selection<-counts
#selection<-position_traits

#

#strip count marker and remnant marks from dataframes
selection$trait <- gsub("_C$" ,"",selection$trait)
hm$trait <- gsub("_C$" ,"",hm$trait)
processed_mapping_df$trait <- gsub("_C$" ,"",processed_mapping_df$trait)
processed_mapping_df$family <- gsub("_C$" ,"",processed_mapping_df$family)
processed_mapping_df$family <- gsub("_$" ,"",processed_mapping_df$family)
processed_mapping_df$family <- gsub("_non-reference(.*)$" ,"",processed_mapping_df$family)

processed_mapping_df<-mutate(processed_mapping_df,ID=paste(trait,peak_id,sep="_"))
copy<-processed_mapping_df

##FIX HERE
processed_mapping_df<-mutate(processed_mapping_df,SNP_col= "NO")
#### NEED TO FIX HERE

count_QTL<-mutate(count_QTL, trait2=gsub("_\\d+$","",trait)) 
#position_QTL<-mutate(position_QTL, trait2=gsub("_\\d+$","",trait))  



#skip the below
#selection <- filter(selection, (trait %in% count_QTL$trait2 |trait %in% position_QTL$trait2))

processed_mapping_df<-filter(processed_mapping_df,CHROM != "MtDNA")
#load("reciprocal_removals.Rda")
#processed_mapping_df<-filter(processed_mapping_df,!(trait %in% reciprocal_removals$TE))


#selection<-filter(selection,!(trait %in% reciprocal_removals$TE))
#selection<-filter(selection,grepl('total',family))
selection<-filter(selection,trait=="ZERO_new_TRANS_CELETC2"|trait=="ZERO_new_TRANS_LINE2C"|trait=="ZERO_new_TRANS_Tc5B"|trait=="ZERO_new_TRANS_WBTransposon00000637"|trait=="ZERO_new_TRANS_NeSL-1"|trait=="ONE_new_TRANS_CELETC2"|trait=="ONE_new_TRANS_LINE2C"|trait=="ONE_new_TRANS_Tc5B"|trait=="ONE_new_TRANS_WBTransposon00000637"|trait=="ONE_new_TRANS_NeSL-1")

#HERE
class_subset<- positions %>% distinct(class,family) %>% select(class,family)
#skip the below
#selection <-merge(selection, class_subset, by="family")


selection<-arrange(selection,family,method)
label<-expression(bold(-log["10"](p)))
i="ONE_new_TRANS_total"
for (i in unique(selection$trait)){
  specific_trait<- processed_mapping_df[processed_mapping_df$trait == i, ]
  empty <-specific_trait[specific_trait$method==NA,]
  #specific_trait_mx <- max(specific_trait$log10p)
  class_TE<-unique(filter(selection,trait==i)$class)
  pvalues<-filter(specific_trait,log10p !="Inf") #
  specific_trait_mx <- max(pvalues$log10p) #
  TE<-specific_trait$family[1]
  rect_data<-filter(specific_trait,SNP_col== "NO")
  plot_method<-unique(filter(selection,trait==i)$method)
  plot_title<-i
  plot_title
  
  ##check for NAs
  #sapply(Mappings, function(x)all(is.na(x)))
  A<- processed_mapping_df %>%
    filter(trait == i)%>%
    ggplot(.)+
    aes(x=POS/1e6,y=log10p)+
    #geom_rect(data=rect_data,mapping=aes(xmin=startPOS/1e6, xmax=endPOS/1e6, ymin=0, ymax= Inf),fill="thistle1", alpha=1)+
    geom_point(aes( color='black',size=.5))+
    
    facet_grid(.~CHROM,scale="free_x",space = "free_x")+scale_color_identity()  +
    ggtitle(plot_title)+
    #geom_hline(aes(yintercept=BF),color="grey60",linetype="dashed")+
    theme(strip.background = element_rect(fill = "white"),
          strip.text.x = element_text(size = 9, colour = "black",face="bold"),
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(color="black", size=0.5, linetype="solid", fill=NA),
          panel.margin = unit(.6, "lines"),
          axis.ticks =element_line(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black"),
          axis.title=element_text(size=9,face="bold"),
          plot.margin=unit(c(.1,.1,.1,.1), "cm"),
          
          #ADDD COLOR HERE STILLLL
          # plot.title = element_text(colour=ifelse(class_TE=="dnatransposon","navy",ifelse(class_TE=="retrotransposon","brown3","darkgoldenrod2"))),
          legend.position=('none'))+
    labs(x="Chromosome Position (Mb)",y=label)+
    scale_y_continuous(expand=c(0,0),limits=c(0,specific_trait_mx+.075*specific_trait_mx),labels = function(x) format(x,width = 4),breaks= pretty_breaks())
  
  plot(A)
  
  
  ggsave(filename=i,dpi=300, width=7.5,height=3.5,units="in")
   
  

}



 