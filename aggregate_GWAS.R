#!/usr/bin/R
# this script plots the aggregate GWAS results (ie: represents each genomic location of a significant mapping with a bar)
# USE: aggregate_GWAS.R

library(ggplot2)
library(stringr)
library(dplyr)
library(grid)

setwd("/Users/kristen/Documents/transposon_figure_data/data")
load("Processed_Transposon_Mappings.Rda")

#pull max position on each chromosome for the phantom points
max_1<-max(processed_mapping_df[(processed_mapping_df$CHROM=="I"), ]$POS)
max_2<-max(processed_mapping_df[(processed_mapping_df$CHROM=="II"), ]$POS)
max_3<-max(processed_mapping_df[(processed_mapping_df$CHROM=="III"), ]$POS)
max_4<-max(processed_mapping_df[(processed_mapping_df$CHROM=="IV"), ]$POS)
max_5<-max(processed_mapping_df[(processed_mapping_df$CHROM=="V"), ]$POS)
max_6<-max(processed_mapping_df[(processed_mapping_df$CHROM=="X"), ]$POS)

##MAYBE NEED TO GET RID OF THIS
processed_mapping_df$trait<- gsub("^ONE_new" ,"new",processed_mapping_df$trait)

# creat family and caller columns
transposon <- stringr::str_split_fixed(processed_mapping_df$trait, "_TRANS_",2)[,2]
processed_mapping_df$family <- transposon
caller <- stringr::str_split_fixed(processed_mapping_df$trait, "_TRANS_",2)[,1]
processed_mapping_df$method <- caller


#remove fraction and movement traits
base_traits<-subset(processed_mapping_df, grepl('_C$', processed_mapping_df$trait))
base_traits<-subset(base_traits,!grepl('^no_', base_traits$trait))
base_traits<-subset(base_traits,!grepl('^ZERO_new', base_traits$trait))
#base_traits <-processed_mapping_df[(processed_mapping_df$method=="absent"| processed_mapping_df$method=="new" |processed_mapping_df$method=="reference"), ]
processed_mapping_df<-base_traits
# write out table of info on each unique peak
peaks<-filter(processed_mapping_df,!is.na(peak_id))

#pull unique combinations of trait and peak id
sites_clean <- distinct(peaks, peak_id, trait)
all_sites_clean<-sites_clean # this df will be used to retain the total counts

# do we only want the distinct ones here?
names(sites_clean)
levels(sites_clean$CHROM)
print(sites_clean$trait)
nrow(sites_clean)

sites_clean<-mutate(sites_clean, traitT=gsub("_C$","",trait))
# add te class info to summarydata(new_TRANS_end_tes will be removed)
classdata <- read.table("CtCp_all_nonredundant.txt",header=TRUE)
names(classdata)<-c("CHROM","start","end","TE","orientation","method","strain","class")
classdata$id<- stringr::str_split_fixed(classdata$TE, regex("_(non-)?reference"),2)[,1]
#remove CHROM and pos info from TE info
classdata$id <- gsub("\\w+_\\d+_" ,"",classdata$id)

classdata<-mutate(classdata, traitT=paste(method,"TRANS",id,sep="_"))
class_subset <- classdata %>% distinct(traitT) %>% select(traitT,class)

#
sites_clean$traitTset<-gsub("_C$","",sites_clean$traitT)
sites_clean <-merge(sites_clean, class_subset, by="traitT")
#after the merge there are less rows because no longer including the "total" counts, only the family ones

#revalue classes
sites_clean$class <- factor(sites_clean$class,
                            levels = c("dnatransposon", "retrotransposon","unknown"),
                            labels = c("DNA Transposon", "Retrotransposon", "Unknown"))
#revalue methods
sites_clean$method <- factor(sites_clean$method,
                             levels = c("new","reference","absent"),
                             labels = c("Insertion", "Reference","Absence"))

a <- ggplot(data = sites_clean, aes(x = POS/1e6, y=value)) #,colour=method
a <- a + geom_segment(aes(x = POS/1e6, y = 1, xend = POS/1e6, yend = 25,color=class))+
  facet_grid(method ~ CHROM,scale="free",space = "free_x")+
 geom_point(data = sites_clean,aes(x=0, y=2),alpha=0) +  #phantom point at x=0
  geom_point(data = subset(sites_clean, CHROM=="I"),aes(x=max_1/1e6, y=2),alpha=0) +
  geom_point(data = subset(sites_clean, CHROM=="II"),aes(x=max_2/1e6, y=2),alpha=0) +
  geom_point(data = subset(sites_clean, CHROM=="III"),aes(x=max_3/1e6, y=2),alpha=0) +
  geom_point(data = subset(sites_clean, CHROM=="IV"),aes(x=max_4/1e6, y=2),alpha=0) +
  geom_point(data = subset(sites_clean, CHROM=="V"),aes(x=max_5/1e6, y=2),alpha=0) +
  geom_point(data = subset(sites_clean, CHROM=="X"),aes(x=max_6/1e6, y=2),alpha=0) +
  labs(x="Chromosome Position (Mb)", y="")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 9, colour = "black",face="bold"),
        panel.margin = unit(.50, "lines"),
        panel.border = element_rect(fill=NA,colour = "black", size=1, linetype="solid"),
        panel.background = element_rect(fill = "white"),
        axis.title=element_text(size=9,face="bold"),
        axis.text.y = element_blank(),
        axis.ticks.x =element_line(colour = "black"),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(colour = "black",size=9),
        legend.title=element_blank(),
        legend.background = element_rect(fill=FALSE),
        legend.key=element_rect(fill=NA),
        legend.position="none",
        legend.text=element_text(size=9))+
  scale_colour_manual(values = c('DNA Transposon' = "navy", 'Retrotransposon' = "brown3", 'Unknown' = "darkgoldenrod2"))+
  scale_y_continuous(expand = c(0,0)) 
a

setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="Aggregate_GWAS.tiff",
       dpi=300,
       width=7.5,
       height=4,
       units="in")

#piRNA
piRNA<-filter(all_sites_clean, CHROM=="IV",POS>=4500000 & POS<=7000000|POS>=13500000 & POS<=17200000)
save(piRNA,file="piRNA_QTL.Rda")


















