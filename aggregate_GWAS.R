#!/usr/bin/R
# this script plots the aggregate GWAS results (ie: represents each genomic lcoation of a significant mapping with a bar)
# USE: aggregate_GWAS.R

library(ggplot2)
library(stringr)
library(dplyr)
library(grid)

load("/Users/kristen/Documents/transposon_figure_data/Processed_Transposon_Mappings.Rda")
names(final_processed_mappings)
setwd("/Users/kristen/Documents/transposon_figure_data")

names(final_processed_mappings)
load('SignificantMappings_Results_Activity.Rda')
names(Mappings)
#pull max position on each chromosome for the phantom points
max_1<-max(Mappings[(Mappings$chr=="I"), ]$pos)
max_2<-max(Mappings[(Mappings$chr=="II"), ]$pos)
max_3<-max(Mappings[(Mappings$chr=="III"), ]$pos)
max_4<-max(Mappings[(Mappings$chr=="IV"), ]$pos)
max_5<-max(Mappings[(Mappings$chr=="V"), ]$pos)
max_6<-max(Mappings[(Mappings$chr=="X"), ]$pos)

transposon <- stringr::str_split_fixed(final_processed_mappings$pheno, "_TRANS_",2)[,2]
final_processed_mappings$family <- transposon
caller <- stringr::str_split_fixed(final_processed_mappings$pheno, "_TRANS_",2)[,1]
final_processed_mappings$method <- caller

base_traits <-final_processed_mappings[(final_processed_mappings$method=="absent"| final_processed_mappings$method=="new" |final_processed_mappings$method=="reference"), ]
final_processed_mappings<-base_traits
# write out table of info on each unique peak
peaks <- final_processed_mappings[final_processed_mappings$peak_id!="NA",]
names(peaks)
#pull unique combinations of pheno and peak id
sites_clean <- distinct(peaks, peak_id, pheno)

# do we only want the distinct ones here?
names(sites_clean)
levels(sites_clean$chr)
print(sites_clean$pheno)
nrow(sites_clean)

# add te class info to summarydata(new_TRANS_end_tes will be removed)
classdata <- read.table("CtCp_all_nonredundant.txt",header=TRUE)
names(classdata)<-c("chr","start","end","TE","support","orientation","method","strain","class")
classdata$id<- stringr::str_split_fixed(classdata$TE, regex("_(non-)?reference"),2)[,1]
classdata<-mutate(classdata, pheno=paste(method,"TRANS",id,sep="_"))
class_subset <- classdata %>% distinct(pheno) %>% select(pheno,class)
print(class_subset$pheno)
sites_clean <-merge(sites_clean, class_subset, by="pheno")
nrow(sites_clean)
#revalue classes
sites_clean$class <- factor(sites_clean$class,
                            levels = c("dnatransposon", "retrotransposon","unknown"),
                            labels = c("DNA Transposon", "Retrotransposon", "Unknown"))

## add in phantom points

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


a <- ggplot(data = sites_clean, aes(x = pos/1e6, y=value)) #,colour=method
a <- a + geom_segment(aes(x = pos/1e6, y = 1, xend = pos/1e6, yend = 25,color=class))+
  facet_grid(method ~ chr,scale="free",space = "free_x",labeller=method_labeller)+
  #phantom point at x=0
 geom_point(data = sites_clean,aes(x=0, y=2),alpha=0) +
  geom_point(data = subset(sites_clean, chr=="I"),aes(x=max_1/1e6, y=2),alpha=0) +
  geom_point(data = subset(sites_clean, chr=="II"),aes(x=max_2/1e6, y=2),alpha=0) +
  geom_point(data = subset(sites_clean, chr=="III"),aes(x=max_3/1e6, y=2),alpha=0) +
  geom_point(data = subset(sites_clean, chr=="IV"),aes(x=max_4/1e6, y=2),alpha=0) +
  geom_point(data = subset(sites_clean, chr=="V"),aes(x=max_5/1e6, y=2),alpha=0) +
  geom_point(data = subset(sites_clean, chr=="X"),aes(x=max_6/1e6, y=2),alpha=0) +
  labs(x="Chromosome Position (Mb)", y="")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 9, colour = "black",face="bold"),
        panel.margin = unit(.6, "lines"),
        panel.border = element_rect(fill=NA,colour = "black"),
        panel.background = element_rect(fill = "white"),
        axis.title=element_text(size=9),
        axis.text.y = element_blank(),
        axis.ticks.x =element_line(colour = "black"),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(colour = "black",size=9),
        legend.title=element_blank(),
        legend.background = element_rect(fill=FALSE),
        legend.key=element_rect(fill=NA),
        legend.text=element_text(size=9))+
  scale_color_manual(values = c("navy", "brown3", "darkgoldenrod2"))+
  scale_y_continuous(expand = c(0,0)) 
a


setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="Aggregate_GWAS.tiff",
       dpi=300,
       width=7.5,
       height=3,
       units="in")
