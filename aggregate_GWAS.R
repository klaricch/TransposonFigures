#!/usr/bin/R
# this script plots the aggregate GWAS results (ie: represents each genomic location of a significant mapping with a bar)
# USE: aggregate_GWAS.R

library(ggplot2)
library(stringr)
library(dplyr)
library(grid)

setwd("/Users/kristen/Documents/transposon_figure_data/data")
load("Processed_Transposon_Mappings.Rda")
load('SignificantMappings_Results_Activity.Rda')

#pull max position on each chromosome for the phantom points
max_1<-max(Mappings[(Mappings$chr=="I"), ]$pos)
max_2<-max(Mappings[(Mappings$chr=="II"), ]$pos)
max_3<-max(Mappings[(Mappings$chr=="III"), ]$pos)
max_4<-max(Mappings[(Mappings$chr=="IV"), ]$pos)
max_5<-max(Mappings[(Mappings$chr=="V"), ]$pos)
max_6<-max(Mappings[(Mappings$chr=="X"), ]$pos)

final_processed_mappings$pheno<- gsub("^ONE_new" ,"new",final_processed_mappings$pheno)

transposon <- stringr::str_split_fixed(final_processed_mappings$pheno, "_TRANS_",2)[,2]
final_processed_mappings$family <- transposon
caller <- stringr::str_split_fixed(final_processed_mappings$pheno, "_TRANS_",2)[,1]
final_processed_mappings$method <- caller


unique(final_processed_mappings$method)

#g<-final_processed_mappings$method
#remove fraction and movement traits
base_traits<-subset(final_processed_mappings, grepl('_C$', final_processed_mappings$pheno))
base_traits<-subset(base_traits,!grepl('^no_', base_traits$pheno))
base_traits<-subset(base_traits,!grepl('^ZERO_new', base_traits$pheno))
#base_traits <-final_processed_mappings[(final_processed_mappings$method=="absent"| final_processed_mappings$method=="new" |final_processed_mappings$method=="reference"), ]
final_processed_mappings<-base_traits
# write out table of info on each unique peak
peaks <- final_processed_mappings[final_processed_mappings$peak_id!="NA",]


#pull unique combinations of pheno and peak id
sites_clean <- distinct(peaks, peak_id, pheno)

# do we only want the distinct ones here?
names(sites_clean)
levels(sites_clean$chr)
print(sites_clean$pheno)
nrow(sites_clean)

sites_clean<-mutate(sites_clean, phenoT=gsub("_C$","",pheno))
# add te class info to summarydata(new_TRANS_end_tes will be removed)
classdata <- read.table("CtCp_all_nonredundant.txt",header=TRUE)
names(classdata)<-c("chr","start","end","TE","orientation","method","strain","class")
classdata$id<- stringr::str_split_fixed(classdata$TE, regex("_(non-)?reference"),2)[,1]
#remove chr and pos info from TE info
classdata$id <- gsub("\\w+_\\d+_" ,"",classdata$id)
print(classdata$id)


classdata<-mutate(classdata, phenoT=paste(method,"TRANS",id,sep="_"))
class_subset <- classdata %>% distinct(phenoT) %>% select(phenoT,class)
#
#
sites_clean <-merge(sites_clean, class_subset, by="phenoT")
nrow(sites_clean)
#revalue classes
sites_clean$class <- factor(sites_clean$class,
                            levels = c("dnatransposon", "retrotransposon","unknown"),
                            labels = c("DNA Transposon", "Retrotransposon", "Unknown"))



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

a<-filter(sites_clean,method=="reference")
unique(sites_clean$method)
a <- ggplot(data = sites_clean, aes(x = pos/1e6, y=value)) #,colour=method
a <- a + geom_segment(aes(x = pos/1e6, y = 1, xend = pos/1e6, yend = 25))+
  facet_grid(method ~ chr,scale="free",space = "free_x",labeller=method_labeller)+
 geom_point(data = sites_clean,aes(x=0, y=2),alpha=0) +  #phantom point at x=0
  geom_point(data = subset(sites_clean, chr=="I"),aes(x=max_1/1e6, y=2),alpha=0) +
  geom_point(data = subset(sites_clean, chr=="II"),aes(x=max_2/1e6, y=2),alpha=0) +
  geom_point(data = subset(sites_clean, chr=="III"),aes(x=max_3/1e6, y=2),alpha=0) +
  geom_point(data = subset(sites_clean, chr=="IV"),aes(x=max_4/1e6, y=2),alpha=0) +
  geom_point(data = subset(sites_clean, chr=="V"),aes(x=max_5/1e6, y=2),alpha=0) +
  geom_point(data = subset(sites_clean, chr=="X"),aes(x=max_6/1e6, y=2),alpha=0) +
  labs(x="Chromosome Position (Mb)", y="")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 9, colour = "black",face="bold"),
        panel.margin = unit(.25, "lines"),
        panel.border = element_rect(fill=NA,colour = "gray48", size=1, linetype="solid"),
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
  #scale_colour_manual(values = c('DNA Transposon' = "navy", 'Retrotransposon' = "brown3", 'Unknown' = "darkgoldenrod2"))+
  #scale_color_manual(values = c("navy", "brown3", "darkgoldenrod2"))+
  scale_y_continuous(expand = c(0,0)) 
a

setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="Aggregate_GWAS.tiff",
       dpi=300,
       width=7.5,
       height=4,
       units="in")



#piRNA
sites_clean <- mutate(sites_clean, MB=pos/(1e6))
piRNA<-filter(sites_clean, chr=="IV",MB>=4.5 & MB<=7|MB>13.5 & MB<17.2)
save(piRNA,file="piRNA_QTL.Rda")




















