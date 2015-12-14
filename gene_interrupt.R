#!/usr/bin/R
# this script plots a histogram showing the genomic features in which TEs are located
# USE: gene_interrupt.R

library(ggplot2)
library(grid)
library(dplyr)

setwd("/Users/kristen/Documents/transposon_figure_data/data")
data <- read.table("essentiality_nonredundant_GO.txt",sep="\t",header=TRUE)
data<-filter(data, Method=="new")
head(data)

method_names <- list(
  'absent'="Absence",
  'new' = "Insertion",
  'reference' = "Reference"
)

method_labeller <- function(variable,value){
  if (variable=='method') {
    return(method_names[value])
  }else {
    return(as.character(value))
  }
}

m <- ggplot(data,aes(x=TE_start/1e6,fill=Region))
m <- m + geom_histogram(bin=.25)+
  facet_grid(.~Chromosome, scale="free", space="free_x",labeller=method_labeller)+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))+
  geom_point(aes(y=45), alpha=0)+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 11, colour = "black", face = "bold"),
        panel.margin = unit(.6, "lines"),
        panel.border = element_rect(fill=NA, colour="black"),
        panel.background = element_blank(),
        axis.title = element_text(size=11,face="bold"),
        axis.text.y = element_text(colour="black", size=11,face="bold"),
        axis.text.x = element_text(colour="black", size=11,face="bold"),
        axis.ticks = element_line(colour="black"))+
        #,
       # legend.position=('none'))+
  labs(x="Chromosome Position (Mb)", y= "Count")+
  scale_fill_manual(values = c('exon' = "darkturquoise", "five_prime_UTR"="brown1","three_prime_UTR"="coral",'intron' = "plum2", 'promoter' = "steelblue4","gene"="lightgoldenrod2","intergenic"="indianred4"))
m
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="Genic_Features.tiff",
       dpi=300,
       width=7,
       height=2.5,
       units="in")

