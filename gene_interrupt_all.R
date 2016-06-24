#!/usr/bin/R
# this script plots a histogram showing the genomic features in which TEs are located
# USE: gene_interrupt.R

library(ggplot2)
library(grid)
library(dplyr)
library(cowplot)

setwd("/Users/kristen/Documents/transposon_figure_data/data")
data <- read.table("essentiality_nonredundant_GO.txt",sep="\t",header=TRUE,stringsAsFactors = F)

# simplify UTRs
data<-mutate(data, Region=ifelse(Region=="three_prime_UTR"|Region=="five_prime_UTR","UTR",Region))




data<-distinct(data,Chromosome, TE_start,Transcript_Name)
# simplify Biotypes
data<-mutate(data,final_bio=ifelse(Region=="intergenic","Intergenic",ifelse(Biotype=="pseudogene"|Biotype=="transposon_pseudogene","Pseudogene","Genic")))
#-split plot: A) intergenic, genic, pseudogene, B) CDS, promoter, intron
#-potential table with pseudogenes for loss of function caused by TE


region_intergenic=filter(data,final_bio=="Intergenic")
region_pseudo=filter(data,final_bio=="Pseudogene")

protein_coding<-filter(data,final_bio=="Genic", Biotype=="protein_coding")

protein_coding<-filter(protein_coding,Region!="exon")
protein_coding<-filter(protein_coding,Region!="gene")
protein_coding$Region <- factor(protein_coding$Region,
                                levels = c("promoter", "CDS","intron","UTR"),
                                labels = c("Promoter", "CDS","Intron","UTR"))
region_promoter=filter(protein_coding,Region=="Promoter")
region_CDS=filter(protein_coding,Region=="CDS")
region_intron=filter(protein_coding,Region=="Intron")
region_utr=filter(protein_coding,Region=="UTR")

pseu_info<-filter(data, Biotype=="pseudogene")
write.table(pseu_info, file="Pseudogene_Table.txt",sep="\t",quote=FALSE,row.names=FALSE)


a <- ggplot(data,aes(x=TE_start/1e6,fill=final_bio))
a <- a + geom_histogram(binwidth=.25)+
  facet_grid(.~Chromosome, scale="free", space="free_x")+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))+
  geom_point(aes(y=30), alpha=0)+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 11, colour = "black", face = "bold"),
        panel.margin = unit(.25, "lines"),
        panel.border = element_rect(fill=NA, colour="black"),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.text=element_text(size=9),
        plot.margin=unit(c(.1,.1,0,.1), "cm"),
        axis.title = element_text(size=9,face="bold"),
        axis.text.y = element_text(colour="black", size=9,face="bold"),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(colour="black", size=11,face="bold"),
        axis.ticks = element_line(colour="black"),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"))+
  labs(x="", y= "Count")+
  scale_fill_manual(values = c('Genic'="gray17",'Intergenic' = "gray60", "Pseudogene"="tan3"))
a


max_y<-ggplot_build(a)$panel$ranges[[1]]$y.range
max_y<-max_y[2]
a<- a + scale_y_continuous(expand = c(0,0),limits=c(0,max_y*1.075))
a



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
        panel.border = element_rect(fill=NA, colour="black"),
        legend.title = element_blank(),
        legend.text=element_text(size=9),
        plot.margin=unit(c(0,.1,.1,.1), "cm"),
        axis.title = element_text(size=9,face="bold"),
        axis.text.y = element_text(colour="black", size=9,face="bold"),
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

all<-plot_grid(a,b,ncol=1 )
all

setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="Genic_Features_All.tiff",
       dpi=300,
       width=7,
       height=3.5,
       units="in")




