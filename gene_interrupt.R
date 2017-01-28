#!/usr/bin/R
# this script plots a histogram showing the genomic features in which TEs are located
# USE: gene_interrupt.R

library(ggplot2)
library(grid)
library(dplyr)
library(cowplot)

setwd("/Users/kristen/Documents/transposon_figure_data/data")
data <- read.table("essentiality_nonredundant_GO.txt",sep="\t",header=TRUE,stringsAsFactors = F)
data<-filter(data, Method=="new")
data$TE<-gsub("_CE$","",data$TE)
data$TE<-gsub("WBTransposon","WBT",data$TE)

# simplify UTRs
data<-mutate(data, Region=ifelse(Region=="three_prime_UTR"|Region=="five_prime_UTR","UTR",Region))




#test<-distinct(data,Chromosome, TE_start)
data<-distinct(data,Chromosome, TE_start,Transcript_Name,.keep_all=TRUE)
# simplify Biotypes
data<-mutate(data,final_bio=ifelse(Region=="intergenic","Intergenic",ifelse(Biotype=="pseudogene"|Biotype=="transposon_pseudogene","Pseudogene","Genic")))
#-split plot: A) intergenic, genic, pseudogene, B) CDS, promoter, intron
#-potential table with pseudogenes for loss of function caused by TE
head(data)
str(data)
test<-distinct(data,Chromosome, TE_start)
test


region_intergenic=filter(data,final_bio=="Intergenic")
region_pseudo=filter(data,final_bio=="Pseudogene")
genic<-filter(data,final_bio=="Genic")

protein_coding<-filter(data,final_bio=="Genic", Biotype=="protein_coding")
test<-distinct(protein_coding,Chromosome, TE_start)
test

protein_coding<-filter(protein_coding,Region!="exon")
protein_coding<-filter(protein_coding,Region!="gene")
protein_coding$Region <- factor(protein_coding$Region,
                                levels = c("promoter", "CDS","intron","UTR"),
                                labels = c("Promoter", "CDS","Intron","UTR"))
region_promoter=filter(protein_coding,Region=="Promoter")
region_CDS=filter(protein_coding,Region=="CDS")
region_intron=filter(protein_coding,Region=="Intron")
region_utr=filter(protein_coding,Region=="UTR")
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
pseu_info<-filter(data, Biotype=="pseudogene")
pseu_info<-select(pseu_info, -Method,-Biotype,-Phenotype,-final_bio)
colnames(pseu_info)<-c("Chromosome","Position","Transposon","Region","Transcript Name", "Gene Name", "GO Annotation")
write.table(pseu_info, file="Pseudogene_Table.txt",sep="\t",quote=FALSE,row.names=FALSE)


a <- ggplot(protein_coding,aes(x=TE_start/1e6,fill=Region))
a <- a + geom_histogram(binwidth=.25)+
  facet_grid(.~Chromosome, scale="free", space="free_x")+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))+
  geom_point(aes(y=26), alpha=0)+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 11, colour = "black", face = "bold"),
        panel.spacing = unit(.25, "lines"),
        panel.border = element_rect(fill=NA, colour="black",size=.5, linetype="solid"),  
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
  labs(x="Chromosome Position (Mb)", y= "Count")+
  scale_fill_manual(values = c('CDS'="orange", 'Intron' = "plum2", 'Promoter' = "cornflowerblue","UTR"="olivedrab3")) 
a

setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="Genic_Features.tiff", dpi=300, width=7, height=3.5, units="in")
ggsave(filename="Genic_Features.png", dpi=300, width=7, height=3.5, units="in")



