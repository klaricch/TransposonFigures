library(cegwas)
library(dplyr)
library(ggplot2)
setwd("/Users/kristen/Documents/transposon_figure_data/data")
load("count_QTL.Rda")
load("Processed_Transposon_Mappings.Rda")
unique(processed_mapping_df$trait)
counts<-subset(processed_mapping_df, grepl("_C$", processed_mapping_df$trait))
unique(counts$trait)
unique(counts$value)

table(counts$value)
df<-read.table("cer_95_set.txt",header=TRUE)

AA<-filter(df, trait=="absent_TRANS_CER1_C")
partial_set <- cegwas_map(AA)
manplot(partial_set)

RR<-filter(df, trait=="reference_TRANS_CER1_C")
partial_set2 <- cegwas_map(RR)
manplot(partial_set2)

###############
###############
###############
setwd ("~/Dropbox/AndersenLab/RCode/GWAS/Data/GenomicTraits/20160115_transposons/")

load("20160115_processed_mapping_df.Rda")
load("201601015_processed_transposons.Rda")
#load("Processed_Transposon_Mappings.Rda")
names(processed_mapping_df)
library(dplyr)
CC<-filter(processed_mapping_df, trait=="absent_TRANS_CER1_C")
DD<-filter(processed_mapping_df, trait=="reference_TRANS_CER1_C")

graphCC<- CC %>%
  ggplot(.)+
  aes(x=POS/1e6,y=log10p)+
  geom_point(aes( color=ifelse(log10p> BF, 'red', 'black')),size=1)+
  facet_grid(.~CHROM,scale="free_x",space = "free_x")+scale_color_identity()+
  geom_hline(aes(yintercept=BF),color="grey60",linetype="dashed")+
  theme(strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color="black", size=0.5, linetype="solid", fill=NA),
        panel.margin = unit(.6, "lines"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.title=element_text(size=9),
        plot.margin=unit(c(.1,.1,-.25,.1), "cm"),
        legend.position=('none'))+
  labs(x="",y="-log10(p)")
graphCC


graphDD<- DD %>%
  ggplot(.)+
  aes(x=POS/1e6,y=log10p)+
  geom_point(aes( color=ifelse(log10p> BF, 'red', 'black')),size=1)+
  facet_grid(.~CHROM,scale="free_x",space = "free_x")+scale_color_identity()+
  geom_hline(aes(yintercept=BF),color="grey60",linetype="dashed")+
  theme(strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 9, colour = "black",face="bold"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color="black", size=0.5, linetype="solid", fill=NA),
        panel.margin = unit(.6, "lines"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.title=element_text(size=9),
        plot.margin=unit(c(.1,.1,-.25,.1), "cm"),
        legend.position=('none'))+
  labs(x="",y="-log10(p)")
graphDD


load("away_counts_phenos.Rda")
library(dplyr)
singles<-read.table("singles.txt",header=TRUE)
t<-filter(processed_mapping_df,trait %in% singles$family)
final_t<-unique(t$trait)
print(final_t)
t<-distinct(t,trait,peakPOS)
save(t,file="pa.Rda")


load("pa.Rda")

bb<-distinct(processed_mapping_df,trait,peakPOS)
aa<-filter(bb,peakPOS=="13458478")
aa<-filter(bb,peakPOS=="5942259")


cc<-filter(processed_mapping_df,trait=="reference_TRANS_CER10.I_CE_C", !is.na(peak_id))
