#!/usr/bin/R
# this script compares the family counts for 15 samples with and without optical duplicates removed
# USE: OD_comparison.R
library(tidyr)
library(dplyr)

setwd("/Users/kristen/Documents/transposon_figure_data/data")

OD<-read.table("OD_FINAL_RESULTS_LF.txt",sep="\t")
#IG<-read.table("OD_original_FINAL_RESULTS_LF.txt",sep="\t")
IG<-read.table("FINAL_RESULTS_LF.txt_IGNORE",sep="\t")
colnames(OD)<-c("date","strain","trait","value")
colnames(IG)<-c("date","strain","trait","value")

OD<-select(OD,-date)
IG<-select(IG,-date)

counts<-merge(OD,IG, by=c("strain","trait"))
colnames(counts)<-c("strain","trait","OD","IG")

counts$family <- stringr::str_split_fixed(counts$trait, "_TRANS_",2)[,2]
counts$method <- stringr::str_split_fixed(counts$trait, "_TRANS_",2)[,1]


#initiate an empty dataframe
corr_df <- data.frame(trait_name=character(),
                      rho=numeric())

#spearman correlation
for (i in unique(counts$trait)){
  print(i)
  countsub<-counts[counts$trait==i,]
  correlation<-cor.test(counts$OD, counts$IG,method="spearman",exact=FALSE)
  rho_value <-correlation$estimate
  corr_df <- rbind(corr_df, data.frame(trait_name = i, rho=rho_value)) 
}

unique(corr_df$rho)
