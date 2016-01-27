#!/usr/bin/R
# ths script prints a table of:
#1) peaks, 
#2) their positions, 
#3) the variance explained for that QTL, 
#4) and the left and right confidence intervals

library(dplyr)
library(tidyr)

setwd("/Users/kristen/Documents/transposon_figure_data/data")
load("/Users/kristen/Documents/transposon_figure_data/data/Processed_Transposon_Mappings.Rda") #processed mappings

#remove fraction and movement traits
processed_mapping_df<-subset(processed_mapping_df,
                                 grepl('^I', processed_mapping_df$trait) |
                                   grepl('^V', processed_mapping_df$trait) |
                                   grepl('^X', processed_mapping_df$trait)|
                                   grepl('_C$', processed_mapping_df$trait))
processed_mapping_df<-subset(processed_mapping_df,!grepl('^no_', processed_mapping_df$trait))

#remove "TRANS" from transposon names
processed_mapping_df$trait <- gsub("_TRANS_" ,"_",processed_mapping_df$trait)

# write out table of info on each unique peak
peaks <- filter(processed_mapping_df, !is.na(peak_id))
peaks <- filter(peaks, !is.na(allele))
#pull unique combinations of trait and peak id
sites_clean <- distinct(peaks, peak_id, trait)
table_info<-select(sites_clean, trait,peak_id,CHROM,POS,var.exp,startPOS,endPOS)
table_info<-arrange(table_info,CHROM,POS,trait)
table_info$var.exp<-round((table_info$var.exp*100),digits=2) # convert variance explained to percentage and round
colnames(table_info) <- c("Trait", "PeakID","Chromosome","BasePosition", "VarianceExplained(%)", "LeftCI", "RightCI")
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
write.table(table_info, file="Peak_Table.txt",sep="\t",quote=FALSE,row.names=FALSE)



