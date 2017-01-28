#!/usr/bin/R
# ths script prints a table of:
#1) peaks, 
#2) their positions, 
#3) the variance explained for that QTL, 
#4) and the left and right confidence intervals

library(dplyr)
library(tidyr)

setwd("/Users/kristen/Documents/transposon_figure_data/data")
load("/Users/kristen/Documents/transposon_figure_data/data/Processed_Transposon_Mappings_2.Rda") #processed mappings
load("count_QTL.Rda")
setwd("/Users/kristen/Documents/transposon_figure_data/figures")

TOI<-read.table("TOI.txt")
colnames(TOI)<-"Trait"


#remove fraction, movement, and position traits
processed_mapping_df<-subset(processed_mapping_df,grepl('_C$', processed_mapping_df$trait))
processed_mapping_df<-subset(processed_mapping_df,!grepl('^no_', processed_mapping_df$trait))

#QTL filter here:
processed_mapping_df$trait<-gsub("_C$","",processed_mapping_df$trait)
count_QTL<-mutate(count_QTL, trait2=gsub("_\\d+$","",trait)) 

#write out trait names to file for GWAS2.Rmd script
tn<-unique(processed_mapping_df$trait)
write.table(tn, file="trait_names.txt",sep="\t",quote=FALSE,row.names=FALSE)

#processed_mapping_df<-filter(processed_mapping_df,(trait %in% count_QTL$trait2))
#no longer doing this step
#processed_mapping_df<-filter(processed_mapping_df,(trait %in% count_QTL$trait2)|grepl("total",trait)) #R

#remove "TRANS" from transposon names
processed_mapping_df$method <- stringr::str_split_fixed(processed_mapping_df$trait, "_TRANS_",2)[,1]
processed_mapping_df$family <- stringr::str_split_fixed(processed_mapping_df$trait, "_TRANS_",2)[,2]
processed_mapping_df$trait <- gsub("_TRANS_" ,"_",processed_mapping_df$trait)

# write out table of info on each unique peak
peaks <- filter(processed_mapping_df, !is.na(peak_id))
peaks <- filter(peaks, !is.na(allele))
#pull unique combinations of trait and peak id
sites_clean <- distinct(peaks, peak_id, trait,.keep_all=TRUE)
sites_clean<-filter(sites_clean,method!="reference")
table_info<-dplyr::select(sites_clean, trait,peak_id,CHROM,POS,var.exp,log10p,startPOS,endPOS,method,family)

table_info$family<-gsub("_C$","",table_info$family)



table_info$family<-gsub("_CE$","",table_info$family)
table_info$family<-gsub("WBTransposon","WBT",table_info$family)
table_info<-mutate(table_info,test=paste(family,ifelse(method=="ZERO_new"|method=="ONE_new"|method=="new","(ins)",ifelse(method=="absent","(AR)",ifelse(method=="cumulative", "(cu)","(ref)"))),sep=""))
unique(table_info$method)
table_info$trait<-table_info$test
table_info<-dplyr::select(table_info,-method,-family,-test)


table_info<-arrange(table_info,CHROM,POS,trait)
table_info$log10p<-round(table_info$log10p,digits=2)
table_info$var.exp<-round((table_info$var.exp*100),digits=2) # convert variance explained to percentage and round
colnames(table_info) <- c("Trait", "PeakID","Chromosome","BasePosition", "VarianceExplained(%)","-log10p", "LeftCI", "RightCI")

table_info<-mutate(table_info,po=ifelse(Trait %in% TOI$Trait, "first","second"))
table_info<-arrange(table_info,po,Trait,Chromosome,BasePosition)
table_info<-dplyr::select(table_info,-po)
test<-filter(table_info,Trait %in% TOI$Trait)
#CHECK THAT  THE  BAD  PEAKS A RE  NOT I N  HERE!
min(test[,5]) #bad
max(test[,5]) #good
order_var=test[,5]
order_var<-sort(order_var)
order_var # take second minimum from here

setwd("/Users/kristen/Documents/transposon_figure_data/figures")
write.table(table_info, file="Peak_Table.txt",sep="\t",quote=FALSE,row.names=FALSE)

#length(unique(table_info$Trait))
#length(unique(table_info$trait))

#testF<-filter(table_info, !grepl("total",Trait))
#length(unique(testF$Trait))
#testT<-filter(table_info, grepl("total",Trait))
#length(unique(testT$Trait))


#length(unique(count_QTL$trait2))
