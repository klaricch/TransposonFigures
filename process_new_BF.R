library(cegwas)
library(dplyr)
library(tidyr)


setwd("/Users/kristen/Documents/transposon_figure_data/data")
load("20160718_complete_mapping_df.Rda")
load("20160718_processed_transposons.Rda")
#load("20160321_complete_mapping_df.Rda")

length(unique(mapping_df$trait))
#map_df <- list()
#for(i in 1:length(mapping_df)){
#  map_df[[i]] <- mapping_df[[i]][[2]]
#}
#map_df<- rbind_all(map_df)

#get_old_bf<-process_mappings(mapping_df,transposon_phenotypes)
#?process_mappings
#CHOICE  BELOW
#processed_mapping_df<-process_mappings(mapping_df,transposon_phenotypes, BF=5,snp_grouping=300)
processed_mapping_df<-process_mappings(mapping_df,transposon_phenotypes,snp_grouping=300)


IDS1<-read.table("key_T_kin_C_matrix_full_id_reduced.txt")
#IDS2<-read.table("key_kin_matrix_full_id.txt")
#IDS3<-read.table("key_T_Full_Results_Activity_id.txt")
colnames(IDS1)<-c("id", "trait")
#colnames(IDS2)<-c("id", "trait")
#colnames(IDS3)<-c("id", "trait")
IDS<-IDS1
#IDS<-rbind(IDS1,IDS2,IDS3)
paper_num<-merge(mapping_df,IDS,by="trait")
length(unique(paper_num$id))
a<-filter(paper_num, grepl("cumulative",id),!grepl("total",id))
unique(a$id)
length(unique(a$id))
b<-filter(paper_num, grepl("total",id))
unique(b$id)
length(unique(b$id))
c<-filter(paper_num, !grepl("total",id),!grepl("cumulative",id))
unique(c$id)
length(unique(c$id))
processed_mapping_df<-merge(processed_mapping_df,IDS, by="trait")


copy<-processed_mapping_df
processed_mapping_df<-dplyr::select(processed_mapping_df, -trait)
names(processed_mapping_df)[names(processed_mapping_df)=="id"] <- "trait"
#save(processed_mapping_df, file="Processed_Transposon_Mappings_with activity.Rda")
processed_mapping_df$trait<-gsub("$","_C",processed_mapping_df$trait)
processed_mapping_df$trait<-gsub("_C_C","_C",processed_mapping_df$trait)
processed_mapping_df$trait<-gsub("coverage_C","coverage",processed_mapping_df$trait)
processed_mapping_df<-subset(processed_mapping_df, !grepl('^no_', processed_mapping_df$trait))
processed_mapping_df<-filter(processed_mapping_df,!grepl("ZERO_new",trait))
save(processed_mapping_df, file="Processed_Transposon_Mappings_2.Rda")

#table raw mapping results
raw_mapping_table<-dplyr::select(processed_mapping_df,marker,trait,log10p)
raw_mapping_table$method <- stringr::str_split_fixed(raw_mapping_table$trait, "_TRANS_",2)[,1]
raw_mapping_table<-filter(raw_mapping_table, method!="reference")
raw_mapping_table$family <- stringr::str_split_fixed(raw_mapping_table$trait, "_TRANS_",2)[,2]
raw_mapping_table$trait <- gsub("_TRANS_" ,"_",raw_mapping_table$trait)
raw_mapping_table$family<-gsub("_C$","",raw_mapping_table$family)
raw_mapping_table<-mutate(raw_mapping_table,test=paste(family,ifelse(method=="ZERO_new"|method=="ONE_new"|method=="new","(ins)",ifelse(method=="absent","(abs)","(ref)")),sep=""))
raw_mapping_table$trait<-raw_mapping_table$test
raw_mapping_table<-dplyr::select(raw_mapping_table,-method,-family,-test)
raw_mapping_table$log10p<-signif(raw_mapping_table$log10p,4)
colnames(raw_mapping_table)<-c("Marker","Trait","Log10p")
#length(unique(raw_mapping_table$Trait))

setwd("/Users/kristen/Documents/transposon_figure_data/figures")
write.table(raw_mapping_table, file="Raw_Mapping_Table.txt",sep="\t",quote=FALSE,row.names=FALSE)

#teest snp groupings
#test<-filter(processed_mapping_df,trait=="ONE_new_TRANS_WBT00000046_C",CHROM=="IV")
#unique(test$peak_id)
#write.table(test, file="test.txt",sep="\t",quote=FALSE,row.names=FALSE)
#test2<-filter(test,peak_id=="2")

#t<-filter(processed_mapping_df,trait=="ONE_new_TRANS_MIRAGE1_C")
#unique(processed_mapping_df$BF)
