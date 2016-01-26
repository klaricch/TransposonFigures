#this script replaces the numeric identifiers with their corresponding trait names
library(dplyr)
setwd("/Users/kristen/Documents/transposon_figure_data/data")
load("20160119_processed_mapping_df.Rda")
unique(processed_mapping_df$BF)
IDS1<-read.table("key_T_kin_C_matrix_full_id.txt")
IDS2<-read.table("key_kin_matrix_full_id.txt")
IDS3<-read.table("key_T_Full_Results_Activity_id.txt")
colnames(IDS1)<-c("id", "trait")
colnames(IDS2)<-c("id", "trait")
colnames(IDS3)<-c("id", "trait")

IDS<-rbind(IDS1,IDS2,IDS3)

processed_mapping_df<-merge(processed_mapping_df,IDS, by="trait")
copy<-processed_mapping_df
processed_mapping_df<-select(processed_mapping_df, -trait)
names(processed_mapping_df)[names(processed_mapping_df)=="id"] <- "trait"
save(processed_mapping_df, file="Processed_Transposon_Mappings_with activity.Rda")

#remove activty traits
#pull out only activity traits
processed_mapping_df<-subset(processed_mapping_df, !grepl('^no_', processed_mapping_df$trait))
save(processed_mapping_df, file="Processed_Transposon_Mappings.Rda")

# check that CER1 traits map correctly
#names(processed_mapping_df)
test<-filter(processed_mapping_df,!is.na(peak_id))
test<-distinct(test,peakPOS,trait)
AA<-filter(test,trait=="absent_TRANS_CER1_C")
RR<-filter(test,trait=="reference_TRANS_CER1_C")
