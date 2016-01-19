#this script replaces the numeric identifiers with their corresponding trait names
library(dplyr)

setwd("/Users/kristen/Documents/transposon_figure_data/data")
load("Processed_Transposon_Mappings_id.Rda")
IDS<-read.table("key_T_kin_C_matrix_full_id.txt")
colnames(IDS)<-c("id", "trait")

processed_mapping_df<-merge(pr_maps,IDS, by="trait")
processed_mapping_df<-select(processed_mapping_df, -trait)
names(processed_mapping_df)[names(processed_mapping_df)=="id"] <- "trait"
save(processed_mapping_df, file="Processed_Transposon_Mappings.Rda")
