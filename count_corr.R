setwd("/Users/kristen/Documents/transposon_figure_data/data")
load("Processed_Transposon_Mappings.Rda")
#keep only count traits
processed_mapping_df<-subset(processed_mapping_df,grepl('_C$', processed_mapping_df$trait))
print(unique(processed_mapping_df$trait))
# shouldn't need this step
processed_mapping_df<-filter(processed_mapping_df, !is.na(strain))
counts<- distinct(processed_mapping_df,trait,strain)
counts$family <- stringr::str_split_fixed(counts$trait, "_TRANS_",2)[,2]
counts$method <- stringr::str_split_fixed(counts$trait, "_TRANS_",2)[,1]


abs<-filter(counts, method=="absent")
refs<-filter(counts, method=="reference")

rename(refs,"ref_value"=value)
