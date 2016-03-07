#this script replaces the numeric identifiers with their corresponding trait names
library(dplyr)
setwd("/Users/kristen/Documents/transposon_figure_data/data")
load("Processed_Transposon_Mappings.Rda")
processed_mapping_df<-filter(processed_mapping_df, trait=="ZERO_new_TRANS_CELETC2_C"|trait=="ZERO_new_TRANS_LINE2C_C"|trait=="ZERO_new_TRANS_Tc5B_C"|trait=="ZERO_new_TRANS_WBTransposon00000637_C"|trait=="ZERO_new_TRANS_NeSL-1_C")
save(processed_mapping_df,file="Processed_Transposon_Mappings_SUBSET2.Rda")

