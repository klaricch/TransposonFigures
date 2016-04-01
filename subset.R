#this script subsets main figure traits
library(dplyr)
setwd("/Users/kristen/Documents/transposon_figure_data/data")
load("Processed_Transposon_Mappings_2.Rda")
#processed_mapping_df<-filter(processed_mapping_df, trait=="ZERO_new_TRANS_CELETC2_C"|trait=="ZERO_new_TRANS_LINE2C_C"|trait=="ZERO_new_TRANS_Tc5B_C"|trait=="ZERO_new_TRANS_WBTransposon00000637_C"|trait=="ZERO_new_TRANS_NeSL-1_C")
processed_mapping_df<-filter(processed_mapping_df,trait=="ONE_new_TRANS_MIRAGE1_C"|
  trait=="absent_TRANS_Tc1_C"|
  trait=="ONE_new_TRANS_WBTransposon00000046_C"|
  trait=="ONE_new_TRANS_CER8-I_CE_C"|
  trait=="ONE_new_TRANS_LINE2A_C"|
  trait=="reference_TRANS_CER4-I_CE_C")
#unique(processed_mapping_df$trait)
save(processed_mapping_df,file="Processed_Transposon_Mappings_SUBSET2.Rda")

