setwd("/Users/kristen/Documents/transposon_figure_data/data")
load("Processed_Transposon_Mappings.Rda")
head(processed_mapping_df)
cn_test<-filter(processed_mapping_df, !is.na(peak_id))
test<-distinct(cn_test,CHROM,trait)
#cn_test<-filter(processed_mapping_df,trait=="ONE_new_TRANS_CELE1_C"|trait=="ONE_new_TRANS_MARINER3_C"|trait=="reference_TRANS_LINE2E_C"|trait=="absent_TRANS_Tc5_C")
cn_test<-filter(processed_mapping_df,trait=="ZERO_new_TRANS_CELETC2_C"|trait=="absent_TRANS_CEMUDR1_C"|trait=="ZERO_new_TRANS_HELITRONY1A_CE_C"|trait=="absent_TRANS_MARINCE1_C")
cn_test<-distinct(cn_test,strain,trait)
cn_test<-select(cn_test,strain,value,trait)
cn_test<-spread(cn_test, trait,value)
cn_test<-filter(cn_test, !is.na(strain))

write.table(cn_test, file="cn_test.txt",sep="\t",quote=FALSE,row.names=FALSE)


100000000/(.25^11)
