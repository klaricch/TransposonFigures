#!/usr/bin/R
# this script runs the fine mappings for the trait listed in "traits_of_interest" and runs snpeff on the genes listed in "genes_of_interest"
# only exons and piRNAs are searched/maintained
# USE: testC.R

library(cegwas)
library(dplyr)
library(tidyr)
setwd("/Users/kristen/Documents/transposon_figure_data/data")

load("/Users/kristen/Dropbox/AndersenLab/LabFolders/Kristen/gene_functions.rda")
load("Processed_Transposon_Mappings.Rda")

#################
#############
########
#####
#
lit=c("WBGene00007444", # TE control genes begin here
      "WBGene00000206",
      "WBGene00009664",
      "WBGene00000041",
      "WBGene00000209",
      "WBGene00010809",
      "WBGene00020757",
      "WBGene00004015",
      "WBGene00020383",
      "WBGene00007000",
      "WBGene00004452",
      "WBGene00008920",
      "WBGene00004450",
      "WBGene00018891",
      "WBGene00019168",
      "WBGene00010627",
      "WBGene00020915",
      "WBGene00014120",
      "WBGene00003508",
      "WBGene00007788",
      "WBGene00007789",
      "WBGene00017075",
      "WBGene00003803",
      "WBGene00003790",
      "WBGene00003902",
      "WBGene00004094",
      "WBGene00022310",
      "WBGene00003499",
      "WBGene00003504",
      "WBGene00003507",
      "WBGene00011323", #piRNA reg genes begin below here
      "WBGene00004178",
      "WBGene00008995",
      "WBGene00001568",
      "WBGene00017549",
      "WBGene00001568",
      "WBGene00012167",
      "WBGene00011061",
      "WBGene00007624",
      "WBGene00011333",
      "WBGene00019862",
      "WBGene00001601",
      "WBGene00001598",
      "WBGene00010280",
      "WBGene00017641",
      "WBGene00006888",
      "WBGene00001214",
      "WBGene00004510",
      "WBGene00004179")

processed_mapping_df$trait <- gsub("_C$" ,"",processed_mapping_df$trait)
total_traits<-filter(processed_mapping_df,grepl('total',trait))
selection<-filter(total_traits,trait=="new_TRANS_total_dnatransposon")


######################################################################################################
######################################################################################################
######################################################################################################

#piRNAs
crsPI <- cegwas::variant_correlation(selection, variant_severity = c("HIGH","MODERATE","LOW"),gene_types = "exon")
p_crsPI <- cegwas::process_correlations(crsPI)
f_crsPI<-filter(p_crsPI,gene_class_description=="21U-RNA")
dPI<-distinct(f_crsPI, abs_spearman_cor, effect, gene_name)
#coding regions
crsEX <- cegwas::variant_correlation(selection, variant_severity = c("HIGH","MODERATE"),gene_types = "exon")
p_crsEX <- cegwas::process_correlations(crsEX)
f_crsEX<-filter(p_crsEX,transcript_biotype=="Coding")
dEX<-distinct(f_crsEX, abs_spearman_cor, effect, gene_name)

#test<-filter(f_crsEX,gene_id %in% lit)


