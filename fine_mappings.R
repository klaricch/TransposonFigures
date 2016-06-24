#!/usr/bin/R
# this script runs the fine mappings for the trait listed in "traits_of_interest" and runs snpeff on the genes listed in "genes_of_interest"
# only exons and piRNAs are searched/maintained
# USE: testC.R
#library(devtools)
#devtools::install_github("AndersenLab/cegwas")
library(cegwas)
library(dplyr)
library(tidyr)
setwd("/Users/kristen/Documents/transposon_figure_data/data")
load("/Users/kristen/Dropbox/AndersenLab/LabFolders/Kristen/gene_functions.rda")
load("Processed_Transposon_Mappings_2.Rda")

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
      "WBGene00003499",#muts begin here
      "WBGene00003504",
      "WBGene00003507",
      "WBGene00011323", 
      "WBGene00004178",#piRNA reg genes begin here
      "WBGene00008995",
      "WBGene00001568",
      "WBGene00017549",
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


traits_of_interest<-c("new_TRANS_total_dnatransposon",
                      "absent_TRANS_total_retrotransposon",
                      "reference_TRANS_total_retrotransposon",
                      "absent_TRANS_total_dnatransposon",
                      "reference_TRANS_total_dnatransposon",
                      "ONE_new_TRANS_MIRAGE1",
                      "absent_TRANS_Tc1",
                      "ONE_new_TRANS_WBTransposon00000046",
                      "ONE_new_TRANS_LINE2A",
                      "reference_TRANS_CER4-I_CE",
                      "reference_TRANS_CER8-I_CE", #remove later
                      "ONE_new_TRANS_HELITRONY1A_CE", #remove later
                      "ONE_new_TRANS_MARINER4" #remove later
                      )
processed_mapping_df$trait
processed_mapping_df$trait <- gsub("_C$" ,"",processed_mapping_df$trait)
selection<-filter(processed_mapping_df,trait %in% traits_of_interest)
######################################################################################################
######################################################################################################
######################################################################################################
#piRNAs
crsPI <- cegwas::variant_correlation(selection, variant_severity = c("HIGH","MODERATE","LOW"),gene_types = "exon")
p_crsPI <- cegwas::process_correlations(crsPI)
f_crsPI<-filter(p_crsPI,gene_class_description=="21U-RNA")
dPI<-distinct(f_crsPI,trait, CHROM,POS,abs_spearman_cor, effect, gene_name)
traitsPI<-distinct(f_crsPI,trait, gene_name)
write.table(traitsPI, "/Users/kristen/Documents/transposon_figure_data/figures/vc_PI.txt", sep="\t",quote=FALSE,row.names=FALSE)
#coding regions
crsEX <- cegwas::variant_correlation(selection, variant_severity = c("HIGH","MODERATE"),gene_types = "exon")
p_crsEX <- cegwas::process_correlations(crsEX)
f_crsEX<-filter(p_crsEX,transcript_biotype=="Coding")
dEX<-distinct(f_crsEX,CHROM,POS,trait, abs_spearman_cor, effect, gene_name)

save(f_crsPI, file="PI.Rda")
save(f_crsEX, file="EX.Rda")
save(dPI, file="dPI.Rda")
save(dEX, file="dEX.Rda")


write.table(dEX, file="/Users/kristen/Documents/transposon_figure_data/figures/dEX_Table.txt",sep="\t",quote=FALSE,row.names=FALSE)


unique(selection$trait)
unique(dEX$trait)

# results: traits of interest individually
total_dna_ins<-filter(dEX, trait=="new_TRANS_total_dnatransposon")
total_retro_abs<-filter(dEX, trait=="absent_TRANS_total_retrotransposon")
total_retro_ref<-filter(dEX, trait=="reference_TRANS_total_retrotransposon")
total_dna_abs<-filter(dEX, trait=="absent_TRANS_total_dnatransposon")
total_dna_ref<-filter(dEX, trait=="reference_TRANS_total_dnatransposon")
cer4<-filter(dEX, trait=="reference_TRANS_CER4-I_CE")
mirage1<-filter(dEX, trait=="ONE_new_TRANS_MIRAGE1")
tc1<-filter(dEX, trait=="absent_TRANS_Tc1")
wb046<-filter(dEX, trait=="ONE_new_TRANS_WBTransposon00000046")
line2a<-filter(dEX, trait=="ONE_new_TRANS_LINE2A")



#remove these later
R_MARINER4<-filter(dEX, trait=="ONE_new_TRANS_MARINER4")
R_HELITRONY1A<-filter(dEX, trait=="ONE_new_TRANS_HELITRONY1A_CE")
R_CER8<-filter(dEX, trait=="reference_TRANS_CER8-I_CE")


# save traits of interest individually
save(total_dna_ins, file="total_dna_ins.Rda")
save(total_retro_abs, file="total_retro_abs.Rda")
save(total_retro_ref, file="total_retro_ref.Rda")
save(total_dna_abs, file="total_dna_abs.Rda")
save(total_dna_ref, file="total_dna_ref.Rda")
save(cer4, file="cer4.Rda")
save(mirage1,file="mirage1.Rda")
save(tc1,file="tc1.Rda")
save(wb046,file="wb046.Rda")
save(line2a,file="line2a.Rda")
save(dPI, file="dPI.Rda")

#REMOVE THESE LATER
save(R_MARINER4,file="R_MARINER4.Rda")
save(R_HELITRONY1A,file="R_HELITRONY1A4.Rda")
save(R_CER8,file="R_CER8.Rda")

#testS<-collect(cegwas::get_db())
    
load("cer4.Rda")
cer4_lit<-filter(cer4,gene_id %in% lit)
load("mirage1.Rda")
mirage1_lit<-filter(mirage1,gene_id %in% lit)
load("tc1.Rda")
tc1_lit<-filter(tc1,gene_id %in% lit)
load("wb046.Rda")
wb046_lit<-filter(wb046,gene_id %in% lit)
load("line2a.Rda")
line2a_lit<-filter(line2a,gene_id %in% lit)

load("total_dna_ins.Rda")
total_dna_ins_lit<-filter(total_dna_ins,gene_id %in% lit)
load("total_retro_abs.Rda")
total_retro_abs_lit<-filter(total_retro_abs,gene_id %in% lit)
load("total_retro_ref.Rda")
total_retro_ref_lit<-filter(total_retro_ref,gene_id %in% lit)
load("total_dna_abs.Rda")
total_dna_abs_lit<-filter(total_dna_abs,gene_id %in% lit)
load("total_dna_ref.Rda")
total_dna_ref_lit<-filter(total_dna_ref,gene_id %in% lit)


#REMOVE THES LATER
R_MARINER4_lit<-filter(R_MARINER4, gene_id %in% lit)
R_HELITRONY1A_lit<-filter(R_HELITRONY1A, gene_id %in% lit)
R_CER8_lit<-filter(R_CER8, gene_id %in% lit)

