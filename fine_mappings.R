#!/usr/bin/R
# this script runs the fine mappings for the trait listed in "traits_of_interest" and runs snpeff on the genes listed in "genes_of_interest"
# only exons and piRNAs are searched/maintained
# USE: testC.R

library(cegwas)
library(dplyr)
library(tidyr)
setwd("/Users/kristen/Documents/transposon_figure_data/data")
#load("Processed_Transposon_Mappings_SUBSET2.Rda")
load("/Users/kristen/Dropbox/AndersenLab/LabFolders/Kristen/gene_functions.rda")
load("Processed_Transposon_Mappings.Rda")

traits_of_interest<-c("ZERO_new_TRANS_CELETC2",
                      "ZERO_new_TRANS_LINE2C",
                      "ZERO_new_TRANS_Tc5B",
                      "ZERO_new_TRANS_WBTransposon00000637",
                      "reference_TRANS_CER8-I_CE",
                      "ZERO_new_TRANS_CER8-LTR_CE",
                      "ZERO_new_TRANS_LINE2G_CE",
                      "ZERO_new_TRANS_NDNAX3_CE",
                      "ZERO_new_TRANS_NeSL-1",
                      "ZERO_new_TRANS_CELE46B",
                      "ZERO_new_TRANS_IR1_CE",
                      "ZERO_new_TRANS_LINE2D",
                      "ZERO_new_TRANS_MIRAGE1",
                      "ZERO_new_TRANS_NPALTA1_CE",
                      "absent_TRANS_Tc1",
                      "ZERO_new_TRANS_Tc3",
                      "ZERO_new_TRANS_Tc5A",
                      "ZERO_new_TRANS_TIR9TA1A_CE",
                      'ZERO_new_TRANS_WBTransposon00000207',
                      "reference_TRANS_WBTransposon00000657",
                      "reference_TRANS_WBTransposon00000661")

genes_of_interest<-c("unc-22",
                     "mut-7",
                     "mut-2",
                     "mut-16",
                     "mut-7",
                     "mut-14",
                     "mut-15",
                     "cdr-2",
                     "ceh-34",
                     "egl-4",
                     "gon-2",
                     "kin-4",
                     "old-1",
                     "pkd-2",
                     "prg-1",
                     "unc-130",
                     "vit-1",
                     "asb-1",
                     "asg-1",
                     "rpl-38",
                     "rpl-36",
                     "npp-4",
                     "pab-1",
                     "ppw-2")

transcripts_of_interest=c("C08F8.2",
                          "F35G12.10",
                          "F43G9.1",
                          "F54H12.1",
                          "K07A12.3",
                          "M01F1.3",
                          "T24C4.1",
                          "T24H7.1",
                          "T09B4.9",
                          "Y71H2AM.23",
                          "C06B8.8",
                          "F17C11.9",
                          "F37C12.4",
                          "F55F8.3",
                          "H06I04.3",
                          "K07C5.4",
                          "W01B11.3",
                          "ZK858.7",
                          "B0379.3",
                          "C28A5.1b",
                          "C28A5.2b",
                          "D2096.8",
                          "F10G8.3",
                          "Y54E5A.4",
                          "Y106G6H.2",
                          "Y110A7A.18",
                          "Y77E11A.7")

processed_mapping_df$trait <- gsub("_C$" ,"",processed_mapping_df$trait)
total_traits<-filter(processed_mapping_df,grepl('total',trait))
good_traits<-filter(processed_mapping_df,trait %in% traits_of_interest) # "fantastic" and "good" traits
selection<-rbind(total_traits,good_traits)

######################################################################################################
######################################################################################################
######################################################################################################

#piRNAs
crsPI <- cegwas::variant_correlation(selection, variant_severity = c("HIGH","MODERATE","LOW"),gene_types = "exon")
p_crsPI <- cegwas::process_correlations(crsPI)
f_crsPI<-filter(p_crsPI,gene_class_description=="21U-RNA")
dPI<-distinct(f_crsPI,trait, abs_spearman_cor, effect, gene_name)

traitsPI<-distinct(f_crsPI,trait, gene_name)
write.table(traitsPI, "/Users/kristen/Documents/transposon_figure_data/figures/vc_PI.txt", sep="\t",quote=FALSE)

#coding regions
crsEX <- cegwas::variant_correlation(selection, variant_severity = c("HIGH","MODERATE"),gene_types = "exon")
p_crsEX <- cegwas::process_correlations(crsEX)
f_crsEX<-filter(p_crsEX,transcript_biotype=="Coding")
dEX<-distinct(f_crsEX,trait, abs_spearman_cor, effect, gene_name)

save(f_crsPI, file="PI.Rda")
save(f_crsEX, file="EX.Rda")
save(dPI, file="dPI.Rda")
save(dEX, file="dEX.Rda")

# genes and transcripts of interest
GOI<-snpeff(genes_of_interest)
TOI<-snpeff(transcripts_of_interest)
GOI<-rbind(GOI,TOI)
save(GOI, file="GOI.Rda")


setwd("/Users/kristen/Documents/transposon_figure_data/data")
load("dPI.Rda")
h<-unique(dPI$trait)
h
processed_mapping_df$trait <- gsub("_C$" ,"",processed_mapping_df$trait)
test<-filter(processed_mapping_df,trait %in% h)
TpiRNA<-filter(test, CHROM=="IV",POS>=4500000 & POS<=7000000|POS>=13500000 & POS<=17200000|
                startPOS>=4500000 & startPOS<=7000000|startPOS>=13500000 & startPOS<=17200000|
                endPOS>=4500000 & endPOS<=7000000|endPOS>=13500000 & endPOS<=17200000)

TpiRNA<-filter(TpiRNA, !is.na(peak_id))
unique(TpiRNA$trait)
test2<-filter(test,CHROM=="IV")
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
load("piRNA_QTL.Rda")
sort(unique(piRNA$trait))
sort(unique(dPI$trait))
######################################################################################################
######################################################################################################
######################################################################################################
#REMOVE ALL OF BELOW LATER


#ttt<-filter(tt,grepl('total',trait))
#ttt<-filter(t,grepl('total',trait))
#sort(unique(tt$gene_name))

GOI<-select(GOI,-region)
unique(GOI$effect)
#i<-filter(GOI, strain=="LKC34", GT=="ALT")
i<-filter(GOI,GT=="ALT")
i<-distinct(i,POS,strain)

table(i$strain)
table(i$gene_name)
table(i$feature_id)

unique(i$gene_name)
unique(i$feature_id)

testS<-collect(cegwas::get_db())
t<-subset(testS,grepl("mut",testS$locus))
unique(t$locus)     




#check igv and snps for 10 insertions in exons in question snapped (SNP eff)

######################################################
######################################################
######################################################
#run snpeff
#filter for piRNAs
#run through variant correlation
#testS<-collect(cegwas::get_db())
#piRNA<-collect(cegwas::get_db() %>% dplyr::filter(biotype=="piRNA"))
#test<-collect(cegwas::get_db() %>% dplyr::filter(grepl("mut", gene_name)))
#filter(testS, grepl("mut",gene_name))

#names(testS)
#head(testS)

#crs_P <- cegwas::variant_correlation(selection, variant_severity = c("HIGH","MODERATE","LOW"),gene_types = c("piRNA"))
#crs_E <- cegwas::variant_correlation(selection, variant_severity = c("HIGH","MODERATE"),gene_types = c("exon"))
#p_crs_P <- cegwas::process_correlations(crs_P)
#p_crs_E <- cegwas::process_correlations(crs_E)
