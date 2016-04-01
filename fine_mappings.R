#!/usr/bin/R
# this script runs the fine mappings for the trait listed in "traits_of_interest" and runs snpeff on the genes listed in "genes_of_interest"
# only exons and piRNAs are searched/maintained
# USE: testC.R

library(cegwas)
library(dplyr)
library(tidyr)
setwd("/Users/kristen/Documents/transposon_figure_data/data")
load("/Users/kristen/Dropbox/AndersenLab/LabFolders/Kristen/gene_functions.rda")
load("Processed_Transposon_Mappings_2.Rda")


#t<-filter(processed_mapping_df,trait=="ONE_new_TRANS_MIRAGE1_C")
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


traits_of_interest<-c("new_TRANS_total_dnatransposon",
                      "absent_TRANS_total_retrotransposon",
                      "reference_TRANS_total_retrotransposon",
                      "absent_TRANS_total_dnatransposon",
                      "reference_TRANS_total_dnatransposon",
                      "ONE_new_TRANS_MIRAGE1",
                      "absent_TRANS_Tc1",
                      "ONE_new_TRANS_WBTransposon00000046",
                      "ONE_new_TRANS_CER8-I_CE",
                      "ONE_new_TRANS_LINE2A",
                      "reference_TRANS_CER4-I_CE")
processed_mapping_df$trait
processed_mapping_df$trait <- gsub("_C$" ,"",processed_mapping_df$trait)
#total_traits<-filter(processed_mapping_df,grepl('total',trait))
selection<-filter(processed_mapping_df,trait %in% traits_of_interest)

unique(selection$trait)
#good_traits<-filter(processed_mapping_df,trait %in% traits_of_interest) # "fantastic" and "good" traits
#selection<-rbind(total_traits,good_traits)

######################################################################################################
######################################################################################################
######################################################################################################
#piRNAs
crsPI <- cegwas::variant_correlation(selection, variant_severity = c("HIGH","MODERATE","LOW"),gene_types = "exon")
p_crsPI <- cegwas::process_correlations(crsPI)
f_crsPI<-filter(p_crsPI,gene_class_description=="21U-RNA")
#load("PI.Rda")
dPI<-distinct(f_crsPI,trait, CHROM,POS,abs_spearman_cor, effect, gene_name)

traitsPI<-distinct(f_crsPI,trait, gene_name)
write.table(traitsPI, "/Users/kristen/Documents/transposon_figure_data/figures/vc_PI.txt", sep="\t",quote=FALSE,row.names=FALSE)

#coding regions
crsEX <- cegwas::variant_correlation(selection, variant_severity = c("HIGH","MODERATE"),gene_types = "exon")
p_crsEX <- cegwas::process_correlations(crsEX)
f_crsEX<-filter(p_crsEX,transcript_biotype=="Coding")
dEX<-distinct(f_crsEX,CHROM,POS,trait, abs_spearman_cor, effect, gene_name)


total_dna_ins<-filter(dEX, trait=="new_TRANS_total_dnatransposon")
total_retro_abs<-filter(dEX, trait=="absent_TRANS_total_retrotransposon")
total_retro_ref<-filter(dEX, trait=="reference_TRANS_total_retrotransposon")
total_dna_abs<-filter(dEX, trait=="absent_TRANS_total_dnatransposon")
total_dna_ref<-filter(dEX, trait=="reference_TRANS_total_dnatransposon")

cer4<-filter(dEX, trait=="reference_TRANS_CER4-I_CE")
mirage1<-filter(dEX, trait=="ONE_new_TRANS_MIRAGE1")
tc1<-filter(dEX, trait=="absent_TRANS_Tc1")
wb046<-filter(dEX, trait=="ONE_new_TRANS_WBTransposon00000046")
cer8<-filter(dEX, trait=="ONE_new_TRANS_CER8-I_CE")
line2a<-filter(dEX, trait=="ONE_new_TRANS_LINE2A")

save(total_dna_ins, file="total_dna_ins.Rda")
save(total_retro_abs, file="total_retro_abs.Rda")

save(total_retro_ref, file="total_retro_ref.Rda")
save(total_dna_abs, file="total_dna_abs.Rda")
save(total_dna_ref, file="total_dna_ref.Rda")

save(cer4, file="cer4.Rda")
save(mirage1,file="mirage1.Rda")
save(tc1,file="tc1.Rda")
save(wb046,file="wb046.Rda")
save(cer8,file="cer8.Rda")
save(line2a,file="line2a.Rda")
save(dPI, file="dPI.Rda")

test<-filter(f_crsEX,gene_id %in% lit)

length(unique(f_crsEX$POS))
length(unique(total_dna_ins$POS))
length(unique(total_retro_abs$POS))
length(unique(cer4$POS))

#test<-filter(total_dna_ins,POS %in% total_retro_abs$POS)
#test<-filter(f_crsEX,POS  %in% total_retro_abs$POS)
#unique(test$trait)
#test2<-filter(cer4,POS  %in% total_retro_abs$POS)



save(f_crsPI, file="PI.Rda")
save(f_crsEX, file="EX.Rda")
save(dPI, file="dPI.Rda")
save(dEX, file="dEX.Rda")

# pull "outlier" strains
outliers<-read.table("outliers_fam_tot.txt",header=TRUE)
outlier_strains<-outliers$strain
outliers_summary<-outliers %>% group_by(strain) %>% summarise(hhhh=toString(trait))

# genes and transcripts of interest
GOI<-snpeff(genes_of_interest)
GOI_pruned<-filter(GOI, strain %in% outlier_strains)

test<-GOI %>% group_by(feature_id) %>% summarise(ref_no=sum(GT=="REF"))
test2<-filter(GOI, feature_id=="B0379.3a")
test3<-unique(test2)

#length(unique(GOI_pruned$strain))
#TOI<-snpeff(transcripts_of_interest)
#GOI<-rbind(GOI,TOI)
save(GOI, file="GOI.Rda")


snpeff("21UR-15703")




test<-snpeff("I:10240031..10240056")
#EF046900

#http://www.wormbase.org/tools/genome/gbrowse/c_elegans_PRJNA13758?name=I:10240031..10240056


#test<-filter(rares,ifelse(minority=="ALT",fail_check != alt_no,fail_check!=ref_no))
#test<-filter(rares,failed_filters != fail_check)
#test<-mutate(test,ID=paste(feature_id,CHROM,POS,sep="_"))
#rares<-mutate(rares,ID=paste(feature_id,CHROM,POS,sep="_"))
#check<-filter(rares,!(ID %in% test$ID))



#unique(rares$GT)
#unique(LIT$majority)
#unique(LIT_pruned$tot_num)
#test2<-filter(LIT_pruned, feature_id=="B0379.3a")
#test3<-distinct(test2)


######################################################################################################
######################################################################################################
######################################################################################################
#REMOVE ALL OF BELOW LATER


#GOI<-select(GOI,-region)
#unique(GOI$effect)
##i<-filter(GOI, strain=="LKC34", GT=="ALT")
#i<-filter(GOI,GT=="ALT")
#<-distinct(i,POS,strain)

#table(i$strain)
#table(i$gene_name)
#table(i$feature_id)

#unique(i$gene_name)
#unique(i$feature_id)

#testS<-collect(cegwas::get_db())
#unique(testS$type_of)
#t<-subset(testS,grepl("21U-RNA",testS$locus))
#t<-subset(testS,grepl("mut",testS$locus))
#unique(t$locus)     

#setwd("/Users/kristen/Documents/transposon_figure_data/data")
#load("dPI.Rda")
#h<-unique(dPI$trait)
#h
#processed_mapping_df$trait <- gsub("_C$" ,"",processed_mapping_df$trait)
#test<-filter(processed_mapping_df,trait %in% h)
#TpiRNA<-filter(test, CHROM=="IV",POS>=4500000 & POS<=7000000|POS>=13500000 & POS<=17200000|
#               startPOS>=4500000 & startPOS<=7000000|startPOS>=13500000 & startPOS<=17200000|
#              endPOS>=4500000 & endPOS<=7000000|endPOS>=13500000 & endPOS<=17200000)

#TpiRNA<-filter(TpiRNA, !is.na(peak_id))
#unique(TpiRNA$trait)
#test2<-filter(test,CHROM=="IV")
#setwd("/Users/kristen/Documents/transposon_figure_data/figures")


#genes_of_interest<-c("unc-22",
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

#transcripts_of_interest=c("C08F8.2",
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


all=c("WBGene00007444",
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
      "WBGene00011323",
      "WBGene00002041",
      "WBGene00007635",
      "WBGene00008296",
      "WBGene00006925",
      "WBGene00000207",
      "WBGene00000210",
      "WBGene00011679",
      "WBGene00012158",
      "WBGene00008133",
      "WBGene00007001",
      "WBGene00003903",
      "WBGene00017641",
      "WBGene00011061",
      "WBGene00018862",
      "WBGene00022877",
      "WBGene00010263",
      "WBGene00004323",
      "WBGene00019971",
      "WBGene00019082",
      "WBGene00014220",
      "WBGene00013256",
      "WBGene00006888",
      "WBGene00002244",
      "WBGene00001598",
      "WBGene00001599",
      "WBGene00001600",
      "WBGene00022587",
      "WBGene00021391",
      "WBGene00001601",
      "WBGene00004969",
      "WBGene00010471",
      "WBGene00010470",
      "WBGene00000412",
      "WBGene00008297",
      "WBGene00010472",
      "WBGene00010473",
      "WBGene00016116",
      "WBGene00013266",
      "WBGene00012987",
      "WBGene00006928",
      "WBGene00006927",
      "WBGene00006929",
      "WBGene00006926",
      "WBGene00006930")

