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
                      "absent_TRANS_total_dnatransposon",
                      "ONE_new_TRANS_total",
                      "absent_TRANS_total_unknown",
                      "absent_TRANS_total",
                      "ONE_new_TRANS_MIRAGE1",
                      "absent_TRANS_Tc1",
                      "ONE_new_TRANS_WBT00000046",
                      "ONE_new_TRANS_LINE2A",
                      "absent_TRANS_CER4-I",
                      "absent_TRANS_WBT00000600",
                      "ONE_new_TRANS_MARINER4",
                      "ONE_new_TRANS_HELITRONY1A",
                      "ONE_new_TRANS_PAL3A",
                      "cumulative_TRANS_MIRAGE1",
                      "cumulative_TRANS_Tc1",
                      "new_TRANS_total_unknown",             
                      "cumulative_TRANS_total_dnatransposon",
                      "cumulative_TRANS_total")


                                       





processed_mapping_df$trait
processed_mapping_df$trait <- gsub("_C$" ,"",processed_mapping_df$trait)
selection<-filter(processed_mapping_df,trait %in% traits_of_interest)
unique(selection$trait)
######################################################################################################
######################################################################################################
######################################################################################################
#piRNAs
crsPI <- cegwas::variant_correlation(selection, quantile_cutoff_high=.5,quantile_cutoff_low=.5,variant_severity = c("HIGH","MODERATE","LOW"),gene_types = "exon",condition_trait = F)
p_crsPI <- cegwas::process_correlations(crsPI)
f_crsPI<-filter(p_crsPI,gene_class_description=="21U-RNA")
dPI<-distinct(f_crsPI,trait, CHROM,POS,abs_spearman_cor, effect, gene_name,.keep_all=TRUE)
traitsPI<-distinct(f_crsPI,trait, gene_name,.keep_all=TRUE)
write.table(traitsPI, "/Users/kristen/Documents/transposon_figure_data/figures/vc_PI.txt", sep="\t",quote=FALSE,row.names=FALSE)
#coding regions
crsEX <- cegwas::variant_correlation(selection, quantile_cutoff_high=.5,quantile_cutoff_low=.5,variant_severity = c("HIGH","MODERATE"),gene_types = "exon",condition_trait = F)
p_crsEX <- cegwas::process_correlations(crsEX)
f_crsEX<-filter(p_crsEX,transcript_biotype=="Coding")
colnames(f_crsEX)

dEX<-distinct(f_crsEX,CHROM,POS,trait, abs_spearman_cor, effect, gene_name,gene_id,molecular_name,.keep_all=TRUE)


sort(unique(dEX$gene_name))
sort(unique(f_crsEX$trait))
unique(dEX$trait)
cegwas::variant_correlation

save(f_crsPI, file="PI.Rda")
save(f_crsEX, file="EX.Rda")
save(dPI, file="dPI.Rda")
save(dEX, file="dEX.Rda")

load("dEX.Rda")
write.table(dEX, file="/Users/kristen/Documents/transposon_figure_data/figures/dEX_Table.txt",sep="\t",quote=FALSE,row.names=FALSE)


unique(selection$trait)
unique(dEX$trait)

# results: traits of interest individually
total_dna_ins<-filter(dEX, trait=="new_TRANS_total_dnatransposon")
total_retro_abs<-filter(dEX, trait=="absent_TRANS_total_retrotransposon")
total_retro_ref<-filter(dEX, trait=="absent_TRANS_total_retrotransposon")
total_dna_abs<-filter(dEX, trait=="absent_TRANS_total_dnatransposon")
total_dna_ref<-filter(dEX, trait=="absent_TRANS_total_dnatransposon")


ONE_new_TRANS_total<-filter(dEX, trait=="ONE_new_TRANS_total")
absent_TRANS_total_unknown<-filter(dEX, trait=="absent_TRANS_total_unknown")
absent_TRANS_total<-filter(dEX, trait=="absent_TRANS_total")
new_TRANS_total_unknown<-filter(dEX, trait=="new_TRANS_total_unknown")                
cumulative_TRANS_total_dnatransposon<-filter(dEX, trait=="cumulative_TRANS_total_dnatransposon")
cumulative_TRANS_total<-filter(dEX, trait=="cumulative_TRANS_total")

#"new_TRANS_total_dnatransposon",
 #                     "absent_TRANS_total_retrotransposon",
  #                    "absent_TRANS_total_dnatransposon",
   #                   "ONE_new_TRANS_total",
    #                  "absent_TRANS_total_unknown",
     #                 "absent_TRANS_total",


cer4<-filter(dEX, trait=="absent_TRANS_CER4-I")
mirage1<-filter(dEX, trait=="ONE_new_TRANS_MIRAGE1")
tc1<-filter(dEX, trait=="absent_TRANS_Tc1")
wb046<-filter(dEX, trait=="ONE_new_TRANS_WBT00000046")
line2a<-filter(dEX, trait=="ONE_new_TRANS_LINE2A")
mar4<-filter(dEX, trait=="ONE_new_TRANS_MARINER4")
wb600<-filter(dEX, trait=="absent_TRANS_WBT00000600")
hel<-filter(dEX, trait=="ONE_new_TRANS_HELITRONY1A")
pal3<-filter(dEX, trait=="ONE_new_TRANS_PAL3A")


cuml_mirage1<-filter(dEX, trait=="cumulative_TRANS_MIRAGE1")
cuml_tc1<-filter(dEX, trait=="cumulative_TRANS_Tc1")



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
save(mar4,file="mar4.Rda")
save(wb600,file="wb600.Rda")
save(hel,file="hel.Rda")
save(pal3,file="pal3.Rda")
save(cuml_mirage1,file="cuml_mirage1.Rda")
save(cuml_tc1,file="cuml_tc1.Rda")
save(dPI, file="dPI.Rda")



save(ONE_new_TRANS_total,file="ONE_new_TRANS_total.Rda")
save(absent_TRANS_total_unknown,file="absent_TRANS_total_unknown.Rda")
save(absent_TRANS_total,file="absent_TRANS_total.Rda")
save(new_TRANS_total_unknown, file="new_TRANS_total_unknown.Rda")                
save(cumulative_TRANS_total_dnatransposon, file="cumulative_TRANS_total_dnatransposon.Rda")
save(cumulative_TRANS_total, file="cumulative_TRANS_total.Rda")


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

load("mar4.Rda")
mar4_lit<-filter(mar4,gene_id %in% lit)
load("wb600.Rda")
wb600_lit<-filter(wb600,gene_id %in% lit)
load("hel.Rda")
hel_lit<-filter(hel,gene_id %in% lit)
load("pal3.Rda")
pal3_lit<-filter(pal3,gene_id %in% lit)
load("cuml_mirage1.Rda")
cuml_mirage1_lit<-filter(cuml_mirage1,gene_id %in% lit)
load("cuml_tc1.Rda")
cuml_tc1_lit<-filter(cuml_tc1,gene_id %in% lit)





load("ONE_new_TRANS_total.Rda")
ONE_new_TRANS_total_lit<-filter(ONE_new_TRANS_total,gene_id %in% lit)
load("absent_TRANS_total_unknown.Rda")
absent_TRANS_total_unknown_lit<-filter(absent_TRANS_total_unknown,gene_id %in% lit)
load("absent_TRANS_total.Rda")
absent_TRANS_total_lit<-filter(absent_TRANS_total,gene_id %in% lit)
load("new_TRANS_total_unknown.Rda")
new_TRANS_total_unknown_lit<-filter(new_TRANS_total_unknown,gene_id %in% lit)                
load("cumulative_TRANS_total_dnatransposon.Rda")
cumulative_TRANS_total_dnatransposon_lit<-filter(cumulative_TRANS_total_dnatransposon,gene_id %in% lit)
load("cumulative_TRANS_total.Rda")
cumulative_TRANS_total_lit<-filter(cumulative_TRANS_total,gene_id %in% lit)




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


load("dEX.Rda")
dEX$method <- stringr::str_split_fixed(dEX$trait, "_TRANS_",2)[,1]
dEX$family <- stringr::str_split_fixed(dEX$trait, "_TRANS_",2)[,2]
dEX<-mutate(dEX,test=paste(family,ifelse(method=="ZERO_new"|method=="ONE_new"|method=="new","(ins)",ifelse(method=="absent","(AR)",ifelse(method=="cumulative","(cu)","(ref)"))),sep=""))
setwd("/Users/kristen/Documents/transposon_figure_data/figures")

TOI<-read.table("TOI.txt")
colnames(TOI)<-"Trait"
dEX_TOI<-filter(dEX, test %in% TOI$Trait||trait %in% traits_of_interest)
unique(dEX_TOI$trait)
dEX_TOI<-arrange(dEX_TOI,test)

dEX_TOI<-dplyr::select(dEX_TOI,test,CHROM, POS,REF,ALT,nt_change,aa_change,
                gene_name, gene_id,molecular_name,effect,strain,GT,startPOS,endPOS,log10p,
                abs_spearman_cor,concise_description,gene_class_description)

dEX_TOI$log10p<-signif(dEX_TOI$log10p,4)
dEX_TOI$abs_spearman_cor<-signif(dEX_TOI$abs_spearman_cor,4)


colnames(dEX_TOI) <-c("Trait","Chromosome","Position","Reference Allele","Alternative Allele",
                      "Nucleotide Change","Amino Acid Change","Gene Name", "Gene ID", "Transcript Name",
                      "Effect","Strain","Genotype","Left Confidence Interval","Right Confidence Interval",
                      "log10p","Absoulte Value Spearman Correlation Coefficient", "Concise Description",
                      "Gene Class Description")

write.table(dEX_TOI, file="dEX_TOI_Table.txt",sep="\t",quote=FALSE,row.names=FALSE)

