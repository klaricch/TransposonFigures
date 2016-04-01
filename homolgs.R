#!/usr/bin/R
# this script runs the snpeff on gnes and their homologs and pulls out the corresponding information for the strain(s) with an insertion
# unto the CDS of that gene (typically has lethal RNAi pehnotype)
# USE: homologs.R

library(cegwas)
library(dplyr)
library(tidyr)
setwd("/Users/kristen/Documents/transposon_figure_data/data")



homologs<-read.table("homologs.txt")
colnames(homologs)<-c("type_of","WB","transcript","gene_name","homolog","homolog_transcript")

hom_list<-list()
for (i in homologs$WB){
  hom_list<-c(hom_list,i)
}
for (i in homologs$homolog){
  hom_list<-c(hom_list,i)
}
str(hom_list)

hom_results<-snpeff(hom_list)
hom_results<-distinct(hom_results)

ins<-read.table("essentiality_nonredundant_cds.txt", sep = "\t")
colnames(ins)<-c(  "chromosome","position","method","TE","CDS","transcript","gene_name","biotype","lethal_pheno","GO","strain_no","strains")
ins_info<-select(ins,strains,transcript)
temp1<-select(homologs,WB,transcript)
temp2<-select(homologs,homolog,homolog_transcript)
colnames(temp2)<-c("WB","transcript")
WB_trans_info<-rbind(temp1,temp2)

full_ins_info<-full_join(ins_info,WB_trans_info,by="transcript")
full_ins_info<-filter(full_ins_info,!is.na(strains))
unique(full_ins_info$transcript)

colnames(full_ins_info)<-c("strains","transcript","gene_id")
hom_results<-left_join(hom_results,full_ins_info, by="gene_id")
save(hom_results, file="hom_results.Rda")
#"strains" column is the strains with an insertion into the CDS of the corresponding gene
class(hom_results$strains)
strains_to_test<-unlist(strsplit(as.character(hom_results$strains), split=", "))
unique(hom_results$transcript)

# homolog pairs
#WBGene00002041  #	WBGene00004969
#WBGene00008296	#	WBGene00008297	
#WBGene00006925 #	WBGene00006926	
#WBGene00014193	#	WBGene00045515
#WBGene00009607	#	WBGene00009609	
#WBGene00004178	#	WBGene00004179	

hom_results_strains<-filter(hom_results, strain %in% strains_to_test)
hom_results_strains<-select(hom_results_strains,CHROM,POS,strain,REF,ALT,a1,a2,GT,FT,FILTER,query,allele,effect,impact,gene_name,gene_id,feature_type,transcript_biotype,exon_intron_rank,nt_change,aa_change,protein_position,distance_to_feature,feature_id,strains,transcript)

compare1<-filter(hom_results_strains,gene_id =="WBGene00002041"  |	gene_id =="WBGene00004969")
strains_to_test<-unlist(strsplit(as.character(compare1$strains), split=", "))
compare1<-filter(compare1, strain %in% strains_to_test)
compare1<-filter(compare1, GT=="ALT")
compare2<-filter(hom_results_strains,gene_id =="WBGene00008296"	|	gene_id =="WBGene00008297")
strains_to_test<-unlist(strsplit(as.character(compare2$strains), split=", "))
compare2<-filter(compare2, strain %in% strains_to_test)
compare2<-filter(compare2, GT=="ALT")
compare3<-filter(hom_results_strains,gene_id =="WBGene00006925" 	| gene_id =="WBGene00006926")	
strains_to_test<-unlist(strsplit(as.character(compare3$strains), split=", "))
compare3<-filter(compare3, strain %in% strains_to_test)
compare3<-filter(compare3, GT=="ALT")
compare4<-filter(hom_results_strains,gene_id =="WBGene00014193"	|	gene_id =="WBGene00045515")
strains_to_test<-unlist(strsplit(as.character(compare4$strains), split=", "))
compare4<-filter(compare4, strain %in% strains_to_test)
compare4<-filter(compare4, GT=="ALT")
compare5<-filter(hom_results_strains,gene_id =="WBGene00009607"	|	gene_id =="WBGene00009609")	
strains_to_test<-unlist(strsplit(as.character(compare5$strains), split=", "))
compare5<-filter(compare5, strain %in% strains_to_test)
compare5<-filter(compare5, GT=="ALT")
compare6<-filter(hom_results_strains,gene_id =="WBGene00004178"  |	gene_id =="WBGene00004179")	
strains_to_test<-unlist(strsplit(as.character(compare6$strains), split=", "))
compare6<-filter(compare6, strain %in% strains_to_test)
compare6<-filter(compare6, GT=="ALT")

save(compare1,compare2,compare3,compare4,compare5,compare6, file="homolog_comparisons.Rda")