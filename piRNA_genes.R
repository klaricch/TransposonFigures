#!/usr/bin/R
# this script runs the fine mappings for the trait listed in "traits_of_interest" and runs snpeff on the genes listed in "genes_of_interest"
# only exons and piRNAs are searched/maintained
# USE: piRNA_genes.R

library(cegwas)
library(devtools)
library(stringr)
#devtools::install_github("AndersenLab/cegwas")
library(dplyr)
library(tidyr)
setwd("/Users/kristen/Documents/transposon_figure_data/data")

piRNA<-scan("piRNA_gene_names.txt",what='character')

PI<-snpeff(piRNA,severity = c("HIGH","MODERATE","LOW"))

save(PI, file="PI_outliers.Rda")
load("PI_outliers.Rda")

PI_table<-filter(PI,gene_id %in% piRNA)
PI_table<-filter(PI,GT=="ALT")
PI_table<-distinct(select(PI_table,CHROM,POS,REF,ALT,strain,gene_id,feature_id))
colnames(PI_table)<-c("Chromosome","Position","Reference Allele","Alternative Allele","Strain","Gene ID", "Transcript Name")

setwd("/Users/kristen/Documents/transposon_figure_data/figures")
out_trait_info<-read.delim("Outlier_family.txt")
setwd("/Users/kristen/Documents/transposon_figure_data/data")
PI_table<-left_join(PI_table,out_trait_info,by="Strain")
write.table(PI_table, file="pi_variants.txt",sep="\t",quote=FALSE,row.names=FALSE)
#PI<-snpeff(piRNA2,variant_severity = c("HIGH","MODERATE","LOW"))
# pull "outlier" strains
#outliers<-read.table("outliers_fam_tot.txt",header=TRUE)
outliers<-read.table("outliers_families_pruned.txt",header=TRUE)


#create family and method columns
outliers<-mutate(outliers,trait=gsub("_C$","",trait))
outliers$family <- stringr::str_split_fixed(outliers$trait, "_TRANS_",2)[,2]
outliers$family<-gsub("_CE$","",outliers$family)
outliers$family<-gsub("WBTransposon","WBT",outliers$family)
outliers$method <- stringr::str_split_fixed(outliers$trait, "_TRANS_",2)[,1]
outliers<-mutate(outliers,trait=paste(family,ifelse(method=="ZERO_new"|method=="ONE_new"|method=="new","(ins)",ifelse(method=="absent","(abs)",ifelse(method=="cumulative","(cu)","(ref)"))),sep=""))
outliers<-arrange(outliers,trait)
outlier_strains<-outliers$strain
outliers_summary<-outliers %>% group_by(strain) %>% summarise(trait_outliers=toString(trait))


#snpeff("WBGene00049455")
#piRNA
#PI<-snpeff(piRNA)

#PI<-snpeff(piRNA,variant_severity = c("HIGH","MODERATE","LOW"))

PI<-distinct(PI,.keep_all=TRUE)
#PI<-filter(PI,FT=="PASS")
# calculate the number of total strains that have the REF and AlT alleles
totals_PI<-PI %>% group_by(feature_id,CHROM,POS) %>% summarise(ref_no=sum(GT=="REF"),alt_no=sum(GT=="ALT"),tot_no=sum(GT=="REF"|GT=="ALT"))
PI<-merge(totals_PI, PI,by=c("feature_id","CHROM","POS"))
# calculate the fraction of total strains that have the REF and AlT alleles
PI<-mutate(PI,frac_ref=ref_no/tot_no)
PI<-mutate(PI,frac_alt=alt_no/tot_no)
PI<-mutate(PI,rare_frac=ifelse(ref_no>=alt_no,frac_alt,frac_ref))
# find whether the REF or ALT allele is the 'minority' allele
PI<-mutate(PI,majority=ifelse(ref_no>=alt_no,"REF","ALT"))
PI<-mutate(PI,minority=ifelse(ref_no<=alt_no,"REF","ALT"))
# scorce a strain as "RARE" is it has the 'minority' allele, else scorece it as "COMMON"
PI<-mutate(PI,score=ifelse(majority!=GT,"RARE","COMMON"))
#only look at strains which are considered outlier strains for particular families or for the total traits 
PI_pruned<-filter(PI, strain %in% outlier_strains)
#pull out only variants in which the minority variant (which is usually the AlT one) , is in 10% or less of the strains (b/c we want ‘rare’ variants)
rares<-filter(PI_pruned,score=="RARE",rare_frac <=.10)
failed<-rares %>% group_by(feature_id,CHROM,POS) %>% summarise(failed_filters=sum(FT!="PASS"))
rares<-merge(rares, failed,by=c("feature_id","CHROM","POS"))
rares<-mutate(rares,fail_check=ifelse(minority=="ALT",alt_no,ref_no))
#rares<-filter(rares,failed_filters != fail_check)
rares<-merge(rares,outliers_summary, by="strain")
rares<-dplyr::select(rares,-region)
#incorporate grantham scores
#gran<-read.table("grantham_scores.txt",sep="\t")
#colnames(gran)<-c("Achange","GS")
#newRow <- data.frame(Achange="STOP",GS="STOP")
#gran <- rbind(gran,newRow) # add stops to grnatham scores
#rares<-mutate(rares,Achange=gsub("\\d+","2",aa_change))
#rares<-mutate(rares,Achange=ifelse(grepl('\\*',Achange),"STOP",gsub("p\\.","",Achange)))
#maintain aa to stop info in separate column
#rares<-mutate(rares,Achange2=gsub("\\d+","2",aa_change))
#rares<-mutate(rares,Achange2=gsub("\\*","STOP",Achange2))
#rares<-mutate(rares,Achange2=gsub("p\\.","",Achange2))

#rares<-left_join(rares, gran, by = "Achange")
save(rares, file="rares_PIs.Rda")

#test<-filter(rares, GS=="STOP")

rares2<-dplyr::select(rares, CHROM, POS,gene_id,gene_name,feature_id,nt_change,strain,trait_outliers)
rares2<-distinct(rares2,.keep_all=TRUE)
#pull out classes of grantham score severity
#stop_change<-filter(rares2, GS=="STOP") 
#raresC<-rares2
#raresC$GS<-as.numeric(raresC$GS)

#skip for now
#raresC<-filter(raresC,!is.na(GS))

#conservative<-filter(raresC, GS<=50)  %>% arrange(CHROM,POS,gene_id,strain)
#moderately_conservative<-filter(raresC, 50<GS & GS<=100) %>% arrange(CHROM,POS,gene_id,strain)
#moderately_radical<-filter(raresC, (GS>100 & GS<=150)) %>% arrange(CHROM,POS,gene_id,strain)
#radical<-filter(raresC, GS>150) %>% arrange(CHROM,POS,gene_id,strain)

#conservative$Severity="Conservative"
#moderately_conservative$Severity="Moderately Conservative"
#moderately_radical$Severity="Moderately Radical"
#radical$Severity="Radical"
#stop_change$Severity="Stop"

#severe_info<-rbind(stop_change,radical,moderately_radical,moderately_conservative,conservative)

rares2<-filter(rares2,gene_id %in% piRNA)
rares2<-arrange(rares2,CHROM,POS,gene_id,strain)
colnames(rares2)<-c("Chromosome", "Position", "Gene ID", "Gene Name", "Transcript Name", "Nucleotide Change","Strain", "Traits For Which Strain is an Outlier")
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
write.table(rares2, file="pi_outliers.txt",sep="\t",quote=FALSE,row.names=FALSE)


unique(PI$feature_id)
