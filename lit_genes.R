#!/usr/bin/R
# this script runs the fine mappings for the trait listed in "traits_of_interest" and runs snpeff on the genes listed in "genes_of_interest"
# only exons and piRNAs are searched/maintained
# USE: testC.R

library(cegwas)
library(devtools)
#devtools::install_github("AndersenLab/cegwas")
library(dplyr)
library(tidyr)
setwd("/Users/kristen/Documents/transposon_figure_data/data")

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

# pull "outlier" strains
#outliers<-read.table("outliers_fam_tot.txt",header=TRUE)
outliers<-read.table("outliers_families_pruned.txt",header=TRUE)
outlier_strains<-outliers$strain
outliers_summary<-outliers %>% group_by(strain) %>% summarise(hhhh=toString(trait))

LIT<-snpeff(lit)
LIT<-distinct(LIT)
#LIT<-filter(LIT,FT=="PASS")
# calculate the number of total strains that have the REF and AlT alleles
totals_LIT<-LIT %>% group_by(feature_id,CHROM,POS) %>% summarise(ref_no=sum(GT=="REF"),alt_no=sum(GT=="ALT"),tot_no=sum(GT=="REF"|GT=="ALT"))
LIT<-merge(totals_LIT, LIT,by=c("feature_id","CHROM","POS"))
# calculate the fraction of total strains that have the REF and AlT alleles
LIT<-mutate(LIT,frac_ref=ref_no/tot_no)
LIT<-mutate(LIT,frac_alt=alt_no/tot_no)
LIT<-mutate(LIT,rare_frac=ifelse(ref_no>=alt_no,frac_alt,frac_ref))
# find whether the REF or ALT allele is the 'minority' allele
LIT<-mutate(LIT,majority=ifelse(ref_no>=alt_no,"REF","ALT"))
LIT<-mutate(LIT,minority=ifelse(ref_no<=alt_no,"REF","ALT"))
# scorce a strain as "RARE" is it has the 'minority' allele, else scorece it as "COMMON"
LIT<-mutate(LIT,score=ifelse(majority!=GT,"RARE","COMMON"))
#only look at strains which are considered outlier strains for particular families or for the total traits 
LIT_pruned<-filter(LIT, strain %in% outlier_strains)
#pull out only variants in which the minority variant (which is usually the AlT one) , is in 10% or less of the strains (b/c we want ‘rare’ variants)
rares<-filter(LIT_pruned,score=="RARE",rare_frac <=.10)
failed<-rares %>% group_by(feature_id,CHROM,POS) %>% summarise(failed_filters=sum(FT!="PASS"))
rares<-merge(rares, failed,by=c("feature_id","CHROM","POS"))
rares<-mutate(rares,fail_check=ifelse(minority=="ALT",alt_no,ref_no))
#rares<-filter(rares,failed_filters != fail_check)
rares<-merge(rares,outliers_summary, by="strain")
rares<-select(rares,-region)
save(rares, file="rares.Rda")

length(unique(rares$strain))
test<-snpeff('prg-1')
test<-filter(test,GT=="ALT")
