#!/usr/bin/R
# this script runs the fine mappings for the trait listed in "traits_of_interest" and runs snpeff on the genes listed in "genes_of_interest"
# only exons and piRNAs are searched/maintained
# USE: testC.R

library(cegwas)
library(devtools)
library(stringr)
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



# pull "outlier" strains
#outliers<-read.table("outliers_fam_tot.txt",header=TRUE)
outliers<-read.table("outliers_families_pruned.txt",header=TRUE)
outliers<-filter(outliers, trait != "coverage")

length(unique(outliers$trait))
unique(outliers$trait)
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


lit
?snpeff
LIT<-snpeff(lit)
LIT<-distinct(LIT,.keep_all=TRUE)
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
rares<-dplyr::select(rares,-region)
#incorporate grantham scores
gran<-read.table("grantham_scores.txt",sep="\t")
colnames(gran)<-c("Achange","GS")
newRow <- data.frame(Achange="STOP",GS="STOP")
gran <- rbind(gran,newRow) # add stops to grnatham scores
rares<-mutate(rares,Achange=gsub("\\d+","2",aa_change))
rares<-mutate(rares,Achange=ifelse(grepl('\\*',Achange),"STOP",gsub("p\\.","",Achange)))
#maintain aa to stop info in separate column
rares<-mutate(rares,Achange2=gsub("\\d+","2",aa_change))
rares<-mutate(rares,Achange2=gsub("\\*","STOP",Achange2))
rares<-mutate(rares,Achange2=gsub("p\\.","",Achange2))

rares<-left_join(rares, gran, by = "Achange")
save(rares, file="rares.Rda")
load("rares.Rda")

test<-filter(rares, GS=="STOP")

rares2<-dplyr::select(rares, CHROM, POS,GS,Achange2,gene_id,gene_name,feature_id,strain,trait_outliers)
rares2<-distinct(rares2,.keep_all=TRUE)
length(unique(rares2$strain))
#pull out classes of grantham score severity
stop_change<-filter(rares2, GS=="STOP") 
raresC<-rares2
raresC$GS<-as.numeric(raresC$GS)

#skip for now
#raresC<-filter(raresC,!is.na(GS))

conservative<-filter(raresC, GS<=50)  %>% arrange(CHROM,POS,gene_id,strain)
moderately_conservative<-filter(raresC, 50<GS & GS<=100) %>% arrange(CHROM,POS,gene_id,strain)
moderately_radical<-filter(raresC, (GS>100 & GS<=150)) %>% arrange(CHROM,POS,gene_id,strain)
radical<-filter(raresC, GS>150) %>% arrange(CHROM,POS,gene_id,strain)

conservative$Severity="Conservative"
moderately_conservative$Severity="Moderately Conservative"
moderately_radical$Severity="Moderately Radical"
radical$Severity="Radical"
stop_change$Severity="Stop"

severe_info<-rbind(stop_change,radical,moderately_radical,moderately_conservative,conservative)
unique(severe_info$strain)
severe_info<-distinct(severe_info,CHROM,POS,gene_id,strain,.keep_all=TRUE)
colnames(severe_info)<-c("Chromosome", "Position", "Grantham Score", "Amino Acid Change", "Gene ID", "Gene Name", "Transcript Name", "Strain", "Traits For Which Strain is an Outlier", "Severity of Amino Acid Change")
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
write.table(severe_info, file="Severe_Table.txt",sep="\t",quote=FALSE,row.names=FALSE)
colnames(outliers_summary)<-c("Strain","Traits")
write.table(outliers_summary, file="Outlier_family.txt",sep="\t",quote=FALSE,row.names=FALSE)
#aa<-snpeff("WBGene00003501")

