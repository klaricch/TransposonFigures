#!/usr/bin/R
# this script checks the LD between peaks of the subsetted traits of interest
#USE: LD_check.R

library(cegwas)
library(dplyr)
library(tidyr)
setwd("/Users/kristen/Documents/transposon_figure_data/data")
load("Processed_Transposon_Mappings_SUBSET2.Rda")

 GWAS_LD <- function(df){
   library('corrplot') #package corrplot   
   sn <- paste(df$CHROM,df$POS,sep="_")
   tg <- data.frame(snps)%>%
     mutate(snp = paste(CHROM,POS,sep="_"))%>%
     filter(snp %in% sn)%>%
     select(-CHROM,-POS)%>%
     gather(strain,geno,-snp)%>%
     spread(snp,geno)
   
   c <- cor(tg[,2:ncol(tg)],method="spearman")
   corrplot(c, method = "circle") #plot matrix 
   return(c)
 }


interest<-unique(processed_mapping_df$trait)
distinct_sites <- distinct(processed_mapping_df,trait,CHROM,peak_id,POS)
distinct_sites<-arrange(distinct_sites, trait,peak_id,log10p)
distinct_sites <- filter(distinct_sites,peak_id!="NA")
distinct_sites <- distinct(distinct_sites,trait,peak_id)


setwd("/Users/kristen/Documents/transposon_figure_data/figures")
for (i in interest) {
  print(i)
  t1 <- dplyr::filter(data.frame(distinct_sites), trait == i)
  sns <- dplyr::filter(t1, aboveBF == 1 )
  
  if (nrow(sns)>1){
  crs <- GWAS_LD(sns)
  } else {
    crs<-0 #single peak, no need to check for LD
  }
  
  crs <- cbind(Row.Names = rownames(crs),crs)

  
  
  name<-paste("LD", i, sep ="_")
  filename<-paste(name,"txt",sep=".")
  write.table(crs, filename, sep="\t",quote=FALSE,row.names=FALSE)
}
