library(cegwas)
library(dplyr)
library(tidyr)
setwd("/Users/kristen/Documents/transposon_figure_data/data")
load("Processed_Transposon_Mappings_2.Rda")
 
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
distinct_sites <- distinct(processed_mapping_df,trait,CHROM,peak_id,POS)
distinct_sites<-arrange(distinct_sites, trait,peak_id,log10p)
distinct_sites <- filter(distinct_sites,peak_id!="NA")
distinct_sites <- distinct(distinct_sites,trait,peak_id)
t1 <- dplyr::filter(data.frame(distinct_sites), trait == "ONE_new_TRANS_LINE2A_C")
sns <- dplyr::filter(t1, aboveBF == 1 )
crs <- GWAS_LD(sns)



