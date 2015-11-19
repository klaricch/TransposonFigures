
setwd("/Users/kristen/Documents/transposon_figure_data/data")
load("Processed_Transposon_Mappings.Rda")

#install.packages("corrplot")
library(dplyr)
library(tidyr)
library(cegwas)

nrow(snps)
nrow(crs)
ncol(crs)
nrow(c)
ncol(c)

GWAS_LD <- function(df){
  library('corrplot') #package corrplot

  
  sn <- paste(df$chr,df$pos,sep="_")
  
  tg <- data.frame(snps)%>%
    mutate(snp = paste(CHROM,POS,sep="_"))%>%
    filter(snp %in% sn)%>%
    select(-CHROM,-POS)%>%
    gather(strain,geno,-snp)%>%
    spread(snp,geno)
  
  
  c <- cor(tg[,2:ncol(tg)])
  corrplot(c, method = "circle") #plot matrix
  
  return(c)
}

t1 <- dplyr::filter(data.frame(final_processed_mappings), pheno == unique(pheno)[50])
sns <- dplyr::filter(t1, aboveBF == 1 )
crs <- GWAS_LD(sns)

# chage the pearson corr to spearman!!!!
View(c)
sn <- paste(t1$chr,t1$pos,sep="_")

tg <- data.frame(snps)%>%
  mutate(snp = paste(CHROM,POS,sep="_"))%>%
  filter(snp %in% sn)%>%
  select(-CHROM,-POS)%>%
  gather(strain,geno,-snp)%>%
  spread(snp,geno)


c <- cor(tg[,2:ncol(tg)])
corrplot(c, method = "circle") #plot matrix

return(c)
