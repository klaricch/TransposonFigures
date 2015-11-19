
GWAS_LD <- function(df){
  library('corrplot') #package corrplot
  
  
  load("~/Dropbox/AndersenLab/LabFolders/Stefan/GWAS/Ancillary/200strain_numeric_SNPS.Rda")
  
  sn <- paste(df$chr,df$pos,sep="_")
  
  tg <- data.frame(snp1)%>%
    mutate(snp = row.names(snp1))%>%
    filter(row.names(snp1) %in% sn)%>%
    gather(strain,geno,-snp)%>%
    spread(snp,geno)
  
  
  c <- cor(tg[,2:ncol(tg)])
  corrplot(c, method = "circle") #plot matrix
}