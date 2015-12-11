
setwd("/Users/kristen/Documents/transposon_figure_data/data")
load("Processed_Transposon_Mappings.Rda")

#install.packages("corrplot")
library(dplyr)
library(tidyr)
library(cegwas)


##################################################################################################
#                               LD
##################################################################################################

GWAS_LD <- function(df){
  library('corrplot') #package corrplot

  
  sn <- paste(df$chr,df$pos,sep="_")
  
  tg <- data.frame(snps)%>%
    mutate(snp = paste(CHROM,POS,sep="_"))%>%
    filter(snp %in% sn)%>%
    select(-CHROM,-POS)%>%
    gather(strain,geno,-snp)%>%
    spread(snp,geno)
  
  
  c <- cor(tg[,2:ncol(tg)],method="pearson")
  test<-cor(tg[,2:ncol(tg)])
  corrplot(c, method = "circle") #plot matrix
  
  return(c)
}


# change the pearson corr to spearman!!!!


##################################################################################################
#                               PXG
##################################################################################################
library(ggplot2)
final_processed_mappings<-subset(final_processed_mappings,
                                 grepl('^I', final_processed_mappings$pheno) |
                                   grepl('^V', final_processed_mappings$pheno) |
                                   grepl('^X', final_processed_mappings$pheno)|
                                   grepl('_C$', final_processed_mappings$pheno))
hm<-final_processed_mappings
hm<-distinct(final_processed_mappings, pheno,strain)

#need to remove the activity traits


gwasPxG <- function(trt){
  #load("~/Dropbox/AndersenLab/RCode/Stefan/good_gwasMappingsINlinkage_phenotypes.Rda")
  
  hm %>%
    filter(pheno==trt)%>%
    ggplot(.)+
    aes(x=allele,y = value,fill=as.factor(allele))+
    geom_boxplot(outlier.shape=NA,size =1)+
    #scale_fill_brewer(palette = "Set2")+
    geom_jitter(size = 3, alpha = .8)  +
    theme_bw()+
    theme(axis.text.x = element_text(size=16, face="bold", color="black"),
          axis.text.y = element_text(size=14, face="bold", color="black"),
          axis.title.x = element_text(size=20, face="bold", color="black"),
          axis.title.y = element_text(size=20, face="bold", color="black",vjust=1),
          strip.text.x = element_text(size=20, face="bold", color="black"),
          strip.text.y = element_text(size=20, face="bold", color="black"),
          plot.title = element_text(size=24, face="bold", vjust=1),
          legend.title = element_text(size=14),
          panel.border = element_rect(size=1, colour = "black"),
          legend.position = "none")+
    labs( x = "Genotype")+
    scale_fill_manual( values = c("darkgoldenrod1", "indianred1"))
}

gwasPxG("ONE_new_TRANS_Tc1A_C")

