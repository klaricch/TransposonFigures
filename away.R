#!/usr/bin/R
# this script checks:
#1) which QTL are located more than 100 SNPs away from the physical location of
#   the TE family that was mapped 
#2) which QTL have medians that are not equal (for both 0/1 and multi count traits)
#   and outputs such QTL to new dataframes
#3) which QTL are not in LD with a QTL within 100 SNPs of a TE position
# NOTE: processes position traits
# USE: away.R

library(dplyr)
library(ggplot2)
library(data.table)
library(grid)
library(stringr)
library(gridExtra)
library(stringr)
library(tidyr)
setwd("/Users/kristen/Documents/transposon_figure_data/data")
load("Processed_Transposon_Mappings.Rda")


# read in positions of all TEs
positions <- read.table("CtCp_all_nonredundant.txt",header=TRUE)
names(positions)<-c("CHROM","start","end","TE","orientation","method","strain","class")
#create TE family column
positions$family<- stringr::str_split_fixed(positions$TE, regex("_(non-)?reference"),2)[,1]
positions$family<- paste(stringr::str_split_fixed(positions$family, "_",4)[,3],stringr::str_split_fixed(positions$family, "_",4)[,4],sep="_")
positions$family <- gsub("_$" ,"",positions$family)
positions$family <- gsub("_non-reference(.*)$" ,"",positions$family)

#pull out only position traits from mappings dataframe
position_traits<-subset(processed_mapping_df,
                        grepl('^I', processed_mapping_df$trait) |
                          grepl('^V', processed_mapping_df$trait) |
                          grepl('^X', processed_mapping_df$trait))
#create family column
position_traits$family  <- paste(stringr::str_split_fixed(position_traits$trait, "_",4)[,3],stringr::str_split_fixed(position_traits$trait, "_",4)[,4],sep="_")
position_traits$family <- gsub("_$" ,"",position_traits$family)
position_traits$family <- gsub("_non-reference(.*)$" ,"",position_traits$family)

#pull out unique phenotypes and peak ids,CHROM,POS
distinct_sites <- distinct(position_traits,trait,CHROM,peak_id,POS)
distinct_sites<-arrange(distinct_sites, trait,peak_id,log10p) #sort by ascending trait and descending pvalue (ascending log10p)
distinct_sites <- filter(distinct_sites,!is.na(peak_id))#remove NA peaks
#distinct_sites <- filter(distinct_sites,!is.na(allele))
distinct_sites <- distinct(distinct_sites,trait,peak_id)

#pull out all unique marker used for the mappings
SNP<- distinct(processed_mapping_df,CHROM,POS)

#initiate list
one <- vector()
two <- vector()
three <- vector()
four <- vector()
five <- vector()
six <- vector()

#create separate lists for each chromosome, with the SNP positions for that chromosome as the values
for(i in 1:nrow(SNP)) {
  row <- SNP[i,]
  CHROM=row$CHROM
  POS=row$POS
  #CHROM=row[,row$CHROM]
  if (CHROM=="I"){
    one<-c(one,POS)
  }
  else if (CHROM=="II"){
    two<-c(two,POS)
  }
  else if (CHROM=="III"){
    three<-c(three,POS)
  }
  else if (CHROM=="IV"){
    four<-c(four,POS)
  }
  else if (CHROM=="V"){
    five<-c(five,POS)
  }
  else if (CHROM=="X"){
    six<-c(six,POS)
  }
  
}

#sort above lists
one<-sort(one)
two<-sort(two)
three<-sort(three)
four<-sort(four)
five<-sort(five)
six<-sort(six)

#find max length of each list
max_one=length(one)
max_two=length(two)
max_three=length(three)
max_four=length(four)
max_five=length(five)
max_six=length(six)

#initiate an empty dataframe
away<- as.data.frame(matrix(0, ncol = 18, nrow = 0))
names(away) <-c(colnames(position_traits))

not_away<- as.data.frame(matrix(0, ncol = 18, nrow = 0))
names(not_away) <-c(colnames(position_traits))

#iterate through unique QTL peaks
for (phenotype in unique(distinct_sites$trait)){
  TE_one <- vector()
  TE_two <- vector()
  TE_three <- vector()
  TE_four <- vector()
  TE_five <- vector()
  TE_six <- vector()
  transposon  <- paste(stringr::str_split_fixed(phenotype, "_",4)[,3],stringr::str_split_fixed(phenotype, "_",4)[,4],sep="_")
  transposon <- gsub("_$" ,"",transposon)
  transposon <- gsub("_non-reference(.*)$" ,"",transposon)
  locations<-filter(positions, family==transposon)
  locations<-distinct(locations,CHROM,start,end,family)
  for(i in 1:nrow(locations)) { #add all positions of a TE to a chromosome list
    row <- locations[i,]
    TE_chr=row$CHROM
    TE_pos=row$start

    if (TE_chr=="I"){
      TE_one<-c(TE_one,TE_pos)
    }
    else if (TE_chr=="II"){
      TE_two<-c(TE_two,TE_pos)
    }
    else if (TE_chr=="III"){
      TE_three<-c(TE_three,TE_pos)
    }
    else if (TE_chr=="IV"){
      TE_four<-c(TE_four,TE_pos)
    }
    else if (TE_chr=="V"){
      TE_five<-c(TE_five,TE_pos)
    }
    else if (TE_chr=="X"){
      TE_six<-c(TE_six,TE_pos)
    }
    
  }
  TE_one<-sort(TE_one)
  TE_two<-sort(TE_two)
  TE_three<-sort(TE_three)
  TE_four<-sort(TE_four)
  TE_five<-sort(TE_five)
  TE_six<-sort(TE_six)
  
  pheno_subset=distinct_sites[distinct_sites$trait==phenotype,] #pull out mappings realted to that phenotype
  #for each row in subset 
  for(i in 1:nrow(pheno_subset)) {
    trig=FALSE #start trig at FALSE, if the TE is within 100 marker of the QTL, change trig to TRUE
    row <- pheno_subset[i,]
    CHROM<-row$CHROM
    bp<-row$POS # the top SNP for the QTL peak
    
    #check which chromosome it's in
    if (CHROM=="I"){
      index <-which(one==bp)
      lower_range=index-100
      upper_range=index+100
      
      #dont go past first or last index index
      if (lower_range<1){lower_range<-1}
      if (upper_range>max_one){upper_range<-max_one}
      lower=one[lower_range]
      upper=one[upper_range]

      for(ind in TE_one){ # for given position where that TE is found
        if(ind>lower & ind<upper){
          trig<-TRUE
          not_away<-rbind(not_away,row)
          break
          }
        #if not----> rbind
      }
    }
    
    
    #TWO
    if (CHROM=="II"){
      index <-which(two==bp)
      lower_range=index-100
      upper_range=index+100
      
      #dont go past first index
      if (lower_range<1){lower_range<-1}
      if (upper_range>max_two){upper_range<-max_two}
      lower=two[lower_range]
      upper=two[upper_range]

      for(ind in TE_two){ # for given position where that TE is found
        if(ind>lower & ind<upper){
          trig<-TRUE
          not_away<-rbind(not_away,row)
          break
          }
        #if not----> rbind
      }
    } 
    
    #THREE
    if (CHROM=="III"){
      index <-which(three==bp)

      lower_range=index-100
      upper_range=index+100
      
      #dont go past first index
      if (lower_range<1){lower_range<-1}
      if (upper_range>max_three){upper_range<-max_three}
      lower=three[lower_range]
      upper=three[upper_range]

      for(ind in TE_three){ # for given position where that TE is found
        if(ind>lower & ind<upper){
          trig<-TRUE
          not_away<-rbind(not_away,row)
          break
          }
        #if not----> rbind
      }
    } 
    
    ##FOUR
    if (CHROM=="IV"){
      index <-which(four==bp)

      lower_range=index-100
      upper_range=index+100
      
      #dont go past first index
      if (lower_range<1){lower_range<-1}
      if (upper_range>max_four){upper_range<-max_four}
      lower=four[lower_range]
      upper=four[upper_range]

      for(ind in TE_four){ # for given position where that TE is found
        if(ind>lower & ind<upper){
          trig<-TRUE
          not_away<-rbind(not_away,row)
          break
          }
        #if not----> rbind
      }
    } 
    
    #FIVE
    if (CHROM=="V"){
      index <-which(five==bp)

      lower_range=index-100
      upper_range=index+100
      
      #dont go past first index
      if (lower_range<1){lower_range<-1}
      if (upper_range>max_five){upper_range<-max_five}
      lower=five[lower_range]
      upper=five[upper_range]

      for(ind in TE_five){ # for given position where that TE is found
        if(ind>lower & ind<upper){
          trig<-TRUE
          not_away<-rbind(not_away,row)
          break
          }
        #if not----> rbind
      }
    } 
    
    #SIX
    if (CHROM=="X"){
      index <-which(six==bp)

      lower_range=index-100
      upper_range=index+100
      
      #dont go past first index
      if (lower_range<1){lower_range<-1}
      if (upper_range>max_six){upper_range<-max_six}
      lower=six[lower_range]
      upper=six[upper_range]

      for(ind in TE_six){ # for given position where that TE is found
        if(ind>lower & ind<upper){
          trig<-TRUE
          not_away<-rbind(not_away,row)
          break
          }
        #if not----> rbind
      }
    } 
    if(trig==FALSE){
      away<-rbind(away,row)#bind with entire row
    }
  }
}
save(away, file="away_phenos.Rda")

away<-mutate(away,ID=paste(trait,peak_id,sep="_"))

##########################################################################################
#                               MEDIAN Differences
##########################################################################################

#get medians of positions traits that are in away df:processed_mapping_df
position_traits<-mutate(position_traits,ID=paste(trait,peak_id,sep="_"))
median_df <- position_traits[position_traits$ID %in% away$ID,]
#unique(median_df$ID)
median_df <- filter(median_df,!is.na(peak_id))
median_df <- filter(median_df,!is.na(allele))
median_df <- distinct(median_df,strain,trait,peak_id)
unique(median_df$ID)
#calculate median for each phenotype allele group
median_df <- median_df %>% group_by(ID,allele) %>% summarise(med=median(value,na.rm=TRUE))
#median_df<-mutate(median_df, Trait=paste(trait,peak_id,sep="_"))
#median_df <- median_df %>% ungroup() %>% select(-trait,-peak_id)
median_df<-spread(median_df, allele,med)
#
colnames(median_df)<-c("trait","ref","alt")
#pull out only those phenos in which the median value of the strains with the ref allele doesn't match the medidan of those with the alt allele
median_df<-filter(median_df,ref!=alt)
#median_df$trait <- gsub("_[0-9]+$" ,"",median_df$trait) # at this point onyl saving trait, acutlly want ID?
position_QTL<-median_df
save(position_QTL, file="position_QTL.Rda")
#save(median_df, file="median_phenos.Rda")

d<-filter(processed_mapping_df,strain=="N2",!is.na(allele))
##########################################################################################
#                               LD check
##########################################################################################
# library(cegwas)
# 
# GWAS_LD <- function(df){
#   library('corrplot') #package corrplot
#   sn <- paste(df$CHROM,df$POS,sep="_")
#   tg <- data.frame(snps)%>%
#     mutate(snp = paste(CHROM,POS,sep="_"))%>%
#     filter(snp %in% sn)%>%
#     select(-CHROM,-POS)%>%
#     gather(strain,geno,-snp)%>%
#     spread(snp,geno)
#   
#   c <- cor(tg[,2:ncol(tg)],method="spearman")
#   corrplot(c, method = "circle") #plot matrix 
#   return(c)
# }
# 
# #initiate an empty dataframe that will cotain the rho value between "away" and "non-away" peaks for a given phenotype
# final_LD<- as.data.frame(matrix(0, ncol = 4, nrow = 0))
# 
# for (i in unique(away$trait)) { # go though away (not median!) dataframe
#   t1 <- dplyr::filter(data.frame(distinct_sites), trait == i )
#   sns <- dplyr::filter(t1, aboveBF == 1 )
#   if (nrow(sns)>1){ #if only one QTL and already in away df, nothing to compare for LD
#     crs <- GWAS_LD(sns)
#     crs_df<-as.data.frame(crs,row.names=rownames(crs))
#     crs_df$first_QTL<-rownames(crs_df)
#     ld<-gather(crs_df,"second_QTL","LD",1:(ncol(crs_df)-1))
#     compareN<- filter(not_away, trait==i) # the QTL with TEs within 100 SNPs
#     compareA<- filter(away, trait==i) # the QTL with TEs at least 100 SNPs away
#     ld$trait<-i
#     in_ld<-filter(ld, second_QTL  %in% compareN$marker, first_QTL %in% compareA$marker, i %in% away$trait)
#     if (nrow(in_ld) >0) {  #this condition may not be met if all snps are in the "away category"
#       final_LD<-rbind(final_LD,in_ld)
#     }
#   }
# }
# 
# final_LD <- mutate(final_LD,abs_rho=abs(LD)) # take the absolute value of Spearman's rho
# final_LD <- mutate(final_LD,ID=paste(trait,first_QTL,sep="_")) #creat ID: trait+SNP
# high_ld<-filter(final_LD,abs_rho>.5)
# away_df<-away
# away_df <- mutate(away_df,ID=paste(trait,marker,sep="_"))
# low_ld<-filter(away_df, !(ID %in% high_ld$ID)) #remove QTLs that were in high with a "position QTL"
# 
# save(low_ld, file="low_ld.Rda")
# 
# #plot histogram of LD correlation values  
# m <- ggplot(final_LD,aes(x=abs_rho))
# m <- m + geom_histogram(bin=.01) +
#   theme_bw()+
#   scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
# m
# 
# setwd("/Users/kristen/Documents/transposon_figure_data/figures")
# ggsave(filename="LD_Histogram.tiff",
#        dpi=300,
#        width=7,
#        height=2.5,
#        units="in")
