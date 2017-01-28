#!/usr/bin/R
# this script checks:
#1) which QTL are located more than 100 SNPs away from the physical location of
#   the TE family that was mapped 
#2) which QTL have medians that are not equal (for both 0/1 and multi count traits)
#   and outputs QTL that do not meet the above requirements to a new dataframe
# can uncomment lines to check for QTL in LD and QTL found on more than 2 chromosomes
# NOTE: processes raw_count traits
# USE: away_counts.R

library(dplyr)
library(ggplot2)
library(data.table)
library(grid)
library(stringr)
library(gridExtra)
library(stringr)
library(tidyr)
setwd("/Users/kristen/Documents/transposon_figure_data/data")
load("Processed_Transposon_Mappings_2.Rda")
#unique(test$trait)
##test<-filter(processed_mapping_df,trait=="ONE_new_TRANS_MIRAGE1_C")
##test<-filter(processed_mapping_df,trait=="ONE_new_TRANS_Tc1_C")
##test<-filter(processed_mapping_df,trait=="ONE_new_TRANS_Tc6_C")
#remove fraction and movement traits
processed_mapping_df<-subset(processed_mapping_df,
                                 grepl('^I', processed_mapping_df$trait) |
                                   grepl('^V', processed_mapping_df$trait) |
                                   grepl('^X', processed_mapping_df$trait)|
                                   grepl('_C$', processed_mapping_df$trait))
processed_mapping_df<-subset(processed_mapping_df,!grepl('^no_', processed_mapping_df$trait))

#pull distinct combos of trait, marker, and strain
processed_mapping_df<- processed_mapping_df %>% distinct(trait,marker,strain,.keep_all=TRUE)

#create family and method columns
processed_mapping_df$family <- stringr::str_split_fixed(processed_mapping_df$trait, "_TRANS_",2)[,2]
processed_mapping_df$method <- stringr::str_split_fixed(processed_mapping_df$trait, "_TRANS_",2)[,1]
#select traits above BF.....this step not needed, double checking everything is above BF
selection<-filter(processed_mapping_df, log10p > BF)
unique(selection$method)
#extract the count base traits
base_traits <-selection[(selection$method=="cumulative"|selection$method=="absent" |selection$method=="reference"|selection$method=="ZERO_new"|selection$method=="ONE_new"), ]
counts<-subset(base_traits, grepl("_C$", base_traits$family))
counts$family <- gsub("_C$" ,"",counts$family)

#length(unique(counts$trait))

#count_df<-filter(counts,!is.na(peak_id))
#count_df<-filter(count_df,!is.na(allele))
#count_df <- distinct(count_df,trait,strain)
##calculate median, min, max for each phenotype allele group
#count_df_temp1 <- count_df %>% group_by(trait) %>% summarise(minimum=min(value,na.rm=TRUE),maximum=max(value,na.rm=TRUE))
#count_df_temp2 <- count_df %>% group_by(trait,allele) %>% summarise(med=median(value,na.rm=TRUE)) #dont actually have to calculate median here b/c do this step later
library(tidyr)
#count_df_temp2<-spread(count_df_temp2, allele,med)
##merge min/max and median dataframes
#count_df<-merge(count_df_temp1, count_df_temp2,by="trait")
#colnames(count_df)<-c("trait","minimum", "maximum","ref","alt")

#count_df_pa<-filter(count_df,minimum=="0",maximum=="1")
test<-filter(processed_mapping_df, !grepl("reference",trait))
length(unique(test$trait))
#get count pa traits where the medians are the same between the ref and alt allele, then remove these from the count_df
#count_df_pa<-filter(count_df_pa,ref==alt)
#count_df<-filter(count_df,!(trait %in% count_df_pa$trait))
all_counts<-counts
#counts<-all_counts

unique(counts$trait)
#get the one occurence count traits from the full dataset
#counts<-filter(counts,(trait %in% count_df_pa$trait))
counts$trait <- gsub("_C$" ,"",counts$trait)


# read in positions of all TEs
positions <- read.table("CtCp_all_nonredundant_reduced.txt")
names(positions)<-c("CHROM","start","end","TE","orientation","method","strain","class")
#create TE family column

positions$family<- stringr::str_split_fixed(positions$TE, regex("_(non-)?reference"),2)[,1]
positions$family<- paste(stringr::str_split_fixed(positions$family, "_",4)[,3],stringr::str_split_fixed(positions$family, "_",4)[,4],sep="_")
positions$family <- gsub("_$" ,"",positions$family)
positions$family <- gsub("_non-reference(.*)$" ,"",positions$family)
positions<-mutate(positions, traitT=paste(method,"TRANS",family,sep="_"))
#class_subset <- positions %>% distinct(traitT) %>% select(traitT,class)
class_subset <- positions %>% distinct(family,class,.keep_all=TRUE) %>% dplyr::select(family,class)

counts<-mutate(counts,method2=gsub("ONE_","",method))
counts<-mutate(counts,method3=gsub("ZERO_","",method2))
counts<-mutate(counts, traitT=paste(method3,"TRANS",family,sep="_"))
#unique(counts$trait)
#unique(counts$family)
#counts<-merge(counts, class_subset, by="family") #only traits not merging here will be the "total" traits (they will drop out)
counts<-dplyr::select(counts, -method2,-method3,-traitT)


# get rid of reciprocal traits according to class
#counts<-filter(counts,ifelse(class=="retrotransposon", method=="absent"| method=="ONE_new",ifelse(class=="dnatransposon",method=="absent"| method=="ONE_new",method=="absent"| method=="ONE_new")))
counts<-filter(counts, method=="absent"| method=="ONE_new"|method=="cumulative")
unique(counts$trait)
# counts<- filter(counts,ifelse(class=="retrotransposon", method=="absent"| method=="ONE_new",ifelse(class=="dnatransposon",method=="absent"| method=="ONE_new",method=="absent"| method=="ONE_new")))
#counts<- filter(counts,ifelse(class=="retrotransposon",method=="ZERO_new"| method=="reference"| method=="ONE_new",ifelse(class=="dnatransposon",method=="ZERO_new"|method=="absent"| method=="ONE_new",method=="ZERO_new"|method=="absent"|method=="reference"|| method=="ONE_new")))
#test1<- filter(counts,class=="retrotransposon",method=="ZERO_new"| method=="reference")
#test2<- filter(counts,class=="dnatransposon",method=="ZERO_new"| method=="absent")
#test3<- filter(counts,class=="unknown",method=="ZERO_new"| method=="reference"| method=="absent")

#length(unique(counts$trait))
#length(unique(x$trait))

#z<-filter(counts,!(counts$trait %in% x$trait))
#zz<-distinct(z,trait,method)
#xx<-distinct(x,trait,method)
#cc<-distinct(counts,trait,method)
#zzz<-filter(x,family=="CEMUDR1")

#nrow(zz)
#nrow(xx)
#nrow(cc)
#nrow(test1)+nrow(test2)+nrow(test3)
#qtest1<-distinct(test1,method,class,trait)
#qtest2<-distinct(test2,method,class,trait)
#qtest3<-distinct(test3,method,class,trait)


#pull out only position traits from mappings dataframe
#position_traits<-subset(processed_mapping_df,
 #                       grepl('^I', processed_mapping_df$trait) |
 #                         grepl('^V', processed_mapping_df$trait) |
  #                        grepl('^X', processed_mapping_df$trait))
#create family column


#position_traits$family  <- paste(stringr::str_split_fixed(position_traits$trait, "_",4)[,3],stringr::str_split_fixed(position_traits$trait, "_",4)[,4],sep="_")
#position_traits$family <- gsub("_$" ,"",position_traits$family)
#position_traits$family <- gsub("_non-reference(.*)$" ,"",position_traits$family)

#pull out unique phenotypes and peak ids,chr,pos
distinct_sites <- distinct(counts,trait,CHROM,peak_id,POS,.keep_all=TRUE)
distinct_sites<-arrange(distinct_sites, trait,peak_id,log10p)
distinct_sites <- filter(distinct_sites,peak_id!="NA")
distinct_sites <- distinct(distinct_sites,trait,peak_id,.keep_all=TRUE)

unique(distinct_sites$trait)
#pull out all unique SNPs used for the mappings
SNP<- distinct(processed_mapping_df,CHROM,POS,.keep_all=TRUE)

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
away<- as.data.frame(matrix(0, ncol = 20, nrow = 0))
names(away) <-c(colnames(counts))

not_away<- as.data.frame(matrix(0, ncol = 20, nrow = 0))
names(not_away) <-c(colnames(counts))


#remove totals
distinct_sites<-filter(distinct_sites,!grepl('total', family))
distinct_sites$family
#iterate through unique QTL peaks
for (phenotype in unique(distinct_sites$trait)){
  TE_one <- vector()
  TE_two <- vector()
  TE_three <- vector()
  TE_four <- vector()
  TE_five <- vector()
  TE_six <- vector()
  transposon  <- stringr::str_split_fixed(phenotype, "_TRANS_",2)[,2]
  transposon <- gsub("_C$" ,"",transposon)
  print(phenotype)
  print(transposon)
  #unique(positions$family)
  locations<-filter(positions, family==transposon)
  # locations<-filter(positions, family=="CElE14B")
  # print(locations$family)
  locations<-distinct(locations,CHROM,start,end,family,.keep_all=TRUE)
  
  
  for(i in 1:nrow(locations)) { #add all positions of a TE to a chromosome list
    row <- locations[i,]
    TE_CHROM=row$CHROM
    TE_pos=row$start
    
    if (TE_CHROM=="I"){
      TE_one<-c(TE_one,TE_pos)
    }
    else if (TE_CHROM=="II"){
      TE_two<-c(TE_two,TE_pos)
    }
    else if (TE_CHROM=="III"){
      TE_three<-c(TE_three,TE_pos)
    }
    else if (TE_CHROM=="IV"){
      TE_four<-c(TE_four,TE_pos)
    }
    else if (TE_CHROM=="V"){
      TE_five<-c(TE_five,TE_pos)
    }
    else if (TE_CHROM=="X"){
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
 # pheno_subset=filter(distinct_sites,trait=="ONE_new_TRANS_CELE14B") #pull out mappings realted to that phenotype
#  lb
#  ub
#  CHROM
#  pheno_subset
#  i=2
#  row
#  i=1
  for(i in 1:nrow(pheno_subset)) {
    trig=FALSE #start trig at FALSE, if the TE is within 100 SNPs of the QTL, change trig to TRUE
    row <- pheno_subset[i,]
    CHROM<-row$CHROM
    bp<-row$POS # the top SNP for the QTL peak
    lb<-row$startPOS #left CI
    ub<-row$endPOS #right CI
    #check which chromosome it's in
    if (CHROM=="I"){
      #index <-which(one==bp)
      lb_index <-which(one==lb)
      ub_index <-which(one==ub)
      lower_range=lb_index-100
      upper_range=ub_index+100
      
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
      lb_index <-which(two==lb)
      ub_index <-which(two==ub)
      lower_range=lb_index-100
      upper_range=ub_index+100
      
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
      #index <-which(three==bp)
      
      lb_index <-which(three==lb)
      ub_index <-which(three==ub)
      lower_range=lb_index-100
      upper_range=ub_index+100
      
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
      # index <-which(four==bp)
      
      lb_index <-which(four==lb)
      ub_index <-which(four==ub)
      lower_range=lb_index-100
      upper_range=ub_index+100
      
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
      # index <-which(five==bp)
      
      lb_index <-which(five==lb)
      ub_index <-which(five==ub)
      lower_range<-lb_index-100
      upper_range<-ub_index+100
      
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
    # lb
    #ub
    #ind
    #SIX
    if (CHROM=="X"){
      #index <-which(six==bp)
      
      lb_index <-which(six==lb)
      ub_index <-which(six==ub)
      lower_range=lb_index-100
      upper_range=ub_index+100
      
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
table_info<-mutate(not_away,test=paste(family,ifelse(method=="ZERO_new"|method=="ONE_new"|method=="new","(ins)",ifelse(method=="absent","(AR)",ifelse(method=="cumulative","(cu)","(ref)"))),sep=""))
setwd("/Users/kristen/Documents/transposon_figure_data/tables")
table_info<-mutate(table_info,close="near")
table_info<-dplyr::select(table_info,test,peak_id,close)
write.table(table_info, file="close_table.txt",sep="\t",quote=FALSE,row.names=FALSE)

away_counts<-away
save(away_counts, file="away_counts_phenos.Rda")
away_counts<-mutate(away_counts,ID=paste(trait,peak_id,sep="_"))


unique(away_counts$trait)
##########################################################################################
#                               MEDIAN Differences
##########################################################################################

#counts <- filter(counts,!is.na(peak_id))
#counts <- filter(counts,!is.na(allele))
#get medians of positions traits that are in away df:processed_mapping_df

#d<-filter(counts,!is.na(peak_id))
#c<-filter(counts,ID=="absent_TRANS_Tc9_2")
#e<-filter(away_counts,ID=="absent_TRANS_Tc9_2")
counts<-mutate(counts,ID=paste(trait,peak_id,sep="_"))



###THIS STEP CHOICE
#median_df<-filter(counts,counts$ID %in% away_counts$ID)
median_df<-counts
#median_df <- counts[counts$ID %in% away_counts$ID]

median_df <- filter(median_df,!is.na(peak_id))
median_df <- filter(median_df,!is.na(allele))
median_df <- distinct(median_df,trait,strain,peak_id,.keep_all=TRUE)

unique(median_df$allele)
#calculate median for each phenotype allele group
median_df <- median_df %>% group_by(ID,allele) %>% summarise(med=median(value,na.rm=TRUE))
median_df<-ungroup(median_df)
#median_df<-mutate(median_df, Trait=paste(trait,peak_id,sep="_"))
#median_df <- median_df %>% ungroup() %>% select(-trait,-peak_id)
median_df<-spread(median_df, allele,med)
colnames(median_df)<-c("trait","ref","alt")
#pull out only those phenos in which the median value of the strains with the ref allele doesn't match the medidan of those with the alt allele

#CHOICE HERE
#median_df<-dplyr::filter(median_df,ref!=alt)
median_match<-dplyr::filter(median_df,ref==alt)
#median_df$trait <- gsub("_[0-9]+$" ,"",median_df$trait)
median_match<-filter(median_match, trait %in% counts$ID)
median_match$family <- stringr::str_split_fixed(median_match$trait, "_TRANS_",2)[,2]
median_match$method <- stringr::str_split_fixed(median_match$trait, "_TRANS_",2)[,1]
median_match$peak_id <- stringr::str_split_fixed(median_match$trait, ".*_",2)[,-1]
median_match$name <- stringr::str_split_fixed(median_match$family, "_\\d+$",2)[,1]
#select traits above BF.....this step not needed, double checking everything is above BF
table_info<-mutate(median_match,test=paste(name,ifelse(method=="ZERO_new"|method=="ONE_new"|method=="new","(ins)",ifelse(method=="absent","(AR)","(ref)")),sep=""))
table_info<-mutate(table_info,median="same")
table_info<-dplyr::select(table_info,test,peak_id,median)
write.table(table_info, file="median_table.txt",sep="\t",quote=FALSE,row.names=FALSE)




#length(unique(away$trait))
#length(unique(median_df$trait))
count_QTL<-median_df
save(count_QTL, file="count_QTL.Rda")
unique(count_QTL$trait)
#median_df_counts<-median_df
#save(median_df_counts, file="median_phenos_counts.Rda")
#all_counts$trait <- gsub("_C$" ,"",all_counts$trait)
##########################################################################################
#                               FEW CHRO
##########################################################################################
# 
# #get medians of positions traits that are in away df:
# chr_df <- all_counts[all_counts$peak_id !="NA",]
# chr_df <- distinct(chr_df,CHROM,trait)
# 
# chr_df<- chr_df%>% group_by(trait) %>% summarise(no_chr=length(unique(CHROM)))
# few_chr<-filter(chr_df,no_chr<=2)
# high_chr<-filter(chr_df,no_chr>2)

test<-filter(positions,family=="CELE14B")

##########################################################################################
#                               MEDIAN Differences ALLL 
##########################################################################################

#get medians of positions traits that are in away df:
#Amedian_df <- all_counts[all_counts$trait %in% away$trait,]

#Amedian_df  <- filter(all_counts,!is.na(peak_id))
#Amedian_df  <- filter(Amedian_df ,!is.na(allele))
#Amedian_df <- distinct(Amedian_df,trait,strain,peak_id)

#calculate median for each phenotype allele group
#Amedian_df <- Amedian_df %>% group_by(trait,allele,peak_id) %>% summarise(med=median(value,na.rm=TRUE))

#Amedian_df<-mutate(Amedian_df, Trait=paste(trait,peak_id,sep="_"))
#Amedian_df <- Amedian_df %>% ungroup() %>% select(-trait,-peak_id)
#Amedian_df<-spread(Amedian_df, allele,med)
#colnames(Amedian_df)<-c("trait","ref","alt")

#pull out only those phenos in which the median value of the strains with the ref allele doesn't match the medidan of those with the alt allele
#Amedian_df<-mutate(Amedian_df,median_diff=abs(ref-alt))
#Amedian_df<-mutate(Amedian_df,diff=ifelse(ref==alt,"SAME","DIFFERENT"))
#ABmedian_df<-Amedian_df
#Amedian_df$trait <- gsub("_C_[0-9]+$" ,"",Amedian_df$trait) # get rid to the "_C"
#ABmedian_df$trait <- gsub("_[0-9]+$" ,"",ABmedian_df$trait) # keep the "_C"
#Amedian_df$trait <- gsub("_C$" ,"",Amedian_df$trait)
#save(Amedian_df, file="Amedian.Rda")

#different<-filter(ABmedian_df, diff == "DIFFERENT")

#unique(different$peak_id)
##########################################################################################
#                               Merge
##########################################################################################

#good_traits<-filter(counts,  (trait %in% away_counts$trait)) #100 SNPs away filter
#good_traits<-filter(good_traits,  (trait %in% median_df_counts$trait))  #median filter

#good_traits<-filter(good_traits,  (trait %in% low_ld_counts$trait)) 

#bad_traits<-filter(counts,  !(trait %in% good_traits$trait)) 
#bad_traits2<-filter(all_counts,  !(trait %in% few_chr$trait)) #chr filter
#bad_traits3<-filter(all_counts,  !(trait %in% different$trait)) #chr filter
#bad_traits2$trait <- gsub("_C$" ,"",bad_traits2$trait)
#bad_traits3$trait <- gsub("_C$" ,"",bad_traits3$trait)
#length(unique(bad_traits2$trait))

#merged <- rbind(bad_traits, bad_traits2)
#merged <- rbind(bad_traits, bad_traits3)

#all_counts$trait <- gsub("_C$" ,"",all_counts$trait)
#counts_to_remove<-filter(all_counts, (trait %in% bad_traits$trait))
#counts_to_remove<-filter(all_counts, (trait %in% merged$trait))
#length(unique(counts_to_remove$trait))
#save(counts_to_remove, file="counts_to_remove.Rda")



##########################################################################################
# #                               LD check
# ##########################################################################################
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
# low_ld_counts<-low_ld
# save(low_ld_counts, file="low_ld_counts.Rda")
# 
# #plot histogram of LD correlation values  
# m <- ggplot(final_LD,aes(x=abs_rho))
# m <- m + geom_histogram(bin=.01) +
#   theme_bw()+
#   scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
# m
# 
# setwd("/Users/kristen/Documents/transposon_figure_data/figures")
# ggsave(filename="LD_Counts_Histogram.tiff",
#        dpi=300,
#        width=7,
#        height=2.5,
#        units="in")