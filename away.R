#!/usr/bin/R
# this script checks which QTL are located more than 100 SNPs away from the physical location of
# the TE family that was mapped and outputs such QTL to a new dataframe
# USE: away.R

library(dplyr)
library(ggplot2)
library(data.table)
library(grid)
library(stringr)
library(gridExtra)
library(stringr)
setwd("/Users/kristen/Documents/transposon_figure_data/data")
load("Processed_Transposon_Mappings.Rda")
load("prune_final_traits_to_this_set.Rda")

final_processed_mappings <- dplyr::filter(final_processed_mappings, (pheno %in% traits))


# read in positions of all TEs
positions <- read.table("CtCp_all_nonredundant.txt",header=TRUE)
names(positions)<-c("chr","start","end","TE","orientation","method","strain","class")
#create TE family column
positions$family<- stringr::str_split_fixed(positions$TE, regex("_(non-)?reference"),2)[,1]
positions$family<- paste(stringr::str_split_fixed(positions$family, "_",4)[,3],stringr::str_split_fixed(positions$family, "_",4)[,4],sep="_")
positions$family <- gsub("_$" ,"",positions$family)
positions$family <- gsub("_non-reference(.*)$" ,"",positions$family)

#pull out only position traits from mappings dataframe
position_traits<-subset(final_processed_mappings,
                        grepl('^I', final_processed_mappings$pheno) |
                          grepl('^V', final_processed_mappings$pheno) |
                          grepl('^X', final_processed_mappings$pheno))
#create family column
position_traits$family  <- paste(stringr::str_split_fixed(position_traits$pheno, "_",4)[,3],stringr::str_split_fixed(position_traits$pheno, "_",4)[,4],sep="_")
position_traits$family <- gsub("_$" ,"",position_traits$family)
position_traits$family <- gsub("_non-reference(.*)$" ,"",position_traits$family)

#pull out unique phenotypes and peak ids
#distinct_sites <- distinct(position_traits,pheno,peak_id)
#distinct_sites <- distinct_sites[distinct_sites$peak_id!="NA",] #remove NA peaks
#dd[with(dd, order(-z, b)), ]

distinct_sites <- distinct(position_traits,pheno,chr,peak_id,pos)
distinct_sites <- distinct_sites[order(pheno,peak_id, -ps),] #sort by ascending pheno and descending pvalue
distinct_sites <- distinct_sites[distinct_sites$peak_id!="NA",] #remove NA peaks
distinct_sites <- distinct(distinct_sites,pheno,peak_id)

#distinct_sites <- distinct_sites[distinct_sites$peak_id!="NA",] #remove NA peaks
#newdata <- mtcars[order(mpg),]
#pull out all unique SNPs used for the mappings
#WANT DESCEDING
SNP<- distinct(final_processed_mappings,chr,pos)


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
  chr=row$chr
  pos=row$pos
  #chr=row[,row$chr]
  if (chr=="I"){
    one<-c(one,pos)
  }
  else if (chr=="II"){
    two<-c(two,pos)
  }
  else if (chr=="III"){
    three<-c(three,pos)
  }
  else if (chr=="IV"){
    four<-c(four,pos)
  }
  else if (chr=="V"){
    five<-c(five,pos)
  }
  else if (chr=="X"){
    six<-c(six,pos)
  }
  
}

one<-sort(one)
two<-sort(two)
three<-sort(three)
four<-sort(four)
five<-sort(five)
six<-sort(six)

#fins max length of each list
max_one=length(one)
max_two=length(two)
max_three=length(three)
max_four=length(four)
max_five=length(five)
max_six=length(six)

#initiate an empty dataframe
away<- as.data.frame(matrix(0, ncol = 27, nrow = 0))
names(away) <-c(colnames(position_traits))


###REMOVE LATER
#distinct_sites<-distinct_sites[distinct_sites$pheno=="II_507265_LINE2C",]
##







#iterate through unique QTL peaks
for (phenotype in unique(distinct_sites$pheno)){
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
  locations<-distinct(locations,chr,start,end,family)
  for(i in 1:nrow(locations)) { #add all positions of a TE to a chromosome list
    row <- locations[i,]
    TE_chr=row$chr
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
  
  pheno_subset=distinct_sites[distinct_sites$pheno==phenotype,] #pull out mappings realted to that phenotype
  #need a for each row in subset here
  for(i in 1:nrow(pheno_subset)) {
    trig=FALSE #start trig at FALSE, if the TE is within 100 SNPs of the QTL, change trig to TRUE
    row <- pheno_subset[i,]
    chr<-row$chr
    bp<-row$pos
    
    #check which chromosome it's in
    if (chr=="I"){
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
          break
          }
        #if not----> rbind
      }
    }
    
    
    #TWO
    if (chr=="II"){
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
          break
          }
        #if not----> rbind
      }
    } 
    
    #THREE
    if (chr=="III"){
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
          break
          }
        #if not----> rbind
      }
    } 
    
    ##FOUR
    if (chr=="IV"){
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
          break
          }
        #if not----> rbind
      }
    } 
    
    #FIVE
    if (chr=="V"){
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
          break
          }
        #if not----> rbind
      }
    } 
    
    #SIX
    if (chr=="X"){
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



######MEDIAN Differences
library(tidyr)
#get medians of positions traits that are in away df:
median_df <- position_traits[position_traits$pheno %in% away$pheno,]
median_df <- median_df[median_df$peak_id !="NA",]
median_df <- distinct(median_df,pheno,strain)
#calcualte median for each phenotype allele group
median_df <- median_df %>% group_by(pheno,allele) %>% summarise(med=median(value,na.rm=TRUE))
median_df<-spread(median_df, allele,med)
colnames(median_df)<-c("pheno","ref","alt")
median_df<-filter(median_df,ref!=alt)
save(median_df, file="median_phenos.Rda")
##########
#CAN REMOVE ALL OF BELOW AFTER CHECKS
nrow(away)
length(unique(away$pheno))
nrow(distinct_sites)
length(unique(distinct_sites$pheno))
length(unique(final_processed_mappings$pheno))
length(unique(position_traits$pheno))

hmmmm<-away[away$pheno=="I_10499869_WBTransposon00000581",]
hmmmm2<-distinct_sites[distinct_sites$pheno=="I_10499869_WBTransposon00000581",]

###################
tail(away)
hmmmm<-away[away$pheno=="II_3406485_MARINER2_non-reference__temp_sr_10950_120.5_-_new_CB4852",]
aaa<-positions[positions$family=="LINE2C",]
bbb<-aaa[aaa$chr=="II",]

indexcc <-which(five==3598054) 
print(indexcc)
print(five[410])
print(five[412])
print(five[310])
print(five[510])
print(length(one))
print(length(two))
print(length(three))
print(length(four))
print(length(five))
print(length(six))
print(TE_five)


indexcc <-which(three==10823037) 
print(indexcc)

indexcc <-which(three==11499581) 
print(indexcc)

print(three[1071])
print(three[1271])
print(three[1069])
print(three[0971])
print(three[1171])




indexcc <-which(two==12171121) 
print(indexcc)

print(two[0903])
print(two[1103])
print(two[1002])
print(two[1004])





print(TE_three)


check<- final_processed_mappings[final_processed_mappings$pheno=="X_9196067_CER10-I_CE",]
check2<-check[check$peak_id!="NA",]









?median
load("20151025_Processed_Transposon_Mappings.Rda")
check<-filter(final_processed_mappings,pheno=="I_10263384_LINE2C",peak_id!="NA")
check<-filter(new_mappings,pheno=="I_10263384_LINE2C",peak_id!="NA")
unique(check$value)
print(check$value)

check2<-distinct(check,pheno,strain,value)
print((check2$value))
unique(check2$value)
