#!/usr/bin/R
# this script:
# calculates the Spearman rank correlations between raw counts traits and the corresponding fraction traits 
# calculates the mean and SD of the min p values for count and fraction traits
# calcuates the average rho and SD between all combinations of counts and activity traits 
# USE: correlations.R

setwd("/Users/kristen/Documents/transposon_figure_data/data")
load("Processed_Transposon_Mappings_with activity.Rda")
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)

################################################################################################################
################################################################################################################
################################################################################################################
# COUNT FRACTION CORRELATION
#add columns for family and method (will only apply to count traits eventually)
processed_mapping_df$family <- stringr::str_split_fixed(processed_mapping_df$trait, "_TRANS_",2)[,2]
processed_mapping_df$method <- stringr::str_split_fixed(processed_mapping_df$trait, "_TRANS_",2)[,1]
#pull out base (counts and fractions if present) traits for families
base_traits <-processed_mapping_df[(processed_mapping_df$method=="absent"|processed_mapping_df$method=="reference"|processed_mapping_df$method=="ZERO_new"|processed_mapping_df$method=="ONE_new"), ]
#remove "NA" peak_ids and alleles
peaks <- filter(base_traits,!is.na(peak_id))
peaks <- filter(peaks,!is.na(allele))
#pull unique combinations of  strain and trait
distinct_sites <- distinct(peaks,strain,trait)
#subset data to get only count traits (will remove fraction traits if present) 
counts<-subset(distinct_sites, grepl("_C$", distinct_sites$trait))
# remove trailing "_C" from family names
counts$family <- gsub("_C$" ,"",counts$family)
#remove "total" traits
counts <- filter(counts, family!="total" ) #remove totals which will only be in counts and not fractions

################################################################################################################
################################################################################################################
################################################################################################################
# COUNT ACTIVITY CORRELATION
#pull out only activity traits
activity_traits<-subset(processed_mapping_df, grepl('^no_', processed_mapping_df$trait))
activity_traits$family <- gsub("_C$" ,"",activity_traits$family)
ac_peaks <- activity_traits[activity_traits$peak_id!="NA",]
#pull unique combinations of strain and trait
activity_traits <- distinct(ac_peaks,strain,trait)
#initiate an empty dataframe
final_corr <- data.frame(transposon=character(),
                      base=character(),
                      activity=character(),
                      rho=numeric())

#correlations for each TE fam, abse, trait, activity trait
for (TE in unique(counts$family)){ #for each count family mapped
  temp_df<-filter(counts, family==TE) #subset the count df for that family
  for (way in unique(temp_df$method)){ #for each method(absent,ref, ONE_new, ZERO_new); no longer matching methods
    count_TE <- filter(counts, family==TE,method==way) # subset the count df for that method
    for (act in unique(activity_traits$method)){ # for each activity trait method
      activity_TE <- filter(activity_traits, family==TE,method==act) # subset the activity df for that method and TE family
      #join dataframes
      all<-merge(count_TE,activity_TE,by="strain")
      if (nrow(all)!=0){ # check if TE and method were in both counts and activity...if not, all will be none (so far only happens with LTRCER1) 
      #spearman correlation
      correlation<-cor.test(all$value.x, all$value.y,method="spearman",exact=FALSE)
      rho_value <-correlation$estimate
      final_corr <- rbind(final_corr, data.frame(transposon = TE, base = way,activity=act,rho=rho_value))
      }
    }
  }
}

#initiate dataframe
combo_corr <- data.frame(method1=character(),
                         method2=character(),
                         mean=numeric(),
                         SD=numeric())

#averages and SD per base/activity combo:
for (M1 in unique(final_corr$base)){
  for (M2 in unique(final_corr$activity)){
    combo<-filter(final_corr, base==M1,activity==M2)
    complete_combo<-combo[complete.cases(combo[,4]),] #exclude NAs when calculating means
    mean_rho<-mean(complete_combo$rho)
    SD_rho<-sd(complete_combo$rho)
    combo_corr <- rbind(combo_corr, data.frame(method1 = M1, method2 = M2,mean=mean_rho,SD=SD_rho))
  }
}

#reverse sort
combo_corr <- combo_corr[order(-combo_corr$mean) , ]
#write output to table
write.table(final_corr, "/Users/kristen/Documents/transposon_figure_data/figures/final_corr.txt", sep="\t")
write.table(combo_corr, "/Users/kristen/Documents/transposon_figure_data/figures/combo_corr.txt", sep="\t")



