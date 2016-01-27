setwd("/Users/kristen/Documents/transposon_figure_data/data")
load("Processed_Transposon_Mappings.Rda")


library(dplyr)
library(tidyr)
# COUNT FRACTION CORRELATION
#add columns for family and method (will only apply to count traits eventually)
processed_mapping_df$family <- stringr::str_split_fixed(processed_mapping_df$trait, "_TRANS_",2)[,2]
processed_mapping_df$method <- stringr::str_split_fixed(processed_mapping_df$trait, "_TRANS_",2)[,1]
#pull out base (counts and fractions if present) traits for families
base_traits <-processed_mapping_df[(processed_mapping_df$method=="absent"|processed_mapping_df$method=="reference"|processed_mapping_df$method=="ZERO_new"|processed_mapping_df$method=="ONE_new"), ]

names(base_traits)
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
counts<-filter(distinct_sites,!grepl('total', family))#remove totals which will only be in counts and not fractions
unique(counts$family)

#
#
#



#initiate an empty dataframe
final_corr <- data.frame(transposon=character(),
                         base=character(),
                         base2=character(),
                         rho=numeric())

#correlations for each TE fam, abse, trait
for (TE in unique(counts$family)){ #for each count family mapped
  temp_df<-filter(counts, family==TE) #subset the count df for that family
  for (way in unique(temp_df$method)){ #for each method(absent,ref, ONE_new, ZERO_new); no longer matching methods
    count_TE <- filter(counts, family==TE,method==way) # subset the count df for that method
    for (act in unique(temp_df$method)){ # for each  trait method
      act_TE <- filter(counts, family==TE,method==act) # subset the activity df for that method and TE family
      #join dataframes
      all<-merge(count_TE,act_TE,by="strain")
      if (nrow(all)!=0){ # check if TE and method were in both counts and activity...if not, all will be none (so far only happens with LTRCER1) 
        #spearman correlation
        correlation<-cor.test(all$value.x, all$value.y,method="spearman",exact=FALSE)
        rho_value <-correlation$estimate
        final_corr <- rbind(final_corr, data.frame(transposon = TE, base = way,base2=act,rho=rho_value))
      }
    }
  }
}


final_corr<-filter(final_corr,base!=base2)

#get unique combinations of final_corr
final_corr<-mutate(final_corr, combo=paste(base,base2,sep=" "))
final_corr$uniqcombo <- sapply(strsplit(final_corr$combo, " "), function(x) 
  paste(sort(x), collapse=" ")) 

final_corr<-distinct(final_corr,transposon,uniqcombo)
final_corr$transposon<-gsub("_C$","",final_corr$transposon)
processed_mapping_df$trait<-gsub("_C$","",processed_mapping_df$trait)

final_corr$base_one <- stringr::str_split_fixed(final_corr$uniqcombo, " ",2)[,1]
final_corr$base_two <- stringr::str_split_fixed(final_corr$uniqcombo, " ",2)[,2]
final_corr<-select(final_corr,-base,-base2,-combo)

novel_TEs<-filter(final_corr,base_one=="ONE_new" & base_two=="ZERO_new")
ref_TEs<-filter(final_corr,base_one=="absent" & base_two=="reference")

final_corr<-filter(final_corr,(base_one=="ONE_new" & base_two=="ZERO_new")|(base_one=="absent" & base_two=="reference"))
thres_set<-filter(final_corr,abs(rho)>=.80)

thres_set<-mutate(thres_set,trait_recip_keep=paste(base_one,"TRANS",transposon, sep="_"))
thres_set<-mutate(thres_set,trait_recip_remove=paste(base_two,"TRANS",transposon, sep="_"))



head(processed_mapping_df)
#initiate an empty dataframe
reciprocal_removals <- data.frame(TE=character())
non_removals <- data.frame(marker=character(),
                           trait=character(),
                           peak_id=character())



for  (i in 1:nrow(thres_set)){
  temp_df<-thres_set[i,]

temp1<-filter(processed_mapping_df,trait==temp_df$trait_recip_keep,!is.na(peak_id))
temp1<-distinct(temp1,marker)
temp2<-filter(processed_mapping_df,trait==temp_df$trait_recip_remove,!is.na(peak_id))
temp2<-distinct(temp2,marker)
if(all(temp1$marker==temp2$marker)==TRUE){reciprocal_removals <- rbind(reciprocal_removals, data.frame(TE = temp_df$trait_recip_remove))} else {
  t1<-select(temp1,marker,trait,peak_id,log10p)
  t2<-select(temp2,marker,trait,peak_id,log10p)
  t3<-full_join(t1, t2, by = "marker")
  t3<-filter(t3, !is.na(log10p.x) & !is.na(log10p.y))
  if (nrow(t3)>0){
    t3<-mutate(t3,dec=ifelse(log10p.x>log10p.y,"WIN",ifelse(log10p.x==log10p.y,"TIE","LOSE")))
  #keep the default recip trait if it WINs or TIEs in number of greater pvalues between shared markers
   if (nrow(t3[t3$dec=="LOSE",]) <= nrow(t3[t3$dec=="WIN",])){
      reciprocal_removals <-rbind(reciprocal_removals, data.frame(TE = temp_df$trait_recip_remove))} else {
        reciprocal_removals <-rbind(reciprocal_removals, data.frame(TE = temp_df$trait_recip_keep))}
  } else{  non_removals<-rbind(non_removals, t1)
           non_removals<-rbind(non_removals, t2)}

}

}


### CHECK TO MAKE SURE ABSENCE COMES BEFORE REFERENCE!!!!

save(reciprocal_removals, file="reciprocal_removals.Rda")

#  t1<-select(temp1,marker,trait,peak_id,log10p)
#t2<-select(temp2,marker,trait,peak_id,log10p)
#non_removals<-rbind(non_removals, t1)
#non_removals<-rbind(non_removals, t2)

test<-filter(processed_mapping_df,trait=="reference_TRANS_MIRAGE1",!is.na(peak_id))
test2<-filter(processed_mapping_df,trait=="absent_TRANS_MIRAGE1",!is.na(peak_id))

