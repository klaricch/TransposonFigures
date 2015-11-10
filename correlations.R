setwd("/Users/kristen/Documents/transposon_figure_data/data")
load("Processed_Transposon_Mappings.Rda")
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)

################################################################################################################
################################################################################################################
################################################################################################################
# COUNT FRACTION CORRELATION
#add columns
final_processed_mappings$family <- stringr::str_split_fixed(final_processed_mappings$pheno, "_TRANS_",2)[,2]
final_processed_mappings$method <- stringr::str_split_fixed(final_processed_mappings$pheno, "_TRANS_",2)[,1]
base_traits <-final_processed_mappings[(final_processed_mappings$method=="absent"| final_processed_mappings$method=="new" |final_processed_mappings$method=="reference"|final_processed_mappings$method=="ZERO_new"|final_processed_mappings$method=="ONE_new"), ]

# TAKE OUT NEW IN ABOVE!!!!!!!

peaks <- base_traits[base_traits$peak_id!="NA",]
#pull unique combinations of peak, strain,and pheno
distinct_sites <- distinct(peaks,strain,pheno)
#split into counts and fractions dataframes
#subset data
counts<-subset(distinct_sites, grepl("_C$", distinct_sites$pheno))
fractions<-subset(distinct_sites, !grepl("_C$", distinct_sites$pheno))
counts$family <- gsub("_C$" ,"",counts$family)


#remove stuff CHECK STUFF OVER HERE
counts <- filter(counts, family!="total" ) #remove totals which will only be in counts and not fractions

#initiate an empty dataframe
corr_df <- data.frame(TE=character(),
                      method=character(),
                      rho=numeric())

for (TE in unique(counts$family)){
  temp_df<-filter(counts, family==TE)
  for (way in unique(temp_df$method)){

    count_TE <- filter(counts, family==TE,method==way)
    fraction_TE <- filter(fractions, family==TE,method==way)
    #join dataframes
    all<-merge(count_TE,fraction_TE,by="strain")
    if (nrow(all)!=0){ # check if TE and method were in both counts and fractions...if not, all will be none (so far only happens with LTRCER1) 
    #spearman correlation
    correlation<-cor.test(all$value.x, all$value.y,method="spearman",exact=FALSE)
    rho_value <-correlation$estimate
    corr_df <- rbind(corr_df, data.frame(TE = TE, method = way,rho=rho_value))
  }
  }
}

m <- ggplot(corr_df, aes(x=rho))
m <- m + geom_bar(binwidth=.001) +
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=9),
        axis.text.x = element_text(colour = "black",size=9),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9))+
  guides(fill=FALSE) +
  labs(x="rho", y="Count")+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))
m
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="Correlation_Counts_Fractions.tiff",dpi=300, width=7.5,height=3.5,units="in")

################################################################################################################
################################################################################################################
################################################################################################################
#MEAN AND SD
#initiate an empty dataframe
count_mins <- data.frame(TE=character(),
                      min_p=numeric())
                      
frac_mins <- data.frame(TE=character(),
                        min_p=numeric())
                                           
#just iterate through disticnt pheno
for (phenotype in unique(base_traits$pheno)){
  temp_df<-filter(base_traits, pheno==phenotype) 
  minimum_p<-min(temp_df$ps)
  if(grepl("_C$", phenotype)==TRUE){ #counts
    count_mins <- rbind(count_mins, data.frame(TE = phenotype, min_p = minimum_p))
  }
  else if (!grepl("_C$", phenotype)==TRUE){ #fractions
    frac_mins <- rbind(frac_mins, data.frame(TE = phenotype, min_p = minimum_p))
   
  }
}

#calculate means and SDs of the minimum P values:
mean_C=mean(count_mins$min_p)
SD_C=sd(count_mins$min_p)
mean_F=mean(frac_mins$min_p)
SD_F=sd(frac_mins$min_p)

print(mean_C)
print(mean_F)
print(SD_C)
print(SD_F)

################################################################################################################
################################################################################################################
################################################################################################################
# COUNT ACTIVITY CORRELATION
activity_traits<-subset(final_processed_mappings, grepl('^(?!absent).*', final_processed_mappings$method,perl=T) &
                          grepl('^(?!reference).*', final_processed_mappings$method,perl=T) &
                          grepl('^(?!ZERO_new).*', final_processed_mappings$method,perl=T) &
                          grepl('^(?!ONE_new).*', final_processed_mappings$method,perl=T) &
                          grepl('^(?!I).*', final_processed_mappings$method,perl=T) &
                          grepl('^(?!V).*', final_processed_mappings$method,perl=T) &
                          grepl('^(?!X).*', final_processed_mappings$method,perl=T) &
                          grepl('^(?!coverage).*', final_processed_mappings$method,perl=T))
activity_traits$family <- gsub("_C$" ,"",activity_traits$family)
ac_peaks <- activity_traits[activity_traits$peak_id!="NA",]
#pull unique combinations of peak, strain,and pheno
activity_traits <- distinct(ac_peaks,strain,pheno)
#initiate an empty dataframe
final_corr <- data.frame(transposon=character(),
                      base=character(),
                      activity=character(),
                      rho=numeric())

#correlations for each TE fam, abse, trait, activity trait
for (TE in unique(counts$family)){
  temp_df<-filter(counts, family==TE)
  for (way in unique(temp_df$method)){ #no longer matching methods
    count_TE <- filter(counts, family==TE,method==way)
    for (act in unique(activity_traits$method)){
      activity_TE <- filter(activity_traits, family==TE,method==act)
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
write.table(combo_corr, "/Users/kristen/Documents/transposon_figure_data/figures/combo_corr.txt", sep="\t")



