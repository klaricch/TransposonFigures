setwd("/Users/kristen/Documents/transposon_figure_data/data")
load("Processed_Transposon_Mappings.Rda")
load("SignificantMappings_Results_Activity.Rda")
library(dplyr)
library(tidyr)
library(data.table)
# TE COUNT FRACTION
head(final_processed_mappings$pheno)
head(final_processed_mappings)

#add columns
transposon <- stringr::str_split_fixed(final_processed_mappings$pheno, "_TRANS_",2)[,2]
final_processed_mappings$family <- transposon
caller <- stringr::str_split_fixed(final_processed_mappings$pheno, "_TRANS_",2)[,1]
final_processed_mappings$method <- caller 

#get only base traits
base_traits <-final_processed_mappings[(final_processed_mappings$method=="absent"| final_processed_mappings$method=="new" |final_processed_mappings$method=="reference"|final_processed_mappings$method=="ZERO_new"|final_processed_mappings$method=="ONE_new"), ]

# calculate smallest p value
CF<-base_traits %>%
  group_by(pheno) %>%
  summarise(Min=min(ps))

#modify/add columns
CF$transposon <- stringr::str_split_fixed(CF$pheno, "_TRANS_",2)[,2]
CF$caller <- stringr::str_split_fixed(CF$pheno, "_TRANS_",2)[,1]
CF$transposon <- gsub("_C$" ,"",CF$transposon)

#subset data
counts<-subset(CF, grepl("_C$", CF$pheno))
fractions<-subset(CF, !grepl("_C$", CF$pheno))

#change column names
setnames(counts, "Min", "C_min")
setnames(fractions, "Min", "F_min")
names(counts)

#merge transposon and caller column
counts<-mutate(counts, type=paste(transposon,caller, sep ="_"))
fractions<-mutate(fractions, type=paste(transposon,caller, sep ="_"))

#join dataframes
min<-cbind(counts,fractions,by="type")

#spearman correlation
cor.test(min$C_min, min$F_min,method="spearman",exact=FALSE)



#####################################################################################################################################
#####################################################################################################################################
#take subset to include only those traits above bonferroni
#?????????

sigs<-final_processed_mappings[final_processed_mappings$aboveBF==1,]
print(sigs$aboveBF)

#add columns
transposon <- stringr::str_split_fixed(sigs$pheno, "_TRANS_",2)[,2]
sigs$family <- transposon
caller <- stringr::str_split_fixed(sigs$pheno, "_TRANS_",2)[,1]
sigs$method <- caller 

#get only base traits
sig_base_traits <-sigs[(sigs$method=="absent"| sigs$method=="new" |sigs$method=="reference"|sigs$method=="ZERO_new"|sigs$method=="ONE_new"), ]

# calculate smallest p value
sig_CF<-sig_base_traits %>%
  group_by(pheno) %>%
  summarise(Min=min(ps))

#modify/add columns
sig_CF$transposon <- stringr::str_split_fixed(sig_CF$pheno, "_TRANS_",2)[,2]
sig_CF$caller <- stringr::str_split_fixed(sig_CF$pheno, "_TRANS_",2)[,1]
sig_CF$transposon <- gsub("_C$" ,"",sig_CF$transposon)

#subset data
sig_counts<-subset(sig_CF, grepl("_C$", sig_CF$pheno))
sig_fractions<-subset(sig_CF, !grepl("_C$", sig_CF$pheno))

#change column names
setnames(sig_counts, "Min", "C_min")
setnames(sig_fractions, "Min", "F_min")
names(sig_counts)

#merge transposon and caller column
sig_counts<-mutate(sig_counts, type=paste(transposon,caller, sep ="_"))
sig_fractions<-mutate(sig_fractions, type=paste(transposon,caller, sep ="_"))

#join dataframes
sig_min<-cbind(sig_counts,sig_fractions,by="type")

#spearman correlation
cor.test(sig_min$C_min, sig_min$F_min,method="spearman",exact=FALSE)




# calculate mean and standard deviation
C_mean<-mean(sig_min$C_min)
C_sd<-sd(sig_min$C_min)
print(C_mean)
print(C_sd)

F_mean<-mean(sig_min$F_min)
F_sd<-sd(sig_min$F_min)
print(F_mean)
print(F_sd)

#####################################################################################################################################
#####################################################################################################################################
#heatmap
activity_traits<-subset(final_processed_mappings, grepl('^(?!absent).*', final_processed_mappings$method,perl=T) &
                          grepl('^(?!reference).*', final_processed_mappings$method,perl=T) &
                          grepl('^(?!ZERO_new).*', final_processed_mappings$method,perl=T) &
                          grepl('^(?!ONE_new).*', final_processed_mappings$method,perl=T) &
                          grepl('^(?!I).*', final_processed_mappings$method,perl=T) &
                          grepl('^(?!V).*', final_processed_mappings$method,perl=T) &
                          grepl('^(?!X).*', final_processed_mappings$method,perl=T) &
                          grepl('^(?!coverage).*', final_processed_mappings$method,perl=T))
head(activity_traits)
# activity traits are already only associated with count, not fraction data

# calculate smallest p value  
AT<-activity_traits %>%
  group_by(pheno) %>%
  summarise(Min=min(ps))

#modify/add columns
AT$transposon <- stringr::str_split_fixed(AT$pheno, "_TRANS_",2)[,2]
AT$caller <- stringr::str_split_fixed(AT$pheno, "_TRANS_",2)[,1]
AT$transposon <- gsub("_C$" ,"",AT$transposon)

#change column names
setnames(AT, "Min", "AT_min")
setnames(AT, "pheno", "pheno2")

#merge transposon and caller column
AT<-mutate(AT, type=paste(transposon,caller, sep ="_"))
print(AT$type)

#join dataframes
heat<-merge(counts,AT,by="transposon",allow.cartesian=TRUE)

#initiate an empty dataframe
corr_df <- data.frame(callerA=character(),
                 callerB=character(),
                 rho=numeric())

#spearman correlation for each combination of caller.x and caller.y
for (i in unique(heat$caller.y)){
  print(i)
  heatsub<-heat[heat$caller.y==i,]
  for (x in unique(heatsub$caller.x)){
    heatdoublesub<-heatsub[heatsub$caller.x==x,]
    correlation<-cor.test(heatdoublesub$C_min, heatdoublesub$AT_min,method="spearman",exact=FALSE)
    rho_value <-correlation$estimate
    corr_df <- rbind(corr_df, data.frame(callerA = x, callerB = i,rho=rho_value))
  } 
}

#generate the heatmap
p <- ggplot(corr_df, aes(callerA, callerB)) + 
  geom_tile(aes(fill = rho),colour = "white") + #white here is the lines separating the tiles
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.ticks = element_blank())+
  scale_fill_gradient(low = "gold", high = "red3")+
  labs(x="Base Trait", y= "Activity Trait")
p

setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="Heat_Map_Activity.tiff",
       dpi=300,
       width=7,
       height=2.5,
       units="in")

