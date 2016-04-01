library(cegwas)
library(dplyr)
library(ggplot2)
setwd("/Users/kristen/Documents/transposon_figure_data/data")

#ZERO_new_TRANS_CELETC2_C
#ZERO_new_TRANS_LINE2C_C
#ZERO_new_TRANS_Tc5B_C
#ZERO_new_TRANS_WBTransposon00000637_C
#ZERO_new_TRANS_NeSL-1_C


df<-read.table("compare_maps_old.txt",header=TRUE)
df<-read.table("compare_maps_new.txt",header=TRUE)
rec<-read.table("recips.txt",header=TRUE)


CELETC2<-filter(df, trait=="ZERO_new_TRANS_CELETC2_C")
temp <- cegwas_map(CELETC2)
manplot(temp)

LINE2C<-filter(df, trait=="ZERO_new_TRANS_LINE2C_C")
temp <- cegwas_map(LINE2C)
manplot(temp)

Tc5B<-filter(df, trait=="ZERO_new_TRANS_Tc5B_C")
temp <- cegwas_map(Tc5B)
manplot(temp)

WBTransposon00000637<-filter(df, trait=="ZERO_new_TRANS_WBTransposon00000637_C")
temp <- cegwas_map(WBTransposon00000637)
manplot(temp)

NeSL<-filter(df, trait=="ZERO_new_TRANS_NeSL-1_C")
temp <- cegwas_map(NeSL)
manplot(temp)

##recips

CEL_rec<-filter(rec, trait=="ONE_new_TRANS_CELETC2_C")
temp <- cegwas_map(CEL_rec)
manplot(temp)



LINE_rec<-filter(rec, trait=="ONE_new_TRANS_LINE2C_C")
temp <- cegwas_map(LINE_rec)
manplot(temp)
#
#
#
#

setwd("/Users/kristen/Documents/transposon_figure_data/manual")
df<-read.table("ONES.txt",header=TRUE)

MIRAGE1<-filter(df, trait=="ONE_new_TRANS_MIRAGE1")
temp1 <- cegwas_map(MIRAGE1)
manplot(temp1)
MIRAGE1
MIRAGE2<-filter(df, trait=="ZERO_new_TRANS_MIRAGE1")
temp2 <- cegwas_map(MIRAGE2)
manplot(temp2)
MIRAGE2
Tc1
Tc1<-filter(df, trait=="ONE_new_TRANS_Tc1")
temp2 <- cegwas_map(Tc1,BF=4)
manplot(temp2)
?cegwas_map
Tc6<-filter(df, trait=="ONE_new_TRANS_Tc6")
temp3 <- cegwas_map(Tc6)
manplot(temp3)
