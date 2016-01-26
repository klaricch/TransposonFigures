#!/usr/bin/R
# this script checks which reciporal traits for each family both mapped and outputs dfs for thos that didn't (these contian "NAs" at some sites)
library(dplyr)
library(tidyr)
setwd("/Users/kristen/Documents/transposon_figure_data/data")
load("count_QTL.Rda")
count_QTL$trait<-gsub("_\\d+$","",count_QTL$trait)
count_QTL$family <- stringr::str_split_fixed(count_QTL$trait, "_TRANS_",2)[,2]
count_QTL$method <- stringr::str_split_fixed(count_QTL$trait, "_TRANS_",2)[,1]
count_QTL<-select(count_QTL,method,family)
count_QTL<-distinct(count_QTL)
test<-spread(count_QTL,method,"h")
absent<-filter(count_QTL,method=="absent")
reference<-filter(count_QTL,method=="reference")
ONE_new<-filter(count_QTL,method=="ONE_new")
ZERO_new<-filter(count_QTL,method=="ZERO_new")

ref<-full_join(absent, reference, by = "family")
novel<-full_join(ONE_new, ZERO_new, by = "family")


refF<-filter(ref,!complete.cases(ref))
novelF<-filter(novel,!complete.cases(novel))
