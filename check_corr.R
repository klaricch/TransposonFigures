#!/usr/bin/R
# this script checks the correlation between ZERO and ONE traits
# USE: check_corr.R
library(dplyr)
library(tidyr)
library(ggplot2)

setwd("/Users/kristen/Documents/transposon_figure_data/data")
data<-read.table("T_kin_C_matrix_full_reduced.txt",header=TRUE)
names(data)
data$family <- stringr::str_split_fixed(data$trait, "_TRANS_",2)[,2]
data$method <- stringr::str_split_fixed(data$trait, "_TRANS_",2)[,1]
data$family
data$method
data<-filter(data,method=="ONE_new" | method=="ZERO_new")
data$family<-gsub('_C$','', data$family)
data<-select(data,-trait)
data<-gather(data,strain,count,AB1:WN2002)

data<-spread(data,method,count)
data<-filter(data,family!="total")
#na<-data %>% group_by(family) %>%  summarise(A=sum(is.na(ONE_new)),B=sum(is.na(ONE_new)))
data<-data %>% group_by(family) %>%  summarise(cor=round((cor.test(ONE_new,ZERO_new,method="spearman",exact=FALSE))$estimate,3))

save(data, file="ONE_ZERO_cor.Rda")
