#!/usr/bin/R
# this script generates a tree based on TEs
setwd("/Users/kristen/Documents/transposon_figure_data/tables/final_tables")
#install.packages("ggdendro")
library(ggdendro)
library(ggplot2)

data<-read.table("Supplemental_Table_S2.txt",header=TRUE,row.names=1) #reduced matrix
dm <- dist(as.matrix(data), method="manhattan")
hc<-hclust(dm)
dend<-ggdendrogram(hc,rotate=TRUE)
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="TE_Tree.tiff",
       dpi=300,
       width=4.5,
       height=7.5)

#dend <- as.dendrogram(hc)
#plot(dend)

