#!/usr/bin/R
# this script plots the distrubition of transposons (totals across all samples) according to their chromosome positions
# USE: TE_density.R

library(ggplot2)
library(grid)

setwd("/Users/kristen/Documents/transposon_figure_data/data")
summarydata <- read.table("CtCp_all_nonredundant.txt",header=TRUE)
names(summarydata)
names(summarydata)<-c("chr","start","end","TE","orientation","method","strain","class")

#3X-BIN .25MB
method_names <- list(
  'absent'="Absence",
  'new'="Insertion",
  'reference'="Reference"
)

method_labeller <- function(variable,value){
  if (variable=='method') {
    return(method_names[value])
  }else {
    return(as.character(value))
  }
}

# Add y coordinates for "phantom" points
names(summarydata)
summarydata$top <- NA
summarydata$top[summarydata$method=="absent"] <- 250
summarydata$top[summarydata$method=="new"] <- 160
summarydata$top[summarydata$method=="reference"] <- 900
levels(summarydata$class)


#revalue classes
summarydata$class <- factor(summarydata$class,
                          levels = c("dnatransposon", "retrotransposon","unknown"),
                          labels = c("DNA Transposon", "Retrotransposon", "Unknown"))

m <- ggplot(summarydata, aes(x=start/1e6,fill=class))
m <-m + geom_bar(binwidth=.25)+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))+
  facet_grid(method ~ chr,scale="free",space = "free_x",labeller=method_labeller)+
  geom_point(data = subset(summarydata, method=="absent"),aes(y=top),alpha=0) +
  geom_point(data = subset(summarydata, method=="new"),aes(y=top),alpha=0) +
  geom_point(data = subset(summarydata, method=="reference"),aes(y=top),alpha=0) +

  labs(x="Chromosome Position (Mb)", y="Number of Transposition Events")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 9, colour = "black",face="bold"),
        panel.margin = unit(.6, "lines"),
        panel.border = element_rect(fill=NA,colour = "black"),
        panel.background = element_blank(),
        axis.title=element_text(size=9),
        axis.text.y = element_text(colour = "black",size=9),
        axis.text.x = element_text(colour = "black",size=9),
        axis.ticks =element_line(colour = "black"),
        legend.title=element_blank(),
        legend.text=element_text(size=8))+
  scale_fill_manual(values = c("navy", "brown3", "darkgoldenrod2"))
m
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="Chromosome_Distribution.tiff",
       dpi=300,
       width=7.5,
       height=3.5,
       units="in")
