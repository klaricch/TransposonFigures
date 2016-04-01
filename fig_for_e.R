#!/usr/bin/R

#Please prepare the following figures for me before next Friday (March 25th). All fonts should be minimum size of 18 and bold. 
#1. TE distribution colored by DNA and retro TEs across the genome (10"wide by 4"high)
#2. Please let me know the number of TEs of each class and rough number of families
#3. Histogram of strain counts for some TEs that map by GWA (5"wide by 5"high)
#4. Manhattan plots for those TE traits that map by GWA (10"wide by 4"high)

#!/usr/bin/R


library(dplyr)
library(ggplot2)
library(data.table)
library(grid)
library(stringr)
library(gridExtra)
library(tidyr)
library(scales)
library(gtable)
library(cowplot)

#FIGURE1
setwd("/Users/kristen/Documents/transposon_figure_data/data")
summarydata <- read.table("CtCp_all_nonredundant.txt")
names(summarydata)
names(summarydata)<-c("chr","start","end","TE","orientation","method","strain","class")

#3X-BIN .25MB
summarydata <- distinct(summarydata, chr,start,method, orientation,class)

# Add y coordinates for "phantom" points
names(summarydata)
summarydata$top <- NA
summarydata$top[summarydata$method=="absent"] <- 6
summarydata$top[summarydata$method=="new"] <- 30
summarydata$top[summarydata$method=="reference"] <- 8
levels(summarydata$class)

#revalue classes
summarydata$class <- factor(summarydata$class,
                            levels = c("dnatransposon", "retrotransposon","unknown"),
                            labels = c("DNA Transposon", "Retrotransposon", "Unknown"))

#revalue methods
summarydata$method <- factor(summarydata$method,
                             levels = c("new","reference","absent"),
                             labels = c("Insertion", "Reference","Absence"))



m <- ggplot(summarydata, aes(x=start/1e6,fill=class))
m <-m + geom_histogram(binwidth=.25)+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))+
  facet_grid(method ~ chr,scale="free",space = "free_x")+
  geom_point(data = subset(summarydata, method=="Absence"),aes(y=top),alpha=0) +
  geom_point(data = subset(summarydata, method=="Insertion"),aes(y=top),alpha=0) +
  geom_point(data = subset(summarydata, method=="Reference"),aes(y=top),alpha=0) +
  
  labs(x="Chromosome Position (Mb)", y="Number of Sites")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 18, colour = "black",face="bold"),
        panel.margin = unit(.25, "lines"),
        panel.margin.y = unit(1, "lines"),
        panel.border = element_rect(fill=NA,colour = "black"),
        panel.background = element_blank(),
        axis.title=element_text(size=18,face="bold"),
        axis.text.y = element_text(colour = "black",size=18,face="bold"),
        axis.text.x=element_blank(),
        #axis.text.x = element_text(colour = "black",size=9),
        axis.ticks =element_line(colour = "black"),
        legend.title=element_blank(),
        legend.position="none",
        legend.key.size=unit(.1,"cm"),
        legend.text=element_text(size=18,,face="bold"))+
  #scale_fill_manual(values = c("navy", "brown3", "darkgoldenrod2"))
  scale_fill_manual(values = c("DNA Transposon" = "navy", "Retrotransposon"="brown3","Unknown"="darkgoldenrod2"))

m
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="fige_Chromosome_Distribution.tiff",
       dpi=300,
       width=10,
       height=5,
       units="in")








#FIGURE2
setwd("/Users/kristen/Documents/transposon_figure_data/data")
summarydata <- read.table("T_kin_C_matrix_full.txt",header=TRUE)
#remove ZERO_new traits
summarydata<-subset(summarydata,!grepl('^ZERO_new', summarydata$trait))
summarydata<-subset(summarydata,!grepl('^coverage', summarydata$trait))
#clean trait names
summarydata$trait <- gsub("_C$" ,"",summarydata$trait)
summarydata$trait <- gsub("^ONE_new" ,"new",summarydata$trait)

classdata<- read.table("CtCp_all_nonredundant.txt")
names(classdata)<-c("chr","start","end","TE","orientation","method","strain","class")

# add te class info to summarydata(new_TRANS_end_tes will be removed)
classdata$family<- stringr::str_split_fixed(classdata$TE, regex("_(non-)?reference"),2)[,1]
classdata$family<- paste(stringr::str_split_fixed(classdata$family, "_",4)[,3],stringr::str_split_fixed(classdata$family, "_",4)[,4],sep="_")
classdata$family <- gsub("_$" ,"",classdata$family)
classdata$family <- gsub("_non-reference(.*)$" ,"",classdata$family)
classdata<-mutate(classdata, trait=paste(method,"TRANS",family,sep="_"))
class_subset <- classdata %>% distinct(family,class) %>% select(family,class)
summarydata$family<- stringr::str_split_fixed(summarydata$trait, "_TRANS_",2)[,2]
summarydata <-merge(summarydata, class_subset, by="family")
summarydata<-select(summarydata, -family)
unique(class_subset$method)


#summarydata <-merge(summarydata, class_subset, by="trait")

#nrow(distinct(classdata,method,TE))
#x<-merge(summarydata, class_subset, by="family")
#Z<-filter(summarydata, !(summarydata$trait %in% x$trait))
#unique(summarydata$trait)
#tail(classdata)
#unique(classdata$family)
#test1<-filter(classdata, family=="CEMUDR2")
#names(summarydata)

summarydata$trait
summarydata<-gather(summarydata, "sample","value",2:(ncol(summarydata)-1))
#summarydata<-rename(summarydata,total_tes=value)
tail(summarydata)

#new column that specifies what caller was used
summarydata$method<- stringr::str_split_fixed(summarydata$trait, "_TRANS_",2)[,1]
#new column that specifies TE family
summarydata$transposon<- stringr::str_split_fixed(summarydata$trait, "_TRANS_",2)[,2]
#summarydata<-filter(summarydata,transposon=="total")

families<-summarydata

#reformat the data
summarydata <- summarydata %>% group_by(sample,class, method) %>% summarise(XX=sum(value,na.rm=TRUE))
hist_data<-summarydata
total_absence<-filter(summarydata,method=="absent")
total_reference<-filter(summarydata,method=="reference")
total_insertion<-filter(summarydata,method=="new")

#D_max_insertions D_max_references D_max_absences



#Histogram Versions
hist_data$method <- factor(hist_data$method,
                           levels = c("new", "reference", "absent"),
                           labels = c("Insertion", "Reference", "Absence"))

hist_data<-filter(hist_data,class=="dnatransposon",method=="Insertion")
#dim(hist_data)
#sub_hist<-hist_data %>% group_by(method) %>% summarize(MM=max(XX))
#hist_data<-merge(hist_data,sub_hist,by="method")

m <- ggplot(hist_data, aes(x=XX,fill=class))
m <-m + geom_histogram(binwidth=2)+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))+
  #facet_wrap(~method,scale="free")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 18, colour = "black",face="bold"),
        panel.margin = unit(.25, "lines"),
        panel.border = element_rect(fill=NA,colour = "black"),
        panel.background = element_blank(),
        axis.title=element_text(size=18,face="bold"),
        axis.text.y = element_text(colour = "black",size=18,face="bold"),
        axis.text.x = element_text(colour = "black",size=18,face="bold"),
        axis.ticks =element_line(colour = "black"),
        legend.title=element_blank(),
        legend.position="none")+
  scale_fill_manual(values = c('dnatransposon' = "navy", "retrotransposon"="brown3","unknown"="goldenrod"))+
  labs(x="Transposon Insertions Per Strain", y="Count",title=" ")
m

#pull out y max of panels
panel1<-filter(ggplot_build(m)$data[[1]],PANEL==1)
max1<-max(panel1$y)

m <- m + geom_point(data = hist_data,aes(y=1.05*max1),alpha=0)

m
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="fige_Histogram_TE_per_Strain.tiff",
       dpi=300,
       width=5,
       height=5,
       units="in")

#FIGURE3

setwd("/Users/kristen/Documents/transposon_figure_data/data")
load("Processed_Transposon_Mappings_SUBSET2.Rda")


unique(processed_mapping_df$trait)
load("count_QTL.Rda")

# pull unique combos, remove strain column(don't need specific strain info at this point)
processed_mapping_df<- processed_mapping_df %>% distinct(trait,marker,strain)

#create family and method columns
processed_mapping_df$family <- stringr::str_split_fixed(processed_mapping_df$trait, "_TRANS_",2)[,2]
processed_mapping_df$method <- stringr::str_split_fixed(processed_mapping_df$trait, "_TRANS_",2)[,1]

#read in position data and create family column
positions <- read.table("CtCp_all_nonredundant.txt")
names(positions)<-c("CHROM","start","end","TE","orientation","method","strain","class")
positions$family<- stringr::str_split_fixed(positions$TE, regex("_(non-)?reference"),2)[,1]
positions$family<- paste(stringr::str_split_fixed(positions$family, "_",4)[,3],stringr::str_split_fixed(positions$family, "_",4)[,4],sep="_")
positions$family <- gsub("_$" ,"",positions$family)
positions$family <- gsub("_non-reference(.*)$" ,"",positions$family)

#select traits above BF.....this step not needed, double checking everything is above BF
selection<-filter(processed_mapping_df, log10p > BF)

#extract the count base traits
base_traits <-selection[(selection$method=="absent"| selection$method=="new" |selection$method=="reference"|selection$method=="ZERO_new"|selection$method=="ONE_new"), ]
counts<-subset(base_traits, grepl("_C$", base_traits$family))
counts$family <- gsub("_C$" ,"",counts$family)

processed_mapping_df <- distinct(select(processed_mapping_df, -strain,-allele,-value))
processed_mapping_df<- processed_mapping_df %>% distinct(trait,marker)

#pull out only position traits from mappings dataframe
position_traits<-subset(selection,
                        grepl('^I', selection$trait) |
                          grepl('^V', selection$trait) |
                          grepl('^X', selection$trait))

#create family column
position_traits$family  <- paste(stringr::str_split_fixed(position_traits$trait, "_",4)[,3],stringr::str_split_fixed(position_traits$trait, "_",4)[,4],sep="_")
position_traits$family <- gsub("_$" ,"",position_traits$family)
position_traits$family <- gsub("_non-reference(.*)$" ,"",position_traits$family)

# add position TRAIT_col family info to processed_mapping_df
processed_mapping_df<-processed_mapping_df %>%mutate(family = ifelse(processed_mapping_df$trait %in% position_traits$trait, (paste(stringr::str_split_fixed(processed_mapping_df$trait, "_",4)[,3],stringr::str_split_fixed(processed_mapping_df$trait, "_",4)[,4],sep="_")), processed_mapping_df$family))

selection<-counts

#strip count marker and remnant marks from dataframes
selection$trait <- gsub("_C$" ,"",selection$trait)
processed_mapping_df$trait <- gsub("_C$" ,"",processed_mapping_df$trait)
processed_mapping_df$family <- gsub("_C$" ,"",processed_mapping_df$family)
processed_mapping_df$family <- gsub("_$" ,"",processed_mapping_df$family)
processed_mapping_df$family <- gsub("_non-reference(.*)$" ,"",processed_mapping_df$family)

processed_mapping_df<-mutate(processed_mapping_df,ID=paste(trait,peak_id,sep="_"))
copy<-processed_mapping_df
processed_mapping_df<-mutate(processed_mapping_df,SNP_col=ifelse(ID %in% count_QTL$trait,"PASS","FAIL" ))

count_QTL<-mutate(count_QTL, trait2=gsub("_\\d+$","",trait)) 
selection <- filter(selection, (trait %in% count_QTL$trait2))

processed_mapping_df<-filter(processed_mapping_df,CHROM != "MtDNA")
class_subset<- positions %>% distinct(class,family) %>% select(class,family)
selection <-merge(selection, class_subset, by="family")
selection<-arrange(selection,class,family,method)

count<-0
for (i in unique(selection$trait)){
  specific_trait<- processed_mapping_df[processed_mapping_df$trait == i, ]
  empty <-specific_trait[specific_trait$method==NA,]
  #specific_trait_mx <- max(specific_trait$log10p)
  class_TE<-unique(filter(selection,trait==i)$class)
  pvalues<-filter(specific_trait,log10p !="Inf") #
  specific_trait_mx <- max(pvalues$log10p) #
  TE<-specific_trait$family[1]
  rect_data<-filter(specific_trait,SNP_col=="PASS")
  plot_title<-gsub(".*_TRANS_","",i)
  plot_title<-gsub("_CE$","",plot_title)
  plot_title<-gsub("WBTransposon","WBT",plot_title)

  A<- processed_mapping_df %>%
    filter(trait == i)%>%
    .[order(.$peak_id,na.last=FALSE),]%>% 
    ggplot(.)+
    aes(x=POS/1e6,y=log10p,fill=BF)+ #fill to get legend
    geom_rect(data=rect_data,mapping=aes(xmin=startPOS/1e6, xmax=endPOS/1e6, ymin=0, ymax= Inf),fill="thistle1", alpha=1) +
    geom_point(aes( color=ifelse(log10p> BF & SNP_col=="PASS", 'red', 'black')),size=1)+
    facet_grid(.~CHROM,scale="free_x",space = "free_x") + #scale_color_identity() +
    geom_hline(aes(yintercept=BF),color="grey60",linetype="dashed")+
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(color="black", size=0.5, linetype="solid", fill=NA),
          panel.margin = unit(.6, "lines"),
          panel.background = element_blank(),
          axis.ticks =element_line(colour = "black"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(colour = "black",size=18,face="bold"),
          axis.title.y = element_text(size=18,face="bold",colour=ifelse(class_TE=="dnatransposon","navy",ifelse(class_TE=="retrotransposon","brown3","darkgoldenrod2"))),
          axis.title=element_text(size=18,face="bold"),
          plot.margin=unit(c(.05,.30,-.5,.30), "cm"),
          legend.title = element_text(size = 18, colour = ifelse(class_TE=="dnatransposon","navy",ifelse(class_TE=="retrotransposon","brown3","darkgoldenrod2")), angle = 270),
          legend.text = element_blank(),
          legend.key.size=unit(0,"cm"),
          legend.key = element_rect(colour = "pink"),
          legend.position=('right'))+
    labs(x="",y="",colour="red")+
    scale_color_identity()+
    scale_fill_continuous(name=plot_title)+
    scale_y_continuous(breaks= pretty_breaks(),expand=c(0,0),limits=c(0,specific_trait_mx+.075*specific_trait_mx),labels = function(x) format(x,width = 4))
  A
  
  if (count==0){B<-A+theme(strip.background = element_rect(fill = "white"),
    strip.text.x = element_text(size = 18, colour = "black",face="bold"));
                  first<-B}
  if (count==1){second<-A}
  if (count==2){third<-A}
  if (count==3){fourth<-A}
  if (count==4){fifth<-A}
  if (count==5){sixth<-A}
  count<-count+1
}

a_all<-plot_grid(first,second,third,fourth,fifth,sixth,ncol=1) #ZER)_new_TRANS_NeSL-1_C no  longer in here so don't need fifth 
label<-expression(bold(-log["10"](p)))

a_all<- a_all + draw_label(label, x = .04, y = 0.5, hjust = .5, vjust = .5,
                           fontfamily = "", fontface = "bold", colour = "black", size = 18,
                           angle = 90, lineheight = 0.9, alpha = 1)

df <- data.frame(1,2)
blank_plot<-ggplot(df,aes(x=1,y=1)) + geom_point(color="white") + theme(axis.line=element_blank(),axis.text =element_blank(),axis.ticks =element_blank(),axis.title =element_blank(),panel.background = element_blank(),panel.grid = element_blank())
a_all<-plot_grid(a_all,blank_plot,ncol=1,rel_heights = c(1, .03))


a_all<- a_all + draw_label("Chromosome Position (Mb)", x = .5, y = 0.020, hjust = .5, vjust = .5,
                           fontfamily = "", fontface = "bold", colour = "black", size = 18,
                           angle = 0, lineheight = 0.9, alpha = 1)

setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="fige_five_trait_QTL.tiff",
       dpi=300,
       width=10,
       height=10,
       units="in")





#FIGURE4


library(dplyr)
library(ggplot2)
library(data.table)
library(grid)
library(stringr)
library(gridExtra)
library(tidyr)
library(scales)
library(gtable)



setwd("/Users/kristen/Documents/transposon_figure_data/data")
load("Processed_Transposon_Mappings_2.Rda")

load("count_QTL.Rda")
pg<-read.table("paragraphs.txt",sep="\t",header=TRUE)

#PxG function
hm<-processed_mapping_df
hm<-distinct(processed_mapping_df, trait,strain,peak_id)
hm$allele <- factor(hm$allele,
                    levels = c(-1,1),
                    labels = c("REF", "ALT"))
c<-filter(processed_mapping_df, log10p> BF)
e<-filter(processed_mapping_df,!is.na(peak_id), trait=="absent_TRANS_CER1_C")
f<-distinct(e,POS,peak_id)

gwasPxG <- function(trt,specific_peak){
  hm %>%
    filter(trait==trt,peak_id==specific_peak,!is.na(allele))%>%
    ggplot(.)+
    aes(x=allele,y = value,fill=as.factor(allele))+
    geom_boxplot(outlier.shape=NA,size =.5,color="gray52")+
    geom_point(size = 1, alpha = .8,position=position_jitter(w=.4,  h=.025),na.rm=TRUE)+
    theme_bw()+
    theme(axis.text.x = element_text(size=18, , color="black",face="bold"),
          axis.text.y = element_text(size=18, , color="black",face="bold"),
          axis.title.x = element_text(size=18, , color="black",face="bold"),
          axis.title.y = element_text(size=18,  color="black",face="bold",vjust=1),
          strip.text.x = element_text(size=18, , color="black",face="bold"),
          strip.text.y = element_text(size=18, , color="black",face="bold"),
          plot.title = element_text(size=18,face="bold",  vjust=1),
          legend.title = element_text(size=18,face="bold"),
          panel.border = element_rect(size=1, colour = "black"),
          plot.margin = unit(c(.05,.05,.05,.05), "cm"),
          legend.position = "none")+
    scale_y_continuous(breaks= pretty_breaks())+
    labs( x = "Genotype",y="Value")+
    
    scale_fill_manual( values = c("darkgray", "burlywood2", "darkolivegreen","black"))
}

# pull unique combos, remove strain column(don't need specific strain info at this point)
#processed_mapping_df <- distinct(select(processed_mapping_df, -strain,-allele,-value))
processed_mapping_df<- processed_mapping_df %>% distinct(trait,marker,strain)


#create family and method columns
processed_mapping_df$family <- stringr::str_split_fixed(processed_mapping_df$trait, "_TRANS_",2)[,2]
processed_mapping_df$method <- stringr::str_split_fixed(processed_mapping_df$trait, "_TRANS_",2)[,1]

#read in position data and create family column
positions <- read.table("CtCp_all_nonredundant.txt")
names(positions)<-c("CHROM","start","end","TE","orientation","method","strain","class")
positions$family<- stringr::str_split_fixed(positions$TE, regex("_(non-)?reference"),2)[,1]
positions$family<- paste(stringr::str_split_fixed(positions$family, "_",4)[,3],stringr::str_split_fixed(positions$family, "_",4)[,4],sep="_")
positions$family <- gsub("_$" ,"",positions$family)
positions$family <- gsub("_non-reference(.*)$" ,"",positions$family)

#select traits above BF.....this step not needed, double checking everything is above BF
selection<-filter(processed_mapping_df, log10p > BF)

#extract the count base traits
base_traits <-selection[(selection$method=="absent"| selection$method=="new" |selection$method=="reference"|selection$method=="ZERO_new"|selection$method=="ONE_new"), ]
counts<-subset(base_traits, grepl("_C$", base_traits$family))
counts$family <- gsub("_C$" ,"",counts$family)

processed_mapping_df <- distinct(select(processed_mapping_df, -strain,-allele,-value))
processed_mapping_df<- processed_mapping_df %>% distinct(trait,marker)

#pull out only position traits from mappings dataframe
position_traits<-subset(selection,
                        grepl('^I', selection$trait) |
                          grepl('^V', selection$trait) |
                          grepl('^X', selection$trait))

#create family column
position_traits$family  <- paste(stringr::str_split_fixed(position_traits$trait, "_",4)[,3],stringr::str_split_fixed(position_traits$trait, "_",4)[,4],sep="_")
position_traits$family <- gsub("_$" ,"",position_traits$family)
position_traits$family <- gsub("_non-reference(.*)$" ,"",position_traits$family)

# add position TRAIT_col family info to processed_mapping_df
processed_mapping_df<-processed_mapping_df %>%mutate(family = ifelse(processed_mapping_df$trait %in% position_traits$trait, (paste(stringr::str_split_fixed(processed_mapping_df$trait, "_",4)[,3],stringr::str_split_fixed(processed_mapping_df$trait, "_",4)[,4],sep="_")), processed_mapping_df$family))

selection<-counts

#strip count marker and remnant marks from dataframes
selection$trait <- gsub("_C$" ,"",selection$trait)
hm$trait <- gsub("_C$" ,"",hm$trait)
processed_mapping_df$trait <- gsub("_C$" ,"",processed_mapping_df$trait)
processed_mapping_df$family <- gsub("_C$" ,"",processed_mapping_df$family)
processed_mapping_df$family <- gsub("_$" ,"",processed_mapping_df$family)
processed_mapping_df$family <- gsub("_non-reference(.*)$" ,"",processed_mapping_df$family)

processed_mapping_df<-mutate(processed_mapping_df,ID=paste(trait,peak_id,sep="_"))
copy<-processed_mapping_df

##FIX HERE
processed_mapping_df<-mutate(processed_mapping_df,SNP_col=ifelse(is.na(peak_id), "NO","PASS"))
#### NEED TO FIX HERE

count_QTL<-mutate(count_QTL, trait2=gsub("_\\d+$","",trait)) 
processed_mapping_df<-filter(processed_mapping_df,CHROM != "MtDNA")
selection<-filter(selection,grepl('total',family))

#HERE
class_subset<- positions %>% distinct(class,family) %>% select(class,family)

selection<-arrange(selection,family,method)
label<-expression(bold(-log["10"](p)))
selection<-filter(selection, trait=="new_TRANS_total_dnatransposon")
for (i in unique(selection$trait)){
  specific_trait<- processed_mapping_df[processed_mapping_df$trait == i, ]
  empty <-specific_trait[specific_trait$method==NA,]
  #specific_trait_mx <- max(specific_trait$log10p)
  class_TE<-unique(filter(selection,trait==i)$class)
  pvalues<-filter(specific_trait,log10p !="Inf") #
  specific_trait_mx <- max(pvalues$log10p) #
  TE<-specific_trait$family[1]
  rect_data<-filter(specific_trait,SNP_col==ifelse(is.na(peak_id), "NO", "PASS"))
  plot_method<-unique(filter(selection,trait==i)$method)
  plot_title<-gsub(".*_TRANS_","",i)
  plot_title<-gsub("_CE$","",plot_title)
  plot_title<-gsub("WBTransposon","WBT",plot_title)
  plot_title<-gsub("total","Total",plot_title)
  plot_title<-gsub("Total$","Total Transposons",plot_title)
  plot_title<-gsub("_"," ",plot_title)
  plot_title<-gsub("retrotransposon","Retrotransposons",plot_title)
  plot_title<-gsub("dnatransposon","DNA Transposons",plot_title)
  plot_title<-paste(plot_title,ifelse(plot_method=="ZERO_new","(ins)",ifelse(plot_method=="absent","(abs)",ifelse(plot_method=="ONE_new", "(ins)",ifelse(plot_method=="new", "(ins)", "(ref)")))),sep=" ")
  plot_title
  
  ##check for NAs
  #sapply(Mappings, function(x)all(is.na(x)))
  A<- processed_mapping_df %>%
    filter(trait == i)%>%
    .[order(.$peak_id,na.last=FALSE),]%>% 
    ggplot(.)+
    aes(x=POS/1e6,y=log10p)+
    geom_rect(data=rect_data,mapping=aes(xmin=startPOS/1e6, xmax=endPOS/1e6, ymin=0, ymax= Inf),fill="thistle1", alpha=1)+
    geom_point(aes( color=ifelse(log10p> BF & SNP_col=="PASS", 'red', 'black')),size=1)+
    
    facet_grid(.~CHROM,scale="free_x",space = "free_x")+scale_color_identity()  +
    ggtitle(plot_title)+
    geom_hline(aes(yintercept=BF),color="grey60",linetype="dashed")+
    theme(strip.background = element_rect(fill = "white"),
          strip.text.x = element_text(size =18, colour = "black",face="bold"),
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(color="black", size=0.5, linetype="solid", fill=NA),
          panel.margin = unit(.6, "lines"),
          plot.title = element_text(size = 18,face="bold"),
          axis.ticks =element_line(colour = "black"),
          axis.text.y = element_text(size=18,colour = "black"),
          axis.text.x = element_text(size=18,colour = "black"),
          axis.title=element_text(size=18,face="bold"),
          plot.margin=unit(c(.1,.1,.1,.1), "cm"),
          
          #ADDD COLOR HERE STILLLL
          # plot.title = element_text(colour=ifelse(class_TE=="dnatransposon","navy",ifelse(class_TE=="retrotransposon","brown3","darkgoldenrod2"))),
          legend.position=('none'))+
    labs(x="Chromosome Position (Mb)",y=label)+
    scale_y_continuous(expand=c(0,0),limits=c(0,specific_trait_mx+.075*specific_trait_mx),labels = function(x) format(x,width = 4),breaks= pretty_breaks())

  plot(A)
  setwd("/Users/kristen/Documents/transposon_figure_data/figures")
  ggsave(filename="fige_total_dna_manh.tiff",dpi=300, width=10,height=4,units="in")
  
  g<-ggplotGrob(A)

  panels <- g$layout$t[grep("panel", g$layout$name)]
  g$heights[panels] <- lapply(c(8,8), unit, "null")
  
  
  
  plist<-vector()
  for (p in rect_data$peak_id){
    plist<-c(plist,p)
  }
  plist<-sort(plist)
  box_list <- lapply(c(plist),FUN=function(x){gwasPxG(i,x)})
  (do.call("grid.arrange", c(box_list, ncol=length(plist))))

}

setwd("/Users/kristen/Documents/transposon_figure_data/figures")
final<-grid.draw(grid.arrange(g,do.call("grid.arrange", c(box_list, ncol=length(plist))),heights=c(.70,.30)))
ggsave(filename="fige_total_dna.tiff",dpi=300, width=7.5,height=3.5,units="in")

d<-grid.arrange(g,do.call("grid.arrange", c(box_list, ncol=length(plist))),heights=c(.70,.30))
ggsave(d,filename="fige_total_dna_manh_box.tiff",dpi=300, width=10,height=6,units="in")


