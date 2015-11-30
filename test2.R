library(dplyr)
library(ggplot2)
library(data.table)
library(grid)
library(stringr)
library(gridExtra)
setwd("/Users/kristen/Documents/transposon_figure_data/data")
load("Processed_Transposon_Mappings.Rda")
load("away_phenos.Rda")

# pull unique combos, remove strain column(don't need specific strain info at this point)
final_processed_mappings <- distinct(select(final_processed_mappings, -strain,-allele,-value))
final_processed_mappings<- final_processed_mappings %>% distinct(pheno,SNPs)

#remove fraction and movement traits
final_processed_mappings<-subset(final_processed_mappings,
                                 grepl('^I', final_processed_mappings$pheno) |
                                   grepl('^V', final_processed_mappings$pheno) |
                                   grepl('^X', final_processed_mappings$pheno)|
                                   grepl('_C$', final_processed_mappings$pheno))

#create family and method columns
final_processed_mappings$family <- stringr::str_split_fixed(final_processed_mappings$pheno, "_TRANS_",2)[,2]
final_processed_mappings$method <- stringr::str_split_fixed(final_processed_mappings$pheno, "_TRANS_",2)[,1]

#set up method labellers
method_names <- list(
  'absent'="absence",
  'new'="insertion",
  'reference'="reference"
)

method_labeller <- function(variable,value){
  if (variable=='method') {
    return(method_names[value])
  }else {
    return(as.character(value))
  }
}


#read in position data and create family column
positions <- read.table("CtCp_all_nonredundant.txt",header=TRUE)
names(positions)<-c("chr","start","end","TE","orientation","method","strain","class")
positions$family<- stringr::str_split_fixed(positions$TE, regex("_(non-)?reference"),2)[,1]
positions$family<- paste(stringr::str_split_fixed(positions$family, "_",4)[,3],stringr::str_split_fixed(positions$family, "_",4)[,4],sep="_")
positions$family <- gsub("_$" ,"",positions$family)
positions$family <- gsub("_non-reference(.*)$" ,"",positions$family)

#select traits above BF.....this step not needed, double checking everything is above BF
selection<-filter(final_processed_mappings, log10p > BF)

#extract the count base traits
base_traits <-selection[(selection$method=="absent"| selection$method=="new" |selection$method=="reference"|selection$method=="ZERO_new"|selection$method=="ONE_new"), ]
counts<-subset(base_traits, grepl("_C$", base_traits$family))
counts$family <- gsub("_C$" ,"",counts$family)

#pull out only position traits from mappings dataframe
position_traits<-subset(selection,
                        grepl('^I', selection$pheno) |
                          grepl('^V', selection$pheno) |
                          grepl('^X', selection$pheno))

#create family column
position_traits$family  <- paste(stringr::str_split_fixed(position_traits$pheno, "_",4)[,3],stringr::str_split_fixed(position_traits$pheno, "_",4)[,4],sep="_")
position_traits$family <- gsub("_$" ,"",position_traits$family)
position_traits$family <- gsub("_non-reference(.*)$" ,"",position_traits$family)

#optional filter for away phenos
position_traits<-position_traits[position_traits$pheno %in% away$pheno,] 

# add position trait family info to final_processed_mappings
final_processed_mappings<-final_processed_mappings %>%mutate(family = ifelse(final_processed_mappings$pheno %in% away$pheno, (paste(stringr::str_split_fixed(final_processed_mappings$pheno, "_",4)[,3],stringr::str_split_fixed(final_processed_mappings$pheno, "_",4)[,4],sep="_")), final_processed_mappings$family))

#bind count and position tratis
selection<-rbind(counts,position_traits)

#strip count marker and remnant marks from dataframes
selection$pheno <- gsub("_C$" ,"",selection$pheno)
final_processed_mappings$pheno <- gsub("_C$" ,"",final_processed_mappings$pheno)
final_processed_mappings$family <- gsub("_C$" ,"",final_processed_mappings$family)
final_processed_mappings$family <- gsub("_$" ,"",final_processed_mappings$family)
final_processed_mappings$family <- gsub("_non-reference(.*)$" ,"",final_processed_mappings$family)

#iterate through the phenotypes and plot the results
for (i in unique(selection$pheno)){
  specific_trait<- final_processed_mappings[final_processed_mappings$pheno == i, ]
  empty <-specific_trait[specific_trait$method==NA,]
  #specific_trait_mx <- max(specific_trait$log10p)
  pvalues<-filter(specific_trait,log10p !="Inf") #
  specific_trait_mx <- max(pvalues$log10p) #
  TE<-specific_trait$family[1]
  ##check for NAs
  #sapply(Mappings, function(x)all(is.na(x)))
  A<- final_processed_mappings %>%
    filter(pheno == i)%>%
    ggplot(.)+
    aes(x=pos/1e6,y=log10p)+
    geom_point(aes( color=ifelse(log10p> BF, 'red', 'black')),size=1)+
    facet_grid(.~chr,scale="free_x",space = "free_x")+scale_color_identity()+
    ggtitle(i)+
    geom_hline(aes(yintercept=BF),color="grey60",linetype="dashed")+
    theme(strip.background = element_rect(fill = "white"),
          strip.text.x = element_text(size = 9, colour = "black",face="bold"),
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(color="black", size=0.5, linetype="solid", fill=NA),
          panel.margin = unit(.6, "lines"),
          axis.ticks =element_line(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black"),
          axis.title=element_text(size=9),
          plot.margin=unit(c(.1,.1,-.25,.1), "cm"),
          legend.position=('none'))+
    labs(x="",y="-log10(p)") #+
  
  #scale_y_continuous(expand=c(0,0),limits=c(0,specific_trait_mx+.075*specific_trait_mx),labels = function(x) format(x,width = 4))
  
  
  # pull out  X maxs of each panel
  panel1<-(ggplot_build(A)$data[[1]])[ggplot_build(A)$data[[1]]$PANEL==1,]
  max1<-(max(panel1$x))
  min1<-(min(panel1$x))
  panel2<-(ggplot_build(A)$data[[1]])[ggplot_build(A)$data[[1]]$PANEL==2,]
  max2<-(max(panel2$x))
  min2<-(min(panel2$x))
  panel3<-(ggplot_build(A)$data[[1]])[ggplot_build(A)$data[[1]]$PANEL==3,]
  max3<-(max(panel3$x))
  min3<-(min(panel3$x))
  panel4<-(ggplot_build(A)$data[[1]])[ggplot_build(A)$data[[1]]$PANEL==4,]
  max4<-(max(panel4$x))
  min4<-(min(panel4$x))
  panel5<-(ggplot_build(A)$data[[1]])[ggplot_build(A)$data[[1]]$PANEL==5,]
  max5<-(max(panel5$x))
  min5<-(min(panel5$x))
  panel6<-(ggplot_build(A)$data[[1]])[ggplot_build(A)$data[[1]]$PANEL==6,]
  max6<-(max(panel6$x))
  min6<-(min(panel6$x))
  
  positions$trait<-paste(positions$method, "TRANS", positions$family, sep="_")
  #traitPositions<-positions[positions$trait==i,]
  traitPositions<-positions[positions$family==TE,]
  
  blank <- data.frame(chr=character(),
                      start=integer(),
                      end=integer(),
                      TE=character(),
                      orientation=character(),
                      method=character(),
                      strain=character(),
                      class=character(),
                      family=character(),
                      trait=character(),
                      stringsAsFactors=FALSE)
  
  blank[1,]<- c( "I", as.integer(10000000),as.integer(10000000),"blank","+","blank","fake","blank","blank","blank")
  blank[2,]<- c( "II", as.integer(10000000),as.integer(10000000),"blank","+","blank","fake","blank","blank","blank")
  blank[3,]<- c( "III", as.integer(10000000),as.integer(10000000),"blank","+","blank","fake","blank","blank","blank")
  blank[4,]<- c( "IV", as.integer(10000000),as.integer(10000000),"blank","+","blank","fake","blank","blank","blank")
  blank[5,]<- c( "V", as.integer(10000000),as.integer(10000000),"blank","+","blank","fake","blank","blank","blank")
  blank[6,]<- c( "X", as.integer(10000000),as.integer(10000000),"blank","+","blank","fake","blank","blank","blank")
  
  traitPositions<-rbind(traitPositions,blank)
  traitPositions$start<-as.integer(traitPositions$start)
  #m <- ggplot(summarydata, aes(x=start/1e6,fill=class))
  #m <-m + geom_bar(binwidth=.25)+
  
  #ggplot(data = combo, aes(x = TEMP_support,y=TELOCATE_support,color=ifelse(method=="absent","darkorange",ifelse(method=="blank","black",ifelse(method=="insertion",""turquoise3"","slateblue1")))))+scale_color_identity()
  m <- ggplot(traitPositions, aes(x=start/1e6,color=ifelse(method=="absent","darkorange",ifelse(method=="blank","black",ifelse(method=="new","turquoise3","slateblue1")))))+scale_color_identity()
  m <-m + geom_bar(data=subset(traitPositions,strain=="fake"), fill="white", colour="white", binwidth=.25)
  m <-m + geom_bar(data=subset(traitPositions,strain!="fake"), binwidth=.25)+
    facet_grid(. ~ chr,scale="free",labeller=method_labeller,drop=FALSE)+
    ggtitle("")+
    geom_point(data = subset(traitPositions, chr=="I"),aes(x=max1,y=0),alpha=0)+
    geom_point(data = subset(traitPositions, chr=="II"),aes(x=max2,y=0),alpha=0) +
    geom_point(data = subset(traitPositions, chr=="III"),aes(x=max3,y=0),alpha=0) +
    geom_point(data = subset(traitPositions, chr=="IV"),aes(x=max4,y=0),alpha=0) +
    geom_point(data = subset(traitPositions, chr=="V"),aes(x=max5,y=0),alpha=0) +
    geom_point(data = subset(traitPositions, chr=="X"),aes(x=max6,y=0),alpha=0) +
    
    geom_point(data = subset(traitPositions, chr=="I"),aes(x=min1,y=0),alpha=0) +
    geom_point(data = subset(traitPositions, chr=="II"),aes(x=min2,y=0),alpha=0) +
    geom_point(data = subset(traitPositions, chr=="III"),aes(x=min3,y=0),alpha=0) +
    geom_point(data = subset(traitPositions, chr=="IV"),aes(x=min4,y=0),alpha=0) +
    geom_point(data = subset(traitPositions, chr=="V"),aes(x=min5,y=0),alpha=0) +
    geom_point(data = subset(traitPositions, chr=="X"),aes(x=min6,y=0),alpha=0)+
    
    labs(x = "Chromosome Position (Mb)", y="Number of Transposition Events")+
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          #strip.text = element_text(size = 9, colour = "black",face="bold"),
          panel.margin = unit(.6, "lines"),
          panel.border = element_rect(fill=NA,colour = "black"),
          panel.background = element_rect(fill = "white"),
          axis.ticks =element_line(colour = "black"),
          axis.title=element_text(size=9),
          axis.text.y = element_text(colour = "black",size=9),
          axis.text.x = element_text(colour = "black",size=9),
          legend.title=element_blank(),
          # legend.position="bottom",
          plot.margin=unit(c(-.25,.1,.1,.1), "cm"),
          legend.position=('none'))
  m
  
  #now can check plot for max value and set y limit to a certain percent above that max value 
  m <- m + scale_y_continuous(expand = c(0,0),limits=c(0,max(ggplot_build(m)$panel$ranges[[1]]$y.range)*1.075)) 
  
  
  library(gtable)
  g1<-ggplotGrob(A)
  g2<-ggplotGrob(m)
  #Bind the tables
  g<-gtable:::rbind_gtable(g1, g2, "first")
  #Remove a row between the plots
  #g <- gtable_add_rows(g, unit(-1,"cm"), pos=nrow(g1))
  #draw
  panels <- g$layout$t[grep("panel", g$layout$name)]
  g$heights[panels] <- lapply(c(40,40), unit, "null")
  grid.newpage()
  grid.draw(g)
  
}




# filter for away traits

test<-position_traits[position_traits$pheno %in% away$pheno,] 
test2<-counts[counts$pheno==NA,]
print(counts$pheno)
