library(dplyr)
library(ggplot2)
library(data.table)
library(grid)
library(stringr)
library(gridExtra)
setwd("/Users/kristen/Documents/transposon_figure_data/data")
load('SignificantMappings_Results_Activity.Rda')


Mappings$family <- stringr::str_split_fixed(Mappings$traits.i., "_TRANS_",2)[,2]
Mappings$method <- stringr::str_split_fixed(Mappings$traits.i., "_TRANS_",2)[,1]


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


positions <- read.table("CtCp_all_nonredundant.txt",header=TRUE)
names(positions)<-c("chr","start","end","TE","orientation","method","strain","class")
positions$family<- stringr::str_split_fixed(positions$TE, regex("_(non-)?reference"),2)[,1]
positions$family<- paste(stringr::str_split_fixed(positions$family, "_",4)[,3],stringr::str_split_fixed(positions$family, "_",4)[,4],sep="_")
positions$family <- gsub("_$" ,"",positions$family)
positions$family <- gsub("_non-reference(.*)$" ,"",positions$family)

selection<-filter(Mappings, -log10(ps) > -log10(.05/8000))
base_traits <-selection[(selection$method=="absent"| selection$method=="new" |selection$method=="reference"|selection$method=="ZERO_new"|selection$method=="ONE_new"), ]
counts<-subset(base_traits, grepl("_C$", selection$family))
counts$family <- gsub("_C$" ,"",counts$family)
#pull out only position traits from mappings dataframe
position_traits<-subset(selection,
                        grepl('^I', selection$traits.i.) |
                          grepl('^V', selection$traits.i.) |
                          grepl('^X', selection$traits.i.))

selection<-rbind(counts,position_traits)
selection$traits.i. <- gsub("_C$" ,"",selection$traits.i.)

selection<-selection[1:5,] #UNCOMMENT LATER

for (i in unique(selection$traits.i.)){
  specific_trait<- Mappings[Mappings$traits.i == i, ]
  empty <-specific_trait[specific_trait$method==NA,]
  specific_trait_mx <- max(-log10(specific_trait$ps))
  
  ##check for NAs
  #sapply(Mappings, function(x)all(is.na(x)))
  A<- Mappings %>%
    filter(traits.i. == i)%>%
    ggplot(.)+
    aes(x=pos/1e6,y=-log10(ps))+
    geom_point(aes( color=ifelse(-log10(ps)> -log10(.05/8000), 'red', 'black')),size=1)+
    facet_grid(.~chr,scale="free_x",space = "free_x")+scale_color_identity()+
    ggtitle(i)+
    geom_hline(aes(yintercept=-log10(.05/8000)),color="grey60",linetype="dashed")+
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
    labs(x="",y="-log10(p)") +
    
    scale_y_continuous(expand=c(0,0),limits=c(0,specific_trait_mx+.075*specific_trait_mx),labels = function(x) format(x,width = 4))
  
  
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
  traitPositions<-positions[positions$trait==i,]
  
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
  
  m <- ggplot(traitPositions, aes(x=start/1e6))
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
          legend.text=element_text(size=8))
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


