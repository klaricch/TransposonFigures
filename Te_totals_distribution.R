#!/usr/bin/R
# this script
# 1) plots no. transposon events vs each other for each possible pairing
# 2) plots total transposons vs strain per insertions, references, and absences
# USE: Te_totals_distribution.R

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
setwd("/Users/kristen/Documents/transposon_figure_data/data")
summarydata <- read.table("T_kin_C_matrix_full.txt",header=TRUE)
#remove ZERO_new traits
summarydata<-subset(summarydata,!grepl('^ZERO_new', summarydata$trait))
summarydata<-subset(summarydata,!grepl('^coverage', summarydata$trait))
#clean trait names
summarydata$trait <- gsub("_C$" ,"",summarydata$trait)
summarydata$trait <- gsub("^ONE_new" ,"new",summarydata$trait)

#new column that specifies what caller was used
summarydata$method<- stringr::str_split_fixed(summarydata$trait, "_TRANS_",2)[,1]
#new column that specifies TE family
summarydata$transposon<- stringr::str_split_fixed(summarydata$trait, "_TRANS_",2)[,2]
summarydata<-filter(summarydata,transposon=="total")

#names(summarydata)
summarydata<-gather(summarydata, "sample","value",2:(ncol(summarydata)-2))
summarydata<-rename(summarydata,total_tes=value)

#reformat the data
total_absence<-filter(summarydata,method=="absent")
total_reference<-filter(summarydata,method=="reference")
total_insertion<-filter(summarydata,method=="new")




#
#
#
#
#remove coverage, etc, traits
#
#
#
#

#SCATTER
initial_merge <- merge(total_absence,total_reference, by="sample")
final_merge <- merge(initial_merge,total_insertion, by="sample")

names(final_merge)<-c("sample", "trait.x",	"method.x",	"transposon.x",	"total_absences",	"trait.y",	"method.y",	"transposon.y",	"total_references",	"trait",	"method",	"transposon",	"total_insertions")

method_names <- list(
  'absent'="Absence",
  'new'="Insertion",
  'reference'="Reference"
)

method_labeller <- function(variable,value){
  return(method_names[value])
}






#TESTS
cor.test(final_merge$total_absences, final_merge$total_insertions,method="spearman",exact=FALSE)
cor.test(final_merge$total_absences, final_merge$total_references,method="spearman",exact=FALSE)
cor.test(final_merge$total_insertions, final_merge$total_references,method="spearman",exact=FALSE)

cor.test(final_merge$total_absences, final_merge$total_insertions,method="pearson",exact=FALSE)
cor.test(final_merge$total_absences, final_merge$total_references,method="pearson",exact=FALSE)
cor.test(final_merge$total_insertions, final_merge$total_references,method="pearson",exact=FALSE)

fit <- lm(total_absences ~ total_insertions, data = final_merge)
print(summary(fit))
fit <- lm(total_absences ~ total_references, data = final_merge)
print(summary(fit))
fit <- lm(total_insertions ~ total_references, data = final_merge)
print(summary(fit))

cor(final_merge$total_absences, final_merge$total_insertions)
cor(final_merge$total_absences, final_merge$total_references)
cor(final_merge$total_insertions, final_merge$total_references)




######
# add te class info to summarydata(new_TRANS_end_tes will be removed)
setwd("/Users/kristen/Documents/transposon_figure_data/data")
classdata<- read.table("CtCp_all_nonredundant.txt",header=TRUE)
names(classdata)<-c("chr","start","end","TE","orientation","method","strain","class")
classdata$id<- stringr::str_split_fixed(classdata$TE, regex("_(non-)?reference"),2)[,1]
classdata$id<- stringr::str_split_fixed(classdata$TE, regex("_(non-)?reference"),2)[,1]
classdata$id<- paste(stringr::str_split_fixed(classdata$id, "_",4)[,3],stringr::str_split_fixed(classdata$id, "_",4)[,4],sep="_")
classdata$id <- gsub("_$" ,"",classdata$id)
classdata$id <- gsub("_non-reference(.*)$" ,"",classdata$id)
classdata<-mutate(classdata, trait=paste(method,"TRANS",id,sep="_"))
class_subset <- classdata %>% distinct(id) %>% select(id,class)
print(class_subset)
df<-t(class_subset)
print(df)
names(summarydata)
names(df)
colnames(df) <- df[1,]
print(df[1,])
names(df)
newdf <-bind_rows(summarydata, class_subset)
print(newdf)
names(newdf)

####

#1 ABSENCE vs INSERTION
#spearman correlation
correlation<-cor.test(final_merge$total_absences, final_merge$total_insertions,method="spearman",exact=FALSE)
rho<-round(correlation$estimate,3)

max_insertions<-max(final_merge$total_insertions)
max_absences<-max(final_merge$total_absences)


la <- paste("italic(rho) == ", rho)
m <- ggplot(final_merge, aes(x=total_insertions, y=total_absences))
m <- m + geom_point(size=1.25) + xlim(0,max_insertions)+ ylim(0,max_insertions)+
  geom_smooth(method="lm",se=FALSE,col="red")+
  geom_abline(slope=1,linetype="dashed",colour="gray52")+
  annotate("text", x=.1*max_insertions, y=.9*max_insertions,label=la,parse=TRUE, colour="red",size=2.5)+
  theme(strip.text.x = element_text(size = 9, colour = "black"),
        strip.background = element_blank(),
        legend.position=c(.90,0.75),
        legend.background = element_rect(fill=FALSE),
        legend.text=element_text(size=9),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=9),
        axis.text.x = element_text(colour = "black",size=9),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9))+
  scale_fill_manual(values = c("darkorange", "turquoise3", "slateblue1")) +
  guides(fill=FALSE) +
  labs(x = "Insertion Events", y = "Absence Events")
m
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="Absence_vs_Insertion.tiff",
       dpi=300,
       width=4,
       height=4,
       units="in")

#2 ABSENCE vs REFERENCE
#spearman correlation
correlation<-cor.test(final_merge$total_absences, final_merge$total_references,method="spearman",exact=FALSE)
rho<-round(correlation$estimate,3)

max_references<-max(final_merge$total_references)
max_absences<-max(final_merge$total_absences)

la <- paste("italic(rho) == ", rho)
m <- ggplot(final_merge, aes(x=total_references, y=total_absences))
m + geom_point(size=1.25) + xlim(0,max_references)+ ylim(0,max_references)+
  geom_smooth(method="lm",se=FALSE,col="red")+
  geom_abline(slope=1,linetype="dashed",colour="gray52")+
  annotate("text", x=.1*max_references, y=.9*max_references,label=la,parse=TRUE, colour="red",size=2.5)+
  theme(strip.text.x = element_text(size = 9, colour = "black"),
        strip.background = element_blank(),
        legend.position=c(.90,0.75),
        legend.background = element_rect(fill=FALSE),
        legend.text=element_text(size=9),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=9),
        axis.text.x = element_text(colour = "black",size=9),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9))+
  scale_fill_manual(values = c("darkorange", "turquoise3", "slateblue1")) +
  guides(fill=FALSE) +
  labs(x = "Reference Events", y = "Absence Events")

ggsave(filename="Absence_vs_Reference.tiff",
       dpi=300,
       width=4,
       height=4,
       units="in")

#3 INSERTION vs REFERENCE
#spearman correlation
correlation<-cor.test(final_merge$total_insertions, final_merge$total_references,method="spearman",exact=FALSE)
rho<-round(correlation$estimate,3)

max_references<-max(final_merge$total_references)
max_insertions<-max(final_merge$total_insertions)
la <- paste("italic(rho) == ", rho)

max_references<-max(final_merge$total_references)
m <- ggplot(final_merge, aes(x=total_references, y=total_insertions))
m + geom_point(size=1.25) + xlim(0,max_references)+ ylim(0,max_references)+
  geom_smooth(method="lm",se=FALSE,col="red")+
  geom_abline(slope=1,linetype="dashed",colour="gray52")+
  annotate("text", x=.1*max_references, y=.9*max_references,label=la,parse=TRUE, colour="red",size=2.5)+
  theme(strip.text.x = element_text(size = 9, colour = "black"),
        strip.background = element_blank(),
        legend.position=c(.90,0.75),
        legend.background = element_rect(fill=FALSE),
        legend.text=element_text(size=9),
        panel.background = element_rect(fill = "white"),
        axis.ticks =element_line(colour = "black"),
        axis.text.y = element_text(colour = "black",size=9),
        axis.text.x = element_text(colour = "black",size=9),
        axis.line=element_line(linetype="solid"),
        axis.title=element_text(size=9))+
  scale_fill_manual(values = c("darkorange", "turquoise3", "slateblue1")) +
  guides(fill=FALSE) +
  labs(x = "Reference Events", y = "Insertion Events")

ggsave(filename="Insertion_vs_Reference.tiff",
       dpi=300,
       width=4,
       height=4,
       units="in")

#TRANSPOSONS vs STRAINS
names(summarydata)
#INSERTIONS
insertions<-summarydata[summarydata$method=="new",]
insertions<-(insertions[ order(insertions$total_tes), ])
#plot(insertions$total_tes~insertions$sample)
#pdf(file = "insertions_per_strain.pdf")
m <- ggplot(insertions, aes(y=reorder(insertions$sample,insertions$total_tes), x=insertions$total_tes)) 
m<- m + geom_point(size=.75) +aes(group=1)+
  theme(axis.text.x = element_text(color="black",size=8),
        axis.text.y = element_text(color="black",size=8),
        axis.title = element_text(color="black",size=9),
        axis.ticks =element_line(colour = "black"))+
  labs(x="Strain", y="Number of Insertion Events")
m
ggsave(filename="Insertions_per_Strain.tiff",
       dpi=300,
       width=7.5,
       height=10,
       units="in")

#ABSENCES
absences<-summarydata[summarydata$method=="absent",]
absences<-(absences[ order(absences$total_tes), ])
#plot(absences$total_tes~absences$sample)
#pdf(file = "absences_per_strain.pdf")
m <- ggplot(absences, aes(y=reorder(absences$sample,absences$total_tes), x=absences$total_tes)) 
m<- m + geom_point(size=.75) +aes(group=1)+
  theme(axis.text.x = element_text(color="black",size=8),
        axis.text.y = element_text(color="black",size=8),
        axis.title = element_text(color="black",size=9),
        axis.ticks =element_line(colour = "black"))+
  labs(x="Strain", y="Number of Absence Events")
m
ggsave(filename="Absences_per_Strain.tiff",
       dpi=300,
       width=7.5,
       height=10,
       units="in")

#REFERENCE
references<-summarydata[summarydata$method=="reference",]
references<-(references[ order(references$total_tes), ])
#plot(references$total_tes~references$sample)
#pdf(file = "references_per_strain.pdf")
m <- ggplot(references, aes(y=reorder(references$sample,references$total_tes), x=references$total_tes)) 
m<- m + geom_point(size=.75) +aes(group=1)+
  theme(axis.text.x = element_text(color="black",size=8),
        axis.text.y = element_text(color="black",size=8),
        axis.title = element_text(color="black",size=9),
        axis.ticks =element_line(colour = "black"))+
  labs(x="Strain", y="Number of Reference Events")
m
ggsave(filename="References_per_Strain.tiff",
       dpi=300,
       width=7.5,
       height=10,
       units="in")
