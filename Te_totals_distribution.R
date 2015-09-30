#!/usr/bin/R
# this script
# 1) plots no. transposon events vs each other for each possible pairing
# 2) plots total transposons vs strain per insertions, references, and absences
# USE: Te_totals_distribution.R

library(ggplot2)
setwd("/Users/kristen/Documents/transposon_figure_data")
summarydata <- read.table("FINAL_RESULTS.txt",header=TRUE)

##SCATTER
nrow(summarydata)
names(summarydata)

#reformat the data
total_absence<-summarydata[summarydata$method=="absent",2:5]
total_reference<-summarydata[summarydata$method=="reference",2:5]
total_insertion=summarydata[summarydata$method=="new",2:5]
total_absence
total_reference
total_insertion
names(total_insertion)
initial_merge <- merge(total_absence,total_reference, by="sample")
final_merge <- merge(initial_merge,total_insertion, by="sample")
names(final_merge)
final_merge

names(final_merge)<- c("sample", "method.x", "end_tes.x", "total_absences", "method.y", "end_tes.y", "total_references", "method", "end_tes", "total_insertions")

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



#1 ABSENCE vs INSERTION
fit <- lm(total_absences ~ total_insertions, data = final_merge)
print(summary(fit))
nums<-summary(fit)
ad_r2<-round(nums$adj.r.squared,digits=3) #adjusted R2
print(ad_r2)

max_insertions<-max(final_merge$total_insertions)
max_absences<-max(final_merge$total_absences)

la <- paste("italic(r)^2 == ", ad_r2)
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
ggsave(filename="absence_vs_insertion.tiff",
       dpi=300,
       width=4,
       height=4,
       units="in")

#2 ABSENCE vs REFERENCE
fit <- lm(total_absences ~ total_references, data = final_merge)
print(lm(total_absences ~ total_references, data = final_merge))
print(summary(fit))
print(fit)
nums<-summary(fit)
ad_r2<-round(nums$adj.r.squared,digits=3) #adjusted R2
print(ad_r2)

max_references<-max(final_merge$total_references)
la <- paste("italic(r)^2 == ", ad_r2)
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

ggsave(filename="absence_vs_reference.tiff",
       dpi=300,
       width=4,
       height=4,
       units="in")

#3 INSERTION vs REFERENCE
fit <- lm(total_insertions ~ total_references, data = final_merge)
nums<-summary(fit)
ad_r2<-round(nums$adj.r.squared,digits=3) #adjusted R2
print(ad_r2)

la <- paste("italic(r)^2 == ", ad_r2)
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

ggsave(filename="insertion_vs_reference.tiff",
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
m <- ggplot(insertions, aes(x=reorder(insertions$sample,insertions$total_tes), y=insertions$total_tes)) 
m<- m + geom_point(size=.75) +aes(group=1)+
   theme(axis.text.x = element_text(angle = 90, hjust = 1, color="black",size=5),
         axis.text.y = element_text(angle = 90, hjust = 1, color="black",size=9),
         axis.ticks =element_line(colour = "black"))+
  labs(x="Strain", y="Number of Insertion Events")
m
ggsave(filename="Insertions_per_Strain.tiff",
       dpi=300,
       width=7.5,
       height=3.5,
       units="in")

#ABSENCES
absences<-summarydata[summarydata$method=="absent",]
absences<-(absences[ order(absences$total_tes), ])
#plot(absences$total_tes~absences$sample)
#pdf(file = "absences_per_strain.pdf")
m <- ggplot(absences, aes(x=reorder(absences$sample,absences$total_tes), y=absences$total_tes)) 
m<- m + geom_point(size=.75) +aes(group=1)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, color="black",size=5),
        axis.text.y = element_text(angle = 90, hjust = 1, color="black",size=9),
        axis.ticks =element_line(colour = "black"))+
  labs(x="Strain", y="Number of Absence Events")
m
ggsave(filename="Absences_per_Strain.tiff",
       dpi=300,
       width=7.5,
       height=3.5,
       units="in")

#REFERENCE
references<-summarydata[summarydata$method=="reference",]
references<-(references[ order(references$total_tes), ])
#plot(references$total_tes~references$sample)
#pdf(file = "references_per_strain.pdf")
m <- ggplot(references, aes(x=reorder(references$sample,references$total_tes), y=references$total_tes)) 
m<- m + geom_point(size=.75) +aes(group=1)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, color="black",size=5),
        axis.text.y = element_text(angle = 90, hjust = 1, color="black",size=9),
        axis.ticks =element_line(colour = "black"))+
  labs(x="Strain", y="Number of Reference Events")
m
ggsave(filename="References_per_Strain.tiff",
       dpi=300,
       width=7.5,
       height=3.5,
       units="in")
