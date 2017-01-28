library(ape)
library(ggplot2)
library(phyloseq)

#source("https://bioconductor.org/biocLite.R")
$biocLite("phyloseq")
#install.packages("phyloseq")
setwd("/Users/kristen/Documents/transposon_figure_data/data")

load_tree <- function(treefile) {
  tree <- read.tree(paste0(treefile))
  tree <- root(tree,outgroup = "QX1211", resolve.root = T)
  treeSegs <- tree_layout(phy_tree(tree), ladderize = T)
  treeSegs$edgeDT <- treeSegs$edgeDT  %>% dplyr::mutate(edge.length =
                                                          ifelse(edge.length < 0, 0, edge.length), xright = xleft + edge.length)
  
  label_data <- treeSegs$edgeDT
  
  edgeMap = aes(x = xleft, xend = xright, y = y, yend = y)
  vertMap = aes(x = x, xend = x, y = vmin, yend = vmax)
  labelMap <- aes(x = xright+0.0001, y = y, label = OTU)
  ggplot(data = treeSegs$edgeDT) + geom_segment(edgeMap) +
    geom_segment(vertMap, data = treeSegs$vertDT) +
    geom_text(labelMap, data = label_data, na.rm = TRUE, hjust = -0.05) +
    theme_nothing() +
    scale_color_manual(values= c("#2c7bb6", "#d7191c")) +
    scale_x_continuous(limits = c(
      min(treeSegs$edgeDT$xleft)-0.15,
      max(treeSegs$edgeDT$xright)+0.15
    ),
    expand = c(0,0))
  
}
TE_tree<-load_tree("te.tre")
setwd("/Users/kristen/Documents/transposon_figure_data/figures")
ggsave(filename="TETREE.tiff",
       dpi=300,
       width=6,
       height=20,
       units="in")

ggsave(filename="TETREE.png",
       dpi=300,
       width=6,
       height=20,
       units="in")

