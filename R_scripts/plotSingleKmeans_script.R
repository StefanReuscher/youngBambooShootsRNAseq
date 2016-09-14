#plot a heatmap showing expression in a single cluster

library(ggplot2)
library(reshape2)
library(scales)

#load cluster IDs and expression data
load("./data/myRNAseqData_bamboo")
load("./data/clustID_bamboo")

fontsize <- 8
saveAs <- "svg"
#saveAs <- "png" #choose output format

#choose k (number of clusters)
k <- "k16"

#choose a cluster for plotting
clusterToPlot <- 13

hmTheme <- theme(axis.text.y = element_blank(), 
                 axis.title.y = element_blank(),
                 axis.title = element_text(size = fontsize, face = "bold"),
                 axis.ticks = element_blank(),
                 panel.background = element_rect(fill = NA),
                 legend.position = "bottom",
                 legend.direction = "horizontal",
                 legend.title = element_text(size = fontsize),
                 legend.text = element_text(size = fontsize))

#plot a heatmap from a single cluster
FC_heatmap_DF <- as.data.frame(myRNASeqData$logFC)
colnames(FC_heatmap_DF) <- c("SAM","YIN","MIN","YNO","MNO")
FC_heatmap_DF <- FC_heatmap_DF[clustID[,k] == clusterToPlot,]
FC_heatmap_DF$ID <- row.names(FC_heatmap_DF)

FC_heatmap_DF <- melt(FC_heatmap_DF, id.vars = "ID")

ggplot(FC_heatmap_DF, aes(x = variable, y = ID, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", high = "red",
                       name = expression(atop(log["2"]-fold~change~relative,
                                              to~mean~expression)),
                       limits = c(-3,3))+# #can be modified for better color range
  xlab("tissue")+
  hmTheme

ggsave(paste("./plots/kmeans_",k,"_ID",clusterToPlot,".",saveAs, sep = ""),
       units = "cm", width = 8, height = 8, dpi = 300)
