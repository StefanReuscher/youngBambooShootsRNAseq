
#plot a heatmap showing average expression in each cluster

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

# define a theme for plotting
hmTheme <- theme(axis.text.y = element_text(size = fontsize, color = "black"),
                 axis.title.x = element_blank(),
                 axis.title = element_text(size = fontsize, face = "bold"),
                 axis.ticks = element_blank(),
                 panel.background = element_rect(fill = NA),
                 legend.position = "bottom",
                 legend.direction = "horizontal",
                 legend.title = element_text(size = fontsize),
                 legend.text = element_text(size = fontsize))

#plot a heatmap with k = 16 and nice layout
FC_heatmap_DF <- as.data.frame(myRNASeqData$logFC)
colnames(FC_heatmap_DF) <- c("SAM","YIN","MIN","YNO","MNO")
FC_heatmap_DF$clustID <- clustID[,k] #select k here
FC_heatmap_DF <- aggregate(FC_heatmap_DF,list(FC_heatmap_DF$clustID),mean)
FC_heatmap_DF <- melt(FC_heatmap_DF[,-1],id.vars = "clustID")

ggplot(FC_heatmap_DF, aes(x = variable, y = clustID))+
  geom_tile(aes(fill = value), color = "black")+
  ylab("Cluster ID")+
  scale_y_continuous(breaks = 1:16)+
  scale_fill_gradient2(low = "blue", high = "red",
                       name = expression(atop(log["2"]-fold~change~relative,
                                              to~mean~expression)))+
  hmTheme

ggsave(paste("./plots/kmeansAll.",saveAs, sep = ""),
       units = "cm", width = 16, height = 8.5, dpi = 300)