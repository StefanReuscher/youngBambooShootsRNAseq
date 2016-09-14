#PCA-related figures from Gamuyao et al. 2016

library(edgeR)
library(ggplot2)

#### prepare Figure 4A: PCA on glm_corrected CPM, no replicates ####
# load glm object
load("./data_in/y_glmFit")

# extract fitted values
glmFitted <- y_glmFit$fitted.values

#average over replicates
mean_glmFitted <- cbind("SAM" = rowMeans(glmFitted[,1:4]),
                        "YIN" = rowMeans(glmFitted[,5:8]),
                        "MIN" = rowMeans(glmFitted[,9:12]),
                        "YNO" = rowMeans(glmFitted[,13:16]),
                        "MNO" = rowMeans(glmFitted[,17:20]))

#remove rows which are all 0
mean_glmFitted <- mean_glmFitted[rowSums(mean_glmFitted) != 0,]

#do PCA
PCA_data <- prcomp(t(mean_glmFitted), scale = TRUE,  scores = TRUE)

#look at scores and % variation explained
summary(PCA_data)

#extract the scores and plot
scores <- as.data.frame(PCA_data$x)
ggplot(data = scores, aes(x = PC1, y = PC2, label = rownames(scores)))+
  geom_hline(yintercept = 0) + #new axis
  geom_vline(xintercept = 0) +
  geom_text(colour = "black", alpha = 1)+
  geom_point(color = "black", shape = 18, size = 3)+
  xlab(paste("PC1", round((PCA_data$sdev^2/sum(PCA_data$sdev^2)*100)[1], 1 ), "%")) + #put % of variation to axis
  ylab(paste("PC2", round((PCA_data$sdev^2/sum(PCA_data$sdev^2)*100)[2], 1 ), "%"))+
  theme(panel.grid.major = element_line(colour="grey75"),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(vjust = 0),
        axis.title.y = element_text(vjust = 0.25),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black", face = "bold"))

#save the plot for publication
ggsave("./plots/Fig4B_PCAnoRepsGLM.svg", 
       units = "cm", width = 12, height = 10, dpi = 300)


#### prepare Figure S1: CPM and GLMcpm ####

#similar analyses are performed for Figures S1A and B

#S1A
load("./data_in/allCountsRefbased")
CPM <- cpm(all_counts)
CPM <- CPM[rowSums(CPM) != 0,]
glmFitted <- glmFitted[rowSums(glmFitted) != 0,]

PCA_data <- prcomp(t(CPM), scale = TRUE,  scores = TRUE)
summary(PCA_data)
scores <- as.data.frame(PCA_data$x)
ggplot(data = scores, aes(x = PC1, y = PC2, label = rownames(scores)))+
  geom_hline(yintercept = 0) + #new axis
  geom_vline(xintercept = 0) +
  geom_text(colour = "black", alpha = 1, size = 2, nudge_y = -8)+
  geom_point(color = "black", shape = 18, size = 3)+
  xlab(paste("PC1", round((PCA_data$sdev^2/sum(PCA_data$sdev^2)*100)[1], 1 ), "%")) + #put % of variation to axis
  ylab(paste("PC2", round((PCA_data$sdev^2/sum(PCA_data$sdev^2)*100)[2], 1 ), "%"))+
  theme(panel.grid.major = element_line(colour="grey75"),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(vjust = 0),
        axis.title.y = element_text(vjust = 0.25),
        axis.text = element_text(size = 8, colour = "black"),
        axis.title = element_text(size = 8, colour = "black", face = "bold"))
ggsave("./plots/SupplFig1A_PCAcpm.svg", 
       units = "cm", width = 15, height = 10, dpi = 300)

#S1B
PCA_data <- prcomp(t(glmFitted), scale = TRUE,  scores = TRUE)
summary(PCA_data)
scores <- as.data.frame(PCA_data$x)
ggplot(data = scores, aes(x = PC1, y = PC2, label = rownames(scores)))+
  geom_hline(yintercept = 0) + #new axis
  geom_vline(xintercept = 0) +
  geom_text(colour = "black", alpha = 1, size = 2, nudge_y = -8)+
  geom_point(color = "black", shape = 18, size = 3)+
  xlab(paste("PC1", round((PCA_data$sdev^2/sum(PCA_data$sdev^2)*100)[1], 1 ), "%")) + #put % of variation to axis
  ylab(paste("PC2", round((PCA_data$sdev^2/sum(PCA_data$sdev^2)*100)[2], 1 ), "%"))+
  theme(panel.grid.major = element_line(colour="grey75"),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(vjust = 0),
        axis.title.y = element_text(vjust = 0.25),
        axis.text = element_text(size = 8, colour = "black"),
        axis.title = element_text(size = 8, colour = "black", face = "bold"))
ggsave("./plots/SupplFig1B_PCAglmFit.svg", 
       units = "cm", width = 15, height = 10, dpi = 300)
