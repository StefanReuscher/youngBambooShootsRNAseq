#perform k-means clustering of raw count data

library(MBCluster.Seq)
library(ggplot2)
library(reshape2)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#use fitted values from the y_glmfit object
load("./data/y_glmFit")
fitVal <- y_glmFit$fitted.values

####kmeans clustering by MBCluster.Seq####

#experiment structure
GeneID <- row.names(fitVal)
Treatment <- c(1,1,1,1,
               2,2,2,2,
               3,3,3,3,
               4,4,4,4,
               5,5,5,5)
#prepare data for clustering
myRNASeqData <- RNASeq.Data(fitVal,Normalize=NULL,Treatment,GeneID)
# save(myRNASeqData, file = "./data_out/myRNAseqData_bamboo") save intermediate if needed

#the number of clusters k was initially set to this
#later 16 appeared to be the best value
for(no_of_cluster in c(8,10,12,14,16,18,20)){
  
  myc0 <- KmeansPlus.RNASeq(myRNASeqData, nK = no_of_cluster,
                            model = "nbinom",
                            print.steps = TRUE)
  
  myCls <- Cluster.RNASeq(data= myRNASeqData, model = "nbinom",
                          centers = myc0$centers, method = "EM")
  
  myTr <- Hybrid.Tree(data = myRNASeqData,
                      cluster0 = myCls$cluster,
                      model = "nbinom")
  #save obejcts
  assign(paste("myc0_k",no_of_cluster, sep = ""), myc0)
  assign(paste("myCls_k",no_of_cluster,sep = ""), myCls)
  assign(paste("myTr_k",no_of_cluster, sep = ""), myTr)
}
#####

# all gene<->clusterID mapping are brought together
clustID <- cbind( "k8" = myCls_k8$cluster,
                 "k10" = myCls_k10$cluster,
                 "k12" = myCls_k12$cluster,
                 "k14" = myCls_k14$cluster,
                 "k16" = myCls_k16$cluster,
                 "k18" = myCls_k18$cluster,
                 "k20" = myCls_k20$cluster)

#cluster IDs are save and used later for plotting
save(clustID, file = "./data_out/clustID_bamboo")

