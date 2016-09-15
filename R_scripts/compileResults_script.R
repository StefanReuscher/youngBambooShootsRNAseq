# compile relevant information into one results object
# this was used as the basis for Table S2 
# need objects created before and will not run on its own

results <- cbind(all_counts,
                 
                 CPM,
                 
                 "SAM_mean" = rowMeans(CPM[,1:4]),
                 "SAM_SD"   = rowSds(CPM[,  1:4]),
                 "YIN_mean" = rowMeans(CPM[,5:8]),
                 "YIN_SD"   = rowSds(CPM[,  5:8]),
                 "MIN_mean" = rowMeans(CPM[,9:12]),
                 "MIN_SD"   = rowSds(CPM[,  9:12]),
                 "YNO_mean" = rowMeans(CPM[,13:16]),
                 "YNO_SD"   = rowSds(CPM[,  13:16]),
                 "MNO_mean" = rowMeans(CPM[,17:20]),
                 "MNO_SD"   = rowSds(CPM[,  17:20]),
                 
                 "SAM_YIN_FC"    = SAM_YIN_results$logFC,
                 "SAM_YIN_FDR_P" = SAM_YIN_results$FDR,
                 "SAM_MIN_FC"    = SAM_MIN_results$logFC,
                 "SAM_MIN_FDR_P" = SAM_MIN_results$FDR,
                 "SAM_YNO_FC"    = SAM_YNO_results$logFC,
                 "SAM_YNO_FDR_P" = SAM_YNO_results$FDR,
                 "SAM_MNO_FC"    = SAM_MNO_results$logFC,
                 "SAM_MNO_FDR_P" = SAM_MNO_results$FDR,
                 
                 "YIN_MIN_FC"    = YIN_MIN_results$logFC,
                 "YIN_MIN_FDR_P" = YIN_MIN_results$FDR,
                 "YIN_YNO_FC"    = YIN_YNO_results$logFC,
                 "YIN_YNO_FDR_P" = YIN_YNO_results$FDR,
                 "YIN_MNO_FC"    = YIN_MNO_results$logFC,
                 "YIN_MNO_FDR_P" = YIN_MNO_results$FDR,
                 
                 "MIN_YNO_FC"    = MIN_YNO_results$logFC,
                 "MIN_YNO_FDR_P" = MIN_YNO_results$FDR,
                 "MIN_MNO_FC"    = MIN_MNO_results$logFC,
                 "MIN_MNO_FDR_P" = MIN_MNO_results$FDR,
                 
                 "YNO_MNO_FC"    = YNO_MNO_results$logFC,
                 "YNO_MNO_FDR_P" = YNO_MNO_results$FDR           
)

# add descriptions to the results object
CDS <- readBStringSet("./data_in/P_heterocycla_v1.0.genemodel-cds-DNA.fa")
CDS_names <- ldply(strsplit(names(CDS), split = "|", fixed = TRUE))
annot <- data.frame("ID" = CDS_names$V1,
                    "fromContig" = substr(CDS_names$V1, 1, 10),
                    "strand" = CDS_names$V7,
                    "start" = CDS_names$V4,
                    "stop" = CDS_names$V5,
                    "Description" = substr(CDS_names$V9,
                                           start = 33, stop = 1000),
                    stringsAsFactors = FALSE
)

results <- cbind(results,
                 annot[annot$ID %in% rownames(results),])

#write the results object, not included, see Table S2
#write.csv(results, file = "./data/bamboo_results_DEG_GLM.txt")

