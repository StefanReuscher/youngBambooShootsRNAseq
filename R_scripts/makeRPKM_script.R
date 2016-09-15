#produce rpkm values to direclty display gene expression

require(edgeR)
require(matrixStats)
require(Biostrings)

#load count data and remove low expressed genes
load("./data_in/allCountsRefbased")
keep <- rowSums(cpm(all_counts) >= 0.5) >= 4 # cpm of at least 0.5 in 4 or more samples
table(keep)
all_counts <- all_counts[keep, ]

# get the transcript length from the bamboo reference CDS
# I did not include that file as it is too big. It is hosted
# at bamboogdb.org
CDS <- readBStringSet("./data_in/P_heterocycla_v1.0.genemodel-cds-DNA.fa")
CDS <- CDS[keep]

# make rpkm
RPKM <- rpkm(all_counts, gene.length = width(CDS))

#sum up the data
RPKM_mean <- data.frame("SAM_mean" = rowMeans(RPKM[,1:4]),
                        "YIN_mean" = rowMeans(RPKM[,5:8]),
                        "MIN_mean" = rowMeans(RPKM[,9:12]),
                        "YNO_mean" = rowMeans(RPKM[,13:16]),
                        "MNO_mean" = rowMeans(RPKM[,17:20]))
RPKM_mean <- round(RPKM_mean,2 )
RPKM_sd <- data.frame("SAM_mean" = rowSds(RPKM[,1:4]),
                        "YIN_mean" = rowSds(RPKM[,5:8]),
                        "MIN_mean" = rowSds(RPKM[,9:12]),
                        "YNO_mean" = rowSds(RPKM[,13:16]),
                        "MNO_mean" = rowSds(RPKM[,17:20]))
RPKM_sd <- round(RPKM_sd, 2)
rownames(RPKM_sd) <- rownames(RPKM_mean)

#save as .txt file, those files are in the repo
write.table(RPKM_mean, file = "./data/RPKM_mean.txt",
            sep = "\t", quote = FALSE)
write.table(RPKM_sd, file = "./data/RPKM_sd.txt",
            sep = "\t", quote = FALSE)
