# read in the raw count data and apply a GLM

library(edgeR)

# load raw counts from tophat2/RSubread
load("./data/allCountsRefbased")

#define the structure of the experiment
replicates <- c(rep(c("A","B","C","D"),5))
tissues <- c(rep("SAM",4),rep("YIN",4),
             rep("MIN",4),rep("YNO",4),
             rep("MNO",4))

# remove low expressed genes based on empirical threshold
keep <- rowSums(cpm(all_counts) >= 0.5) >= 4 # cpm of at least 0.5 in 4 or more samples
table(keep) #how many to keep
all_counts <- all_counts[keep, ]

#create a DGE list and normfactors
y <- DGEList(all_counts, group = tissues)
y <- calcNormFactors(y)
y$samples #check samples are correct

#make a design matrix using both replicates and tissues as additive factors
design <- model.matrix(~replicates+tissues)
rownames(design) <- colnames(y)

#estimate the average, trended and tagwise dispersion
y <- estimateGLMCommonDisp(y, design, verbose = TRUE)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
plotBCV(y) #diagnostic plot

#fit genewise linear models, using design
y_glmFit <- glmFit(y, design)
save(y_glmFit, file = "./data_out/y_glmFit")

#the y_glmFit object contains all information and is used in downstream analyses
