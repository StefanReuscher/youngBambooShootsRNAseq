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
#save(y_glmFit, file = "./data/y_glmFit") save if needed


# make pairwise combinations of tissues
SAM_YIN_LRT <- glmLRT(y_glmFit,  contrast = c(0,0,0,0,0,1,-1,0))
SAM_MIN_LRT <- glmLRT(y_glmFit,  coef = 6) #because MIN is the intercept column and I dont know another synthax
SAM_YNO_LRT <- glmLRT(y_glmFit,  contrast = c(0,0,0,0,0,1,0,-1))
SAM_MNO_LRT <- glmLRT(y_glmFit,  contrast = c(0,0,0,0,-1,1,0,0))

YIN_MIN_LRT <- glmLRT(y_glmFit,  coef = 7)
YIN_YNO_LRT <- glmLRT(y_glmFit,  contrast = c(0,0,0,0,0,0,1,-1))
YIN_MNO_LRT <- glmLRT(y_glmFit,  contrast = c(0,0,0,0,-1,0,1,0))

MIN_YNO_LRT <- glmLRT(y_glmFit,  contrast = c(0,0,0,0,0,0,0,-1))
MIN_MNO_LRT <- glmLRT(y_glmFit,  contrast = c(0,0,0,0,-1,0,0,0))

YNO_MNO_LRT <- glmLRT(y_glmFit,  contrast = c(0,0,0,0,-1,0,0,1))

#save all LRT objects for DEG analysis in the app
save(list = grep("LRT", ls(), value = TRUE), file = "./data/LRT_objects")

#create results tables for each comparison
SAM_YIN_results <- topTags(SAM_YIN_LRT, n = NULL)$table
SAM_MIN_results <- topTags(SAM_MIN_LRT, n = NULL)$table
SAM_YNO_results <- topTags(SAM_YNO_LRT, n = NULL)$table
SAM_MNO_results <- topTags(SAM_MNO_LRT, n = NULL)$table

YIN_MIN_results <- topTags(YIN_MIN_LRT, n = NULL)$table
YIN_YNO_results <- topTags(YIN_YNO_LRT, n = NULL)$table
YIN_MNO_results <- topTags(YIN_MNO_LRT, n = NULL)$table

MIN_YNO_results <- topTags(MIN_YNO_LRT, n = NULL)$table
MIN_MNO_results <- topTags(MIN_MNO_LRT, n = NULL)$table

YNO_MNO_results <- topTags(YNO_MNO_LRT, n = NULL)$table

#order the results by tag (transcript) ID to later build the results object
SAM_YIN_results <- SAM_YIN_results[order(row.names(SAM_YIN_results)),]
SAM_MIN_results <- SAM_MIN_results[order(row.names(SAM_MIN_results)),]
SAM_YNO_results <- SAM_YNO_results[order(row.names(SAM_YNO_results)),]
SAM_MNO_results <- SAM_MNO_results[order(row.names(SAM_MNO_results)),]

YIN_MIN_results <- YIN_MIN_results[order(row.names(YIN_MIN_results)),]
YIN_YNO_results <- YIN_YNO_results[order(row.names(YIN_YNO_results)),]
YIN_MNO_results <- YIN_MNO_results[order(row.names(YIN_MNO_results)),]

MIN_YNO_results <- MIN_YNO_results[order(row.names(MIN_YNO_results)),]
MIN_MNO_results <- MIN_MNO_results[order(row.names(MIN_MNO_results)),]

YNO_MNO_results <- YNO_MNO_results[order(row.names(YNO_MNO_results)),]
