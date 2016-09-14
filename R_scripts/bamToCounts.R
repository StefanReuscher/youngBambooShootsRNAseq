#run once
source("https://bioconductor.org/biocLite.R")
biocLite("Rsubread")
biocLite("Biostrings")

library(Rsubread)
library(Biostrings)

## Transcript quantification from bamboo ref-based transcriptome

# featureCounts() needs a dataframe with feature annotations
# build that from the CDS fasta headers,
CDS <- readDNAStringSet("./P_heterocycla_v1.0.genemodel-cds-DNA.fa")

CDSnames <- strsplit(names(CDS),"\\|")
CDSnames <- data.frame(matrix(unlist(CDSnames),nrow = length(CDS), byrow = TRUE), stringsAsFactors = FALSE)

bambooAnn <- data.frame("GeneID" = CDSnames$X1,
                        "Chr" = substr(CDSnames$X1,1,10),
                        "Start" = CDSnames$X4,
                        "End" = CDSnames$X5,
                        "Strand" = CDSnames$X7,
                        stringsAsFactors = FALSE)

# correct errors from input files manually
bambooAnn[12440, "End"] <- 385813
bambooAnn[12971, "End"] <- 564134

#sample names
tissues <- c("SAM_A","SAM_B","SAM_C","SAM_D",
             "YIN_A","YIN_B","YIN_C","YIN_D",
             "MIN_A","MIN_B","MIN_C","MIN_D",
             "YNO_A","YNO_B","YNO_C","YNO_D",
             "MNO_A","MNO_B","MNO_C","MNO_D")

# get counts, a bit not elegant, why copy every file, but works fine
for(i in tissues){
  system2("cp", args = c(
    paste("/pathToAssembly/tophat2_",i,"/accepted_hits.bam", sep = ""),
    ".")
  )
  assign(paste(i,"_counts", sep = ""),
         featureCounts(files = "accepted_hits.bam",
                       annot.ext = bambooAnn, isPairedEnd = TRUE) #quantifying transcript here
  )
  system2("rm", args = "accepted_hits.bam")
}

# make a results matrix
all_counts <- data.frame(row.names = bambooAnn$GeneID)
for (i in tissues){
  temp_counts <- get(paste(i,"_counts", sep = ""))
  all_counts <- cbind(all_counts,
        temp_counts$counts)
}
colnames(all_counts) <- tissues

#save an R objects with counts
save(all_counts, file = "./allCountsRefbased")
