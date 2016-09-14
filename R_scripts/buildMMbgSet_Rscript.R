
#build a new MAPMAN ontology for bamboo and create a background set for
#enrichment analysis

library(data.table)

#need to manually replace and correct some typos in the input file
#then read in the rice mappings
#I will not publishing original MM data on github, there might be legal issues.
#You have to live without that file. Just use the bg_MM_PH object.
all_MM <- read.delim("PATH",stringsAsFactors = FALSE)

#then prepare all MM identifier and a background  set for enrichment analysis
all_MM$IDENTIFIER <- toupper(all_MM$IDENTIFIER)

temp <- tstrsplit(all_MM$BINCODE,".",fixed = TRUE, fill = "undefined")
all_MM$bin1 <- as.factor(temp[[1]])
all_MM$bin2 <- as.factor(paste(all_MM$bin1,".",temp[[2]],sep = ""))
all_MM$bin3 <- as.factor(paste(all_MM$bin2,".",temp[[3]],sep = ""))
all_MM$bin4 <- as.factor(paste(all_MM$bin3,".",temp[[4]],sep = ""))
all_MM$bin5 <- as.factor(paste(all_MM$bin4,".",temp[[5]],sep = ""))
all_MM$bin6 <- as.factor(paste(all_MM$bin5,".",temp[[6]],sep = ""))
all_MM$bin7 <- as.factor(paste(all_MM$bin6,".",temp[[7]],sep = ""))
rm(temp)
temp <- tstrsplit(all_MM$NAME,".",fixed = TRUE, fill = "undefined")
all_MM$name1 <- as.factor(temp[[1]])
all_MM$name2 <- as.factor(paste(all_MM$name1,".",temp[[2]],sep = ""))
all_MM$name3 <- as.factor(paste(all_MM$name2,".",temp[[3]],sep = ""))
all_MM$name4 <- as.factor(paste(all_MM$name3,".",temp[[4]],sep = ""))
all_MM$name5 <- as.factor(paste(all_MM$name4,".",temp[[5]],sep = ""))
all_MM$name6 <- as.factor(paste(all_MM$name5,".",temp[[6]],sep = ""))
all_MM$name7 <- as.factor(paste(all_MM$name6,".",temp[[7]],sep = ""))
rm(temp)

all_MM$bin_name1 <- as.factor(paste(all_MM$bin1, all_MM$name1, sep = "_"))
all_MM$bin_name2 <- as.factor(paste(all_MM$bin2, all_MM$name2, sep = "_"))
all_MM$bin_name3 <- as.factor(paste(all_MM$bin3, all_MM$name3, sep = "_"))
all_MM$bin_name4 <- as.factor(paste(all_MM$bin4, all_MM$name4, sep = "_"))
all_MM$bin_name5 <- as.factor(paste(all_MM$bin5, all_MM$name5, sep = "_"))
all_MM$bin_name6 <- as.factor(paste(all_MM$bin6, all_MM$name6, sep = "_"))
all_MM$bin_name7 <- as.factor(paste(all_MM$bin7, all_MM$name7, sep = "_"))

#then use the all_MM set based on rice to make a bamboo bg set

#read in a lookup table OS <-> PH
PH_OS_IDS <- read.delim("./data/PH_OS_IDs.txt", stringsAsFactors = FALSE) 

#remove those Os codes not represented in MM
PH_OS_IDS <- PH_OS_IDS[PH_OS_IDS$OS_ID %in% all_MM$IDENTIFIER,]

bg_MM_PH <- data.frame(NULL)
for(rowcount in 1:nrow(PH_OS_IDS)) {
  
  bg_MM_PH <- rbind(bg_MM_PH,
                cbind("PH_ID" = PH_OS_IDS$PH_ID[rowcount],
                      all_MM[all_MM$IDENTIFIER == PH_OS_IDS$OS_ID[rowcount],]))
  
}
colnames(bg_MM_PH)[4] <- "Os_IDENTIFIER"

#save the bg set for use in the shiny app
save(bg_MM_PH, file = "./data/bg_MM_PH")

#can also be saved for use outside R
#write.table(bg_MM_PH, file = "./data/bg_MM_PH.txt",
#            sep = "\t", row.names = FALSE)

