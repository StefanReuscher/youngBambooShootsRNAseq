

####read data for expression plot####
rpkm_mean <- read.delim("./data/RPKM_mean.txt", sep = "\t", stringsAsFactors = FALSE)
no_of_genes <- nrow(rpkm_mean)
rpkm_df <- melt(rpkm_mean)
rpkm_df$tissue <- factor(c(rep("SAM", no_of_genes),
                             rep("YIN", no_of_genes),
                             rep("MIN", no_of_genes),
                             rep("YNO", no_of_genes),
                             rep("MNO", no_of_genes)),
                           levels = c("SAM","YIN","MIN","YNO","MNO"))
rpkm_df$variable <- NULL
rpkm_sd <- read.delim("data/RPKM_sd.txt", sep = "\t", stringsAsFactors = FALSE)
rpkm_sd <- melt(rpkm_sd)
rpkm_df$sd <- rpkm_sd$value
####

####read data for annotations####
annotations <- read.delim("./data/bamboo_annotations.txt", sep = "\t", stringsAsFactors = FALSE)
####

####theme for plotting####
theme_custom <- theme_grey() +
  theme(axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 14, colour = "black", face = "bold"),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.ticks =element_line(colour = "black"),
        #axis.title.x = element_text(vjust = -0.5),
        #axis.title.y = element_text(vjust = 1.5),
        panel.grid.major = element_line(colour="grey95"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5,linetype="solid"),
        strip.text = element_text(size = 14, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = "bottom"
  )
####


#####prepare MM table#####
MM_table_option <- list(
  columnDefs = list(list(targets = c(2,3,4,5,6)-1,
                         searchable = FALSE)),
  pageLength = 10)
load("./data/bg_MM_PH")
####

####read DEG related objects####
load("./data/LRT_objects")
####

####read cluster related objects####
load("./data/myRNAseqData_bamboo")
load("./data/clustID_bamboo")
####

####shiny server start####
shinyServer(function(input, output) {
#image of bamboo shoot
  output$bambooShootSection <- renderImage({
    
    list(src = "www/bambooShootSampling.png",
         contentType = "image/png",
         alt = "image file missing",
         id = "sidebarImage"
    )
  },
  deleteFile = FALSE)
  
  #expression plot
  output$expressionPlot <- renderPlot({
    
    if (input$genesToPlot == "")
      return(NULL)
    if (!is.numeric(input$numFacetCol)|input$numFacetCol == 0)
      return(NULL)
    if (input$password != "SK1"){
      return(NULL)
    }else{
      
      rpkmToPlot <- rpkm_df[rpkm_df$ID %in% unlist(strsplit(input$genesToPlot,"\n")),]
      
      rpkmToPlot <- rpkmToPlot[rpkmToPlot$tissue %in% input$whichTissue,]
      if(input$plotSort == "inOrder"){
        rpkmToPlot$ID <- factor(rpkmToPlot$ID, unlist(strsplit(input$genesToPlot,"\n")))
      }
      if (input$facetswitch == "transcript") { 
        ggplot(data=rpkmToPlot, aes(x=tissue,
                                    y=value))+
          geom_errorbar(aes(ymin = value, ymax = value + sd),
                        width = 0.2, colour = "grey25")+
          geom_bar(stat = "identity")+
          facet_wrap(~ID, scales = "free", ncol = as.integer(input$numFacetCol))+
          ylab("RPKM")+
          expand_limits(y = 0)+
          theme_custom
      } else {
        ggplot(data=rpkmToPlot, aes(x=ID,
                                    y=value))+
          geom_errorbar(aes(ymin = value, ymax = value + sd),
                        width = 0.2, colour = "black")+
          geom_bar(stat = "identity")+
          facet_wrap(~tissue, scales = "free_x", ncol = as.integer(input$numFacetCol))+
          ylab("RPKM")+
          expand_limits(y = 0)+
          theme_custom+
          theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
      }
    }
  })
  
  ####table below the single gene plot with annotations
  output$expr_transcript_info <- renderDataTable({
    if (input$password != "SK1") return(NULL)
    annotations[annotations$PH_ID %in% unlist(strsplit(input$genesToPlot,"\n")),]},
    options = list(paging = FALSE)
  )
  #####
  
  ####RPKM data for selected genes/DEGs/clusters
  output$expression_table <- renderDataTable({
    if(input$password != "SK1") return(NULL)
    
    rpkm <- read.delim("data/RPKM_mean.txt",
                       stringsAsFactors = FALSE, sep = "\t")
    rpkmToTable <- rpkm[rpkm$ID %in% genesToMM(),]
    rpkmToTable
    })
  
  #obects needs for DEG analyses
  LRTtoPlot <- reactive({get(paste(input$tissueComp,"_LRT", sep = ""))})
  
  DEGstoPlot <- reactive({rownames(LRTtoPlot())[as.logical(decideTestsDGE(
    LRTtoPlot(), p.value = input$minP, lfc = input$minFC))]})
  
  DEGstoPlotUp <- reactive({rownames(LRTtoPlot())[decideTestsDGE(
    LRTtoPlot(), p.value = input$minP, lfc = input$minFC) == 1]})
  
  DEGstoPlotDown <- reactive({rownames(LRTtoPlot())[decideTestsDGE(
    LRTtoPlot(), p.value = input$minP, lfc = input$minFC) == -1]})
  
  #smear plot
  output$DEG_smearplot <- renderPlot({
    if (input$password != "SK1") return(NULL)
    
    plotSmear(LRTtoPlot(), de.tags = DEGstoPlot(), cex = 1, col=rgb(0.5, 0.5, 0.5, 0.3))
    abline(h = c(input$minFC,-input$minFC), col = "blue", lty = 3)
  })
  
  #DEGs numbers
  output$numberofDEGs <- renderText({
    
    numberUp <- table(decideTestsDGE(LRTtoPlot(),
                                     p.value = input$minP,
                                     lfc = input$minFC))["1"]
    numberDown <- table(decideTestsDGE(LRTtoPlot(),
                                       p.value = input$minP,
                                       lfc = input$minFC))["-1"]
    
    paste("There are ",numberUp," transcripts with significantly higher levels and ", 
          numberDown," transcript with significantly lower levels in the ",
          substr(input$tissueComp,1,3)," compared to the ", substr(input$tissueComp,5,7),
          " (log2-fold change >= ",input$minFC,
          " and corrected P-value <= ",input$minP,").",
          sep = "")
  })
  
  ####table DEG plot####
  output$DEG_transcript_info <- renderDataTable({
    if (input$password != "SK1") return(NULL)
    
    DEG_annotations <- annotations[annotations$PH_ID %in% DEGstoPlot(),]
    DEG_annotations$log2FC <- LRTtoPlot()$table[DEGstoPlot(),"logFC"]
    DEG_annotations$FDR_P <- topTags(LRTtoPlot(), n = NULL)$table[DEGstoPlot(),"FDR"]
    DEG_annotations[order(DEG_annotations$FDR_P, decreasing = FALSE),]},
    options = list(paging = TRUE))
  #####
  
  ####plot all clusters
  output$allCluster_heatmap <- renderPlot({
    if (input$password != "SK1"){
      return(NULL)
    }else{
      FC_heatmap_DF <- as.data.frame(myRNASeqData$logFC)
      colnames(FC_heatmap_DF) <- c("SAM","YIN","MIN","YNO","MNO")
      FC_heatmap_DF$clustID <- clustID[,input$k] #select k here
      FC_heatmap_DF <- aggregate(FC_heatmap_DF,list(FC_heatmap_DF$clustID),mean)
      FC_heatmap_DF <- melt(FC_heatmap_DF[,-1],id.vars = "clustID")
      
      ggplot(FC_heatmap_DF, aes(x = variable, y = clustID))+
        geom_tile(aes(fill = value), color = "black")+
        xlab("tissue")+ 
        ylab("cluster ID")+
        #scale_x_discrete(labels = c("0 h","1 h","3 h","6 h","12h","24 h"))+
        #scale_y_discrete(breaks = seq(from = 1, to = max(clustID)))+
        scale_fill_gradient2(low = "blue", high = "red",
                             name = expression(atop(average~log["2"]~FC,
                                                    relative~to~row~mean)))+
        theme(axis.text = element_text(size = 14, color = "black"),
              axis.title = element_text(size = 14, face = "bold"),
              axis.ticks = element_blank(),
              panel.background = element_rect(fill = NA),
              legend.title = element_text(size = 12),
              legend.text = element_text(size = 12))
      
    }
  })
  #####
  
  clustID_genesToPlot <- reactive({
    myRNASeqData$GeneID[clustID[,input$k] == input$whichCluster]
  })
  
  #####plot selected cluster#####
  output$selCluster_heatmap <- renderPlot({
    if (input$password != "SK1") return(NULL)
    
    selCluster_DF <- rpkm_mean[rpkm_mean$ID %in% clustID_genesToPlot(),2:6]
    selCluster_DF <- as.data.frame(t(selCluster_DF))
    selCluster_DF <- cbind("tissue" = c("SAM","YIN","MIN","YNO","MNO"),
                            selCluster_DF)
    
    selCluster_DF <- melt(selCluster_DF)
    selCluster_DF$tissue <- factor(selCluster_DF$tissue,
                                    levels = c("SAM","YIN","MIN","YNO","MNO"))
    
    ggplot(data = selCluster_DF, aes(x = tissue, y = log10(value+0.01), group = tissue))+
      #geom_line(alpha = 0.5)+
      geom_boxplot()+
      #scale_x_discrete(labels = c("SAM","YIN","MIN","YNO","MNO"))+
      xlab("tissue")+
      ylab(expression(log["10"]~RPKM))+
      #ylab(expression(average~log["2"]~FC~relative~to~row~mean))+
      ggtitle(paste("Cluster ID: ",input$whichCluster,
                    " Number of genes: ",length(levels(selCluster_DF$variable))))+
      theme(plot.title = element_text(size = 14, face = "bold"),
            axis.text = element_text(size = 14, color = "black"),
            axis.title = element_text(size = 14, face = "bold"),
            axis.line = element_line(colour = "black", size = 0.5),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(colour = "black"),
            panel.border = element_rect(fill= NA, colour="grey60"),
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            strip.text = element_text(size = 14))
  })
  
  #####table for selected cluster
  output$selClust_info <- renderDataTable({
    if (input$password != "SK1") return(NULL)
    
    selClust_annotations <- annotations[annotations$PH_ID %in% clustID_genesToPlot(),]
    selClust_annotations <- cbind(selClust_annotations,
                                  round(rpkm_mean[rpkm_mean$ID %in% clustID_genesToPlot(),2:6],
                                        2)
                                  )
    },
    options = list(paging = TRUE)
  )
  
  
  #function to select genes in the MAPMAN analysis
  genesToMM <- reactive({
    switch(input$typeOfAnalysis,
           "sg" = unlist(strsplit(input$genesToPlot,"\n")),
           "deg" = DEGstoPlot(), #only used for expression table
           "kg" = clustID_genesToPlot()) 
  })
  
####MM table headers####
  output$MMinfoUp1 <- renderText({
    paste("Overrepresented MM bins among genes significantly HIGHER expressed in ",
          substr(input$tissueComp,1,3)," compared to ",substr(input$tissueComp,5,7),
          ".", sep = "")
  }) #cant reuse same output
  output$MMinfoDown1 <- renderText({
    paste("Overrepresented MM bins among genes significantly LOWER expressed in ",
          substr(input$tissueComp,1,3)," compared to ",substr(input$tissueComp,5,7),
          ".", sep = "")
  })
  output$MMinfoUp2 <- renderText({
    paste("Overrepresented MM bins among genes significantly HIGHER expressed in ",
          substr(input$tissueComp,1,3)," compared to ",substr(input$tissueComp,5,7),
          ".", sep = "")
  })
  output$MMinfoDown2 <- renderText({
    paste("Overrepresented MM bins among genes significantly LOWER expressed in ",
          substr(input$tissueComp,1,3)," compared to ",substr(input$tissueComp,5,7),
          ".", sep = "")
  })
  output$MMinfoUp3 <- renderText({
    paste("Overrepresented MM bins among genes significantly HIGHER expressed in ",
          substr(input$tissueComp,1,3)," compared to ",substr(input$tissueComp,5,7),
          ".", sep = "")
  })
  output$MMinfoDown3 <- renderText({
    paste("Overrepresented MM bins among genes significantly LOWER expressed in ",
          substr(input$tissueComp,1,3)," compared to ",substr(input$tissueComp,5,7),
          ".", sep = "")
  })
#####
  
####MM tables####
  MMfun_lev1 <- reactive({
    if (input$password != "SK1") return(NULL)
    #if (genesToMM() == "") return(NULL)
    
    genesToMM <- genesToMM()
    set_MM <- bg_MM_PH[bg_MM_PH$PH_ID %in% genesToMM,]
    testDF <- data.frame( "BinSet" = table(set_MM$bin_name1),
                          "BinBg" = as.vector(table(bg_MM_PH$bin_name1)),
                          "NotBinSet" =  nrow(set_MM) - as.vector(table(set_MM$bin_name1)),
                          "NotBinBg" = nrow(bg_MM_PH) - as.vector(table(bg_MM_PH$bin_name1)))
    testDF$OR_FC <- round((testDF[,2] / nrow(set_MM)) / (testDF[,3] / nrow(bg_MM_PH)), 2)
    testDF <- testDF[testDF$BinSet.Freq > 0,] #remove bins with 0 genes
    testDF <- testDF[!grepl(pattern = "undefined", testDF$BinSet.Var1),] #remove undefined sub bins
    
    FET_P <- NULL
    for(i in 1:nrow(testDF)) {
      
      BinSet <- testDF[i,2]
      BinBg <- testDF[i,3]
      NotBinSet <- testDF[i,4]
      NotBinBg <- testDF[i,5]
      
      FET_P[i] <- fisher.test(matrix(c(BinSet,BinBg,NotBinSet,NotBinBg),nrow = 2),
                              alternative = "greater")$p.value
    }
    testDF$FET_P <- round(FET_P, 3)
    testDF$FDR_P <- round(p.adjust(FET_P, method = "BH"), 3)
    colnames(testDF) <- c("MAPMAN level 1 bin",
                          "Genes in set", "Genes in total",
                          "Other genes in set", "Other genes in total",
                          "Fold enrichement", "P-value (FET)", "corrected P-value (BH)")
    return(testDF[order(testDF[,8]),])
  }) #not really shorter than single functions
  MMfun_lev2 <- reactive({
    if (input$password != "SK1") return(NULL)
    #if (genesToMM() == "") return(NULL)
    
    genesToMM <- genesToMM()
    set_MM <- bg_MM_PH[bg_MM_PH$PH_ID %in% genesToMM,]
    testDF <- data.frame( "BinSet" = table(set_MM$bin_name2),
                          "BinBg" = as.vector(table(bg_MM_PH$bin_name2)),
                          "NotBinSet" =  nrow(set_MM) - as.vector(table(set_MM$bin_name2)),
                          "NotBinBg" = nrow(bg_MM_PH) - as.vector(table(bg_MM_PH$bin_name2)))
    testDF$OR_FC <- round((testDF[,2] / nrow(set_MM)) / (testDF[,3] / nrow(bg_MM_PH)), 2)
    testDF <- testDF[testDF$BinSet.Freq > 0,] #remove bins with 0 genes
    testDF <- testDF[!grepl(pattern = "undefined", testDF$BinSet.Var1),] #remove undefined sub bins
    
    FET_P <- NULL
    for(i in 1:nrow(testDF)) {
      
      BinSet <- testDF[i,2]
      BinBg <- testDF[i,3]
      NotBinSet <- testDF[i,4]
      NotBinBg <- testDF[i,5]
      
      FET_P[i] <- fisher.test(matrix(c(BinSet,BinBg,NotBinSet,NotBinBg),nrow = 2),
                              alternative = "greater")$p.value
    }
    testDF$FET_P <- round(FET_P, 3)
    testDF$FDR_P <- round(p.adjust(FET_P, method = "BH"), 3)
    colnames(testDF) <- c("MAPMAN level 2 bin",
                          "Genes in set", "Genes in total",
                          "Other genes in set", "Other genes in total",
                          "Fold enrichement", "P-value (FET)", "corrected P-value (BH)")
    return(testDF[order(testDF[,8]),])
  })
  MMfun_lev3 <- reactive({
    if (input$password != "SK1") return(NULL)
    #if (genesToMM() == "") return(NULL)
    
    genesToMM <- genesToMM()
    set_MM <- bg_MM_PH[bg_MM_PH$PH_ID %in% genesToMM,]
    testDF <- data.frame( "BinSet" = table(set_MM$bin_name3),
                          "BinBg" = as.vector(table(bg_MM_PH$bin_name3)),
                          "NotBinSet" =  nrow(set_MM) - as.vector(table(set_MM$bin_name3)),
                          "NotBinBg" = nrow(bg_MM_PH) - as.vector(table(bg_MM_PH$bin_name3)))
    testDF$OR_FC <- round((testDF[,2] / nrow(set_MM)) / (testDF[,3] / nrow(bg_MM_PH)), 2)
    testDF <- testDF[testDF$BinSet.Freq > 0,] #remove bins with 0 genes
    testDF <- testDF[!grepl(pattern = "undefined", testDF$BinSet.Var1),] #remove undefined sub bins
    
    FET_P <- NULL
    for(i in 1:nrow(testDF)) {
      
      BinSet <- testDF[i,2]
      BinBg <- testDF[i,3]
      NotBinSet <- testDF[i,4]
      NotBinBg <- testDF[i,5]
      
      FET_P[i] <- fisher.test(matrix(c(BinSet,BinBg,NotBinSet,NotBinBg),nrow = 2),
                              alternative = "greater")$p.value
    }
    testDF$FET_P <- round(FET_P, 3)
    testDF$FDR_P <- round(p.adjust(FET_P, method = "BH"), 3)
    colnames(testDF) <- c("MAPMAN level 3 bin",
                          "Genes in set", "Genes in total",
                          "Other genes in set", "Other genes in total",
                          "Fold enrichement", "P-value (FET)", "corrected P-value (BH)")
    return(testDF[order(testDF[,8]),])
  })
  
  output$MM_lev1_sg <-renderDataTable({MMfun_lev1()}, options = MM_table_option)
  output$MM_lev1_DEGup <- renderDataTable({
    if (input$password != "SK1") return(NULL)
    
    genesToMM <- DEGstoPlotUp()
    set_MM <- bg_MM_PH[bg_MM_PH$PH_ID %in% genesToMM,]
    testDF <- data.frame( "BinSet" = table(set_MM$bin_name1),
                          "BinBg" = as.vector(table(bg_MM_PH$bin_name1)),
                          "NotBinSet" =  nrow(set_MM) - as.vector(table(set_MM$bin_name1)),
                          "NotBinBg" = nrow(bg_MM_PH) - as.vector(table(bg_MM_PH$bin_name1)))
    testDF$OR_FC <- round((testDF[,2] / nrow(set_MM)) / (testDF[,3] / nrow(bg_MM_PH)), 2)
    testDF <- testDF[testDF$BinSet.Freq > 0,] #remove bins with 0 genes
    testDF <- testDF[!grepl(pattern = "undefined", testDF$BinSet.Var1),] #remove undefined sub bins
    
    FET_P <- NULL
    for(i in 1:nrow(testDF)) {
      
      BinSet <- testDF[i,2]
      BinBg <- testDF[i,3]
      NotBinSet <- testDF[i,4]
      NotBinBg <- testDF[i,5]
      
      FET_P[i] <- fisher.test(matrix(c(BinSet,BinBg,NotBinSet,NotBinBg),nrow = 2),
                              alternative = "greater")$p.value
    }
    testDF$FET_P <- round(FET_P, 3)
    testDF$FDR_P <- round(p.adjust(FET_P, method = "BH"), 3)
    colnames(testDF) <- c("MAPMAN level 1 bin",
                          "Genes in set", "Genes in total",
                          "Other genes in set", "Other genes in total",
                          "Fold enrichement", "P-value (FET)", "corrected P-value (BH)")
    testDF[order(testDF[,8]),]
  }, options = MM_table_option)
  output$MM_lev1_DEGdown <- renderDataTable({
    if (input$password != "SK1") return(NULL)
    
    genesToMM <- DEGstoPlotDown()
    
    set_MM <- bg_MM_PH[bg_MM_PH$PH_ID %in% genesToMM,]
    testDF <- data.frame( "BinSet" = table(set_MM$bin_name1),
                          "BinBg" = as.vector(table(bg_MM_PH$bin_name1)),
                          "NotBinSet" =  nrow(set_MM) - as.vector(table(set_MM$bin_name1)),
                          "NotBinBg" = nrow(bg_MM_PH) - as.vector(table(bg_MM_PH$bin_name1)))
    testDF$OR_FC <- round((testDF[,2] / nrow(set_MM)) / (testDF[,3] / nrow(bg_MM_PH)), 2)
    testDF <- testDF[testDF$BinSet.Freq > 0,] #remove bins with 0 genes
    testDF <- testDF[!grepl(pattern = "undefined", testDF$BinSet.Var1),] #remove undefined sub bins
    
    FET_P <- NULL
    for(i in 1:nrow(testDF)) {
      
      BinSet <- testDF[i,2]
      BinBg <- testDF[i,3]
      NotBinSet <- testDF[i,4]
      NotBinBg <- testDF[i,5]
      
      FET_P[i] <- fisher.test(matrix(c(BinSet,BinBg,NotBinSet,NotBinBg),nrow = 2),
                              alternative = "greater")$p.value
    }
    testDF$FET_P <- round(FET_P, 3)
    testDF$FDR_P <- round(p.adjust(FET_P, method = "BH"), 3)
    colnames(testDF) <- c("MAPMAN level 1 bin",
                          "Genes in set", "Genes in total",
                          "Other genes in set", "Other genes in total",
                          "Fold enrichement", "P-value (FET)", "corrected P-value (BH)")
    testDF[order(testDF[,8]),]
    
  }, options = MM_table_option)
  output$MM_lev1_kmeans <-renderDataTable({MMfun_lev1()}, options = MM_table_option)
  
  output$MM_lev2_sg <-renderDataTable({MMfun_lev2()}, options = MM_table_option)
  output$MM_lev2_DEGup <- renderDataTable({
    if (input$password != "SK1") return(NULL)
    
    genesToMM <- DEGstoPlotUp()
    set_MM <- bg_MM_PH[bg_MM_PH$PH_ID %in% genesToMM,]
    testDF <- data.frame( "BinSet" = table(set_MM$bin_name2),
                          "BinBg" = as.vector(table(bg_MM_PH$bin_name2)),
                          "NotBinSet" =  nrow(set_MM) - as.vector(table(set_MM$bin_name2)),
                          "NotBinBg" = nrow(bg_MM_PH) - as.vector(table(bg_MM_PH$bin_name2)))
    testDF$OR_FC <- round((testDF[,2] / nrow(set_MM)) / (testDF[,3] / nrow(bg_MM_PH)), 2)
    testDF <- testDF[testDF$BinSet.Freq > 0,] #remove bins with 0 genes
    testDF <- testDF[!grepl(pattern = "undefined", testDF$BinSet.Var1),] #remove undefined sub bins
    
    FET_P <- NULL
    for(i in 1:nrow(testDF)) {
      
      BinSet <- testDF[i,2]
      BinBg <- testDF[i,3]
      NotBinSet <- testDF[i,4]
      NotBinBg <- testDF[i,5]
      
      FET_P[i] <- fisher.test(matrix(c(BinSet,BinBg,NotBinSet,NotBinBg),nrow = 2),
                              alternative = "greater")$p.value
    }
    testDF$FET_P <- round(FET_P, 3)
    testDF$FDR_P <- round(p.adjust(FET_P, method = "BH"), 3)
    colnames(testDF) <- c("MAPMAN level 2 bin",
                          "Genes in set", "Genes in total",
                          "Other genes in set", "Other genes in total",
                          "Fold enrichement", "P-value (FET)", "corrected P-value (BH)")
    testDF[order(testDF[,8]),]
  }, options = MM_table_option)
  output$MM_lev2_DEGdown <- renderDataTable({
    if (input$password != "SK1") return(NULL)
    
    genesToMM <- DEGstoPlotDown()
    
    set_MM <- bg_MM_PH[bg_MM_PH$PH_ID %in% genesToMM,]
    testDF <- data.frame( "BinSet" = table(set_MM$bin_name2),
                          "BinBg" = as.vector(table(bg_MM_PH$bin_name2)),
                          "NotBinSet" =  nrow(set_MM) - as.vector(table(set_MM$bin_name2)),
                          "NotBinBg" = nrow(bg_MM_PH) - as.vector(table(bg_MM_PH$bin_name2)))
    testDF$OR_FC <- round((testDF[,2] / nrow(set_MM)) / (testDF[,3] / nrow(bg_MM_PH)), 2)
    testDF <- testDF[testDF$BinSet.Freq > 0,] #remove bins with 0 genes
    testDF <- testDF[!grepl(pattern = "undefined", testDF$BinSet.Var1),] #remove undefined sub bins
    
    FET_P <- NULL
    for(i in 1:nrow(testDF)) {
      
      BinSet <- testDF[i,2]
      BinBg <- testDF[i,3]
      NotBinSet <- testDF[i,4]
      NotBinBg <- testDF[i,5]
      
      FET_P[i] <- fisher.test(matrix(c(BinSet,BinBg,NotBinSet,NotBinBg),nrow = 2),
                              alternative = "greater")$p.value
    }
    testDF$FET_P <- round(FET_P, 3)
    testDF$FDR_P <- round(p.adjust(FET_P, method = "BH"), 3)
    colnames(testDF) <- c("MAPMAN level 2 bin",
                          "Genes in set", "Genes in total",
                          "Other genes in set", "Other genes in total",
                          "Fold enrichement", "P-value (FET)", "corrected P-value (BH)")
    testDF[order(testDF[,8]),]
    
  }, options = MM_table_option)
  output$MM_lev2_kmeans <-renderDataTable({MMfun_lev2()}, options = MM_table_option)

  output$MM_lev3_sg <-renderDataTable({MMfun_lev3()}, options = MM_table_option)
  output$MM_lev3_DEGup <- renderDataTable({
    if (input$password != "SK1") return(NULL)
    
    genesToMM <- DEGstoPlotUp()
    set_MM <- bg_MM_PH[bg_MM_PH$PH_ID %in% genesToMM,]
    testDF <- data.frame( "BinSet" = table(set_MM$bin_name3),
                          "BinBg" = as.vector(table(bg_MM_PH$bin_name3)),
                          "NotBinSet" =  nrow(set_MM) - as.vector(table(set_MM$bin_name3)),
                          "NotBinBg" = nrow(bg_MM_PH) - as.vector(table(bg_MM_PH$bin_name3)))
    testDF$OR_FC <- round((testDF[,2] / nrow(set_MM)) / (testDF[,3] / nrow(bg_MM_PH)), 2)
    testDF <- testDF[testDF$BinSet.Freq > 0,] #remove bins with 0 genes
    testDF <- testDF[!grepl(pattern = "undefined", testDF$BinSet.Var1),] #remove undefined sub bins
    
    FET_P <- NULL
    for(i in 1:nrow(testDF)) {
      
      BinSet <- testDF[i,2]
      BinBg <- testDF[i,3]
      NotBinSet <- testDF[i,4]
      NotBinBg <- testDF[i,5]
      
      FET_P[i] <- fisher.test(matrix(c(BinSet,BinBg,NotBinSet,NotBinBg),nrow = 2),
                              alternative = "greater")$p.value
    }
    testDF$FET_P <- round(FET_P, 3)
    testDF$FDR_P <- round(p.adjust(FET_P, method = "BH"), 3)
    colnames(testDF) <- c("MAPMAN level 3 bin",
                          "Genes in set", "Genes in total",
                          "Other genes in set", "Other genes in total",
                          "Fold enrichement", "P-value (FET)", "corrected P-value (BH)")
    testDF[order(testDF[,8]),]
  }, options = MM_table_option)
  output$MM_lev3_DEGdown <- renderDataTable({
    if (input$password != "SK1") return(NULL)
    
    genesToMM <- DEGstoPlotDown()
    
    set_MM <- bg_MM_PH[bg_MM_PH$PH_ID %in% genesToMM,]
    testDF <- data.frame( "BinSet" = table(set_MM$bin_name3),
                          "BinBg" = as.vector(table(bg_MM_PH$bin_name3)),
                          "NotBinSet" =  nrow(set_MM) - as.vector(table(set_MM$bin_name3)),
                          "NotBinBg" = nrow(bg_MM_PH) - as.vector(table(bg_MM_PH$bin_name3)))
    testDF$OR_FC <- round((testDF[,2] / nrow(set_MM)) / (testDF[,3] / nrow(bg_MM_PH)), 2)
    testDF <- testDF[testDF$BinSet.Freq > 0,] #remove bins with 0 genes
    testDF <- testDF[!grepl(pattern = "undefined", testDF$BinSet.Var1),] #remove undefined sub bins
    
    FET_P <- NULL
    for(i in 1:nrow(testDF)) {
      
      BinSet <- testDF[i,2]
      BinBg <- testDF[i,3]
      NotBinSet <- testDF[i,4]
      NotBinBg <- testDF[i,5]
      
      FET_P[i] <- fisher.test(matrix(c(BinSet,BinBg,NotBinSet,NotBinBg),nrow = 2),
                              alternative = "greater")$p.value
    }
    testDF$FET_P <- round(FET_P, 3)
    testDF$FDR_P <- round(p.adjust(FET_P, method = "BH"), 3)
    colnames(testDF) <- c("MAPMAN level 3 bin",
                          "Genes in set", "Genes in total",
                          "Other genes in set", "Other genes in total",
                          "Fold enrichement", "P-value (FET)", "corrected P-value (BH)")
    testDF[order(testDF[,8]),]
    
  }, options = MM_table_option)
  output$MM_lev3_kmeans <-renderDataTable({MMfun_lev3()}, options = MM_table_option)
#####  

})

