####load packages####
library(shiny)
library(ggplot2)
library(reshape2)
library(grid)
require(edgeR)


shinyUI(fluidPage(
  
  includeCSS("styles.css"),
  
  titlePanel("Visualize bamboo RNAseq data from different tissues of the young shoot apex."),
  
  sidebarPanel(width = 3,
    textInput("password", label = "What is our favourite gene ?",value = ""),
    tags$i("Ask for our favourite gene at reuscher@agr.nagoya-u.ac.jp to gain access."),
    tags$hr(),
    radioButtons("typeOfAnalysis", "What kind of analyses should be performed ?",
                 choices = list("single gene-based" = "sg",
                                "DEG-based" = "deg",
                                "clustering-based" = "kg"),
                 selected = "deg"),
    
    conditionalPanel(condition = "input.typeOfAnalysis == 'sg'",
                     tags$strong("Which genes should be plotted ?"),
                     tags$div(),
                     tags$textarea(id = "genesToPlot",
                                   #placeholder = "find an some examples at the bottom",
                                   #"example:\nPhet014978.000\nPhet022707.000\nPhet024203.000\nPhet028726.000\nPhet018447.001\nPhet025628.000\nPhet035437.000\nPhet010198.000\nPhet009893.000\nPhet026624.001"),
                                   "example:\nPH01000000G0240\nPH01000000G0290\nPH01000000G0330\nPH01000000G0350"),
                                   
                     selectInput("facetswitch",
                                 label = "What should be in a single panel ?",
                                 choices = list("A single transcript" = "transcript",
                                                "One of the five analyzed tissues." = "tissue")),
                     checkboxGroupInput("whichTissue", label = "Which tissues should be plotted ?",
                                        choices = list("shoot apical meristem" = "SAM",
                                                       "young internode"       = "YIN",
                                                       "mature internode"      = "MIN",
                                                       "young node"            = "YNO",
                                                       "mature node"           = "MNO"),
                                        selected = c("SAM","YIN","MIN","YNO","MNO")),
                     numericInput("numFacetCol", label = "How many panel next to each other ?", value = "2"),
                     selectInput("plotSort",
                                 label = "How should the expression plots be sorted ?",
                                 choices = list("by input order" = "inOrder",
                                                "alphanumeric order" = "alphaOrder"))
    ),
    
    conditionalPanel(condition = "input.typeOfAnalysis == 'deg'",
                     selectInput("tissueComp",
                                 "Choose a tissue comparison:",
                                 choices =  list("SAM vs. YIN" = "SAM_YIN",
                                                 "SAM vs. MIN" = "SAM_MIN",
                                                 "SAM vs. YNO" = "SAM_YNO",
                                                 "SAM vs. MNO" = "SAM_MNO",
                                                 
                                                 "YIN vs. MIN" = "YIN_MIN",
                                                 "YIN vs. YNO" = "YIN_YNO",
                                                 "YIN vs. MNO" = "YIN_MNO",
                                                 
                                                 "MIN vs. YNO" = "MIN_YNO",
                                                 "MIN vs. MNO" = "MIN_MNO",
                                                 
                                                 "YNO vs. MNO" = "YNO_MNO")),
                     column(6,numericInput("minFC", "|log2 FC| >=", "1")),
                     column(6,numericInput("minP", "FDR_P <=", "0.05"))
    ),
    
    conditionalPanel(condition = "input.typeOfAnalysis == 'kg'",
                     tags$h3("Clustering by k-means"),
                     selectInput("k",
                                 "Number of k-means clusters (precalculated):",
                                 choices = list( "8" = "k8",
                                                "10" = "k10",
                                                "12" = "k12",
                                                "14" = "k14",
                                                "16" = "k16",
                                                "18" = "k18",
                                                "20" = "k20"),
                                 selected = "k8"),
                     em("Not sure which is best yet."),
                     tags$hr(),
                     numericInput("whichCluster", label = "Which cluster do you want to examine ?",
                                  value = "5"),
                     em("Use the cluster ID from the heatmap as input")
    ),
    
    imageOutput("bambooShootSection", inline = TRUE),
    p("SAM = shoot apical meristem"),
    p("YIN = young internode"),
    p("YNO = young node"),
    p("MIN = mature internode"),
    p("MNO = mature node")
  ),
  mainPanel(
    tabsetPanel(
      tabPanel("single gene Plot", plotOutput("expressionPlot", height = 750),
               dataTableOutput("expr_transcript_info")),
      
      tabPanel("DEG plot and table",
               column(6, plotOutput("DEG_smearplot")),
               column(6,h3(textOutput("numberofDEGs"))),
               column(12, dataTableOutput("DEG_transcript_info"))),
      
      tabPanel("K-means cluster",
               column(6, p(strong("Heatmap of all clusters"), align = "center"),
                      plotOutput("allCluster_heatmap")),
               column(6,p(strong("Detailed view of the selected cluster."), align = "center"),
                      plotOutput("selCluster_heatmap")
               ),
               hr(),
               dataTableOutput("selClust_info")
               
      ),
      
      tabPanel("MAPMAN level 1",
               conditionalPanel(condition = "input.typeOfAnalysis == 'sg'",
                                dataTableOutput("MM_lev1_sg")),
               
               conditionalPanel(condition = "input.typeOfAnalysis == 'deg'",
                                h4(textOutput("MMinfoUp1")),
                                dataTableOutput("MM_lev1_DEGup"),
                                h4(textOutput("MMinfoDown1")),
                                dataTableOutput("MM_lev1_DEGdown")),
               
               conditionalPanel(condition = "input.typeOfAnalysis == 'kg'",
                                dataTableOutput("MM_lev1_kmeans"))
               ),
      
      tabPanel("MAPMAN level 2",
               conditionalPanel(condition = "input.typeOfAnalysis == 'sg'",
                                dataTableOutput("MM_lev2_sg")),
               
               conditionalPanel(condition = "input.typeOfAnalysis == 'deg'",
                                h4(textOutput("MMinfoUp2")),
                                dataTableOutput("MM_lev2_DEGup"),
                                h4(textOutput("MMinfoDown2")),
                                dataTableOutput("MM_lev2_DEGdown")),
               
               conditionalPanel(condition = "input.typeOfAnalysis == 'kg'",
                                dataTableOutput("MM_lev2_kmeans"))
               ),
       
      tabPanel("MAPMAN level 3",
               conditionalPanel(condition = "input.typeOfAnalysis == 'sg'",
                                dataTableOutput("MM_lev3_sg")),
               
               conditionalPanel(condition = "input.typeOfAnalysis == 'deg'",
                                h4(textOutput("MMinfoUp3")),
                                dataTableOutput("MM_lev3_DEGup"),
                                h4(textOutput("MMinfoDown3")),
                                dataTableOutput("MM_lev3_DEGdown")),
               
               conditionalPanel(condition = "input.typeOfAnalysis == 'kg'",
                                dataTableOutput("MM_lev3_kmeans"))
               ),
      tabPanel("data table",
               h4("This shows RPKM values from the selected analysis"),
               dataTableOutput("expression_table"))
    )
  )
)
)
