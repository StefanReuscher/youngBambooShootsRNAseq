# Transcriptome analysis of young bamboo shoots
<p align="center">
<img align="center" src="https://github.com/StefanReuscher/youngBambooShootsRNAseq/blob/master/www/take_kamon_alu.png" width="500">
</p>

This repository contains the source code and most data to reproduce the analyses presented in:  
*Gamuyao et al. 2016, Running title: Bamboo hormones and transcriptomics.*

It also contains the source code for a R/shiny application that allows to do your own analyses.

### Quick start
To just start the data mining application you could head over to [shinyapps.io](https://reuscher.shinyapps.io/Bamboo_RNAseq_datamining/) and start. To start the app in this repo run from within R:  
```R
library(shiny)
shiny::runGitHub("StefanReuscher/youngBambooShootsRNAseq")
```

It might take a while to download all the data, so I suggest to clone the repo first and then run the shiny app locally.


### Further information
If you are interested in how we analyzed the data from start to finish please have a look at the wiki. It also contains more detailed information on installing, running and understanding the app.
