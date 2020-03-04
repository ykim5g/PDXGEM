library(shiny)
library(DT)
library( pheatmap)
library( randomForest)
library( car )
library(genefilter)
library(siggenes)
library(multtest)

require(dplyr)

homedir <- "./"
setwd(homedir)
#source("./library/library-2016-01-25_new.R")
#source("./library/library_appendix_2017-08-10.R")

drug.name <- "paclitaxel"

### Load Patient Data Sets
load("./Combo Chemo Data/breast/hess133.total.RData")

probe.u133a <- rownames(hess133)
  hess133.response
hess133.clinic <- data.frame(response=hess133.response, row.names=colnames(hess133))
hess133.clinic$is.valid=1

load("./Combo Chemo Data/breast/hess100.total.RData")
head( cbind( rownames(hess100.clinic), colnames(hess100), names(hess100.response)) )
hess100.clinic$is.valid <- 1
rownames(hess100.clinic) <- colnames(hess100)
names(hess100.response) <- colnames(hess100)



## array data
load(file="./RData/GSE78806/PDX.array.total.RData" )
  dim(PDX.array) ## 54,675 x 661
  head( PDX.clinic )
  table( PDX.clinic$passage )
  cbind( colnames(PDX.array), PDX.clinic$geo_accession )

### Making Mouse model number and IDs
xx <- unlist( lapply( strsplit(PDX.clinic$title, "-"), function(x)x[1]) )
id.short <- which(nchar(xx)<4)
xx[id.short]

xx[id.short] <- paste("0", xx[id.short], sep="")
xx <- paste("X-", xx, sep="")
#table( table(PDX.clinic$Model))
PDX.clinic$Model <- xx

id.set <- vector()
unique.model <- unique(PDX.clinic$Model)
for( s in unique.model )
{
  #s <- unique.model[1]
  x <- subset(PDX.clinic, Model==s )
  if( max(x$passage)>=4 ) id.set <- c(id.set, rownames(x[which.max(x$passage), ]) )
}

id.set
length(id.set) ## 117
PDX.clinic[ id.set, ]$passage
PDX.array <- PDX.array[, id.set]
PDX.clinic <- PDX.clinic[id.set, ]

##### clinic data
PDX.curve <- read.table(file="./RData/GSE78806/PCT.curve.txt", sep="\t", header=T, as.is=T )
head( PDX.curve)
  dim(PDX.curve)
  length( table(PDX.curve$Treatment))
with(subset(PDX.curve, Treatment.type=="combo"), length(table(Treatment)))
with(subset(PDX.curve, Treatment.type=="single"), length(table(Treatment)))

## No duplication
table(PDX.clinic$Model) ## X-2353, X-2967
table(table(PDX.clinic$Model)) ## 117
subset( PDX.clinic, Model %in% c("X-2353", "X-2967") )




###### a function to generate EuroPDX dataset for a given drug
###### Generate PDX Discovery data set

load("./RData/rma.RData")



########### Novartis Data  ###############
get.PDX.data <- function( compound.name ){
  
  #compound.name <- "paclitaxel"
  a <- PDX.curve
  #a$Treatment[ grep(compound.name, a$Treatment) ]
  #names( table( a$Treatment) )
  
  ## subset of PDX.curve(drug response) on a given compound name
  PDX.paclitaxel <- subset( a, Treatment==compound.name)
  #dim( PDX.paclitaxel) ## 67 x 11
  #head( PDX.paclitaxel)
  
  ## Match mouse models of clinic data to expression data   
  id.match <- match(PDX.clinic$Model, PDX.paclitaxel$Model)
  #id.match
  #table( PDX.curve$Model[ is.na(match(PDX.curve$Model, PDX.clinic$Model))] )
  PDX.clinic <- cbind( PDX.clinic, PDX.paclitaxel[id.match, ] )
  
  ## Get Best response
  PDX.clinic$response <- PDX.clinic$BestResponse
  
  ## valid data
  id.valid <- which( !is.na(PDX.clinic$response) )
  
  probe.init <- rownames(rma) ; #rm( rma)
  #probe.init <- 1:nrow(PDX.array)
  
  PDX.expr.drug <- PDX.array[probe.init, id.valid]
  PDX.clinic.drug <- PDX.clinic[id.valid, ]
  
  ### % change in tumor volume (Best Average Response)
  PDX.drug.response <- PDX.clinic.drug$BestResponse
  names(PDX.drug.response) <- colnames(PDX.expr.drug)
  dim(PDX.expr.drug) # 22215 x 21 (54675 x 21)
  dim(PDX.clinic) # 117 x 71
  
  return( list( expr=PDX.expr.drug, drug.response=PDX.drug.response, clinic=PDX.clinic.drug ) )
}




get.CCCset <- function(origin="Breast"){
  filename=paste("./RData/CCCset/", origin, "-CCC.rds", sep="")
  cat(filename, "\n")
  x <- readRDS( filename )
  x <- list(expr=as.matrix(x))
  return (x)
}




source("./library/get.PDX.data.R")
source("./library/scaleExpressionData.v02.R")
source("./library/scaleExpressionData.R")
source("./library/get.dataSet.R")


###### Gene annotation of U133p2 to add gene symbol to individual probeset
###saveRDS(u133p2.gene.symbol,  file="./RData/u133p2.gene.symbol.name.rds")
##u133p2.gene.symbol <- readRDS(file="./RData/u133p2.gene.symbol.rds")
u133p2.gene.symbol <- readRDS(file="./RData/u133p2.gene.symbol.name.rds")



# Define server logic for random distribution application
shinyServer(function(input, output, session) {
  
  # Reactive expression to generate the requested distribution. This is 
  # called whenever the inputs change. The output renderers defined 
  # below then all used the value computed from this expression
  
 
  
  get.testData <- reactive({
    #get.dataSet <- function(data.name="Hess-133", compound.name="paclitaxel"){
    get.dataSet(data.name=input$testData, compound.name=input$drug ) 
  })
  
  
  
  get.validData <- reactive({
    #get.dataSet <- function(data.name="Hess-133", compound.name="paclitaxel"){
    get.dataSet(data.name=input$validData, compound.name=input$drug ) 
  })
  
 
  
  
  output$drug.name <- renderPrint({
    input$drug
  })
  
  
  
  # Generate Makrer Set Menu
  output$uiDrug <- renderUI({
    #x <- getTestData(input$MarkerSet)
    x <- input$MarkerSet
    drug.available <- c("paclitaxel") 
    
    if(x=="Novartis"){
      drug.euroPDX <- names( table( PDX.curve$Treatment) )
      drug.available <- c("paclitaxel", "5FU", "cetuximab", "gemcitabine-50mpk", "erlotinib", "trastuzumab",
                          "LDK378",
                          sort(drug.euroPDX))
    }else if(x=="EuroPDX-Breast"){
      x <- readRDS(file="./RData/EuroPDX_Breast/EuroPDX_Breast_drugSelection.rds")
      drug.available <- c("Paclitaxel", "5-Fluorouracil", "Cyclophosphamide", x$Drug)
    }  
    
    ##drug.available <- names( table( a$Treatment) )
    ### Send this input UI to a client
    selectInput(inputId="drug", label="Anti-Cancer Agent", 
                choices=c("else", drug.available ), selected="else") #selected="else")
  })
  
  
  
  ### getMarkerSet
  get.MarkerSet <- reactive({
    #print( input$BiomarkerDiscoveryMethod)
    #print( input$MarkerSet )
    #print( input$drug )    
    
    compound.name <- input$drug
    
    ## allow to update a list of drugs when markerset is 'CCLE'.     
    if(input$MarkerSet=="Novartis"){
      cat("replace data set with ", input$MarkerSet, "\n")
      x <- get.PDX.data( compound.name=compound.name)
      ##list(expr=x$expr, drug.response=x$drug.response, clinic=x$clinic)
      #rm(x); gc()
    }else if(input$MarkerSet=="EuroPDX-Breast"){
      x <- get.PDX.data.Bruna( compound.name=compound.name)
    }
    
    filename<- paste("./RData/rds/MarkerSet", input$MarkerSet, input$drug, ".rds", sep="_")
    saveRDS(x, file=filename )
    cat("filename of markerset is ", filename, "\n")
    x
  })
  
  
  
  # Generate Makrer Set Menu
  output$uiTissue <- renderUI({
    #x <- getTestData(input$MarkerSet)
    x <- get.MarkerSet()
    tissueSite <- unique(x$clinic$tissue)
    ### Send this input UI to a client
    selectInput(inputId="tissueSite", label="Select PDX tumor sites being used for Biomarker Discovery", 
                choices=c("all", tissueSite, "else" ), selected="all")
  })
  
  
  
  ### Draw Drug Response Profiles..
  output$step1 <- renderPlot({
    if( input$drug!='else'){
      ## update markerset with the reactivated one (get.MarkerSet)
      x <- get.MarkerSet()
      # x <- get.PDX.data("paclitaxel")
      
      if( !(input$tissueSite %in% c("all", "else")) ){
        id <- which( x$clinic$tissue==input$tissueSite)
      }else if( input$tissueSite %in% c("all") ) {
        id <- 1:nrow(x$clinic)
        #id <- which( x$clinic$tissue==input$tissueSite)
      }
      
      markerset <- x$expr[,id]
      markerset.response <- x$drug.response[id]
      markerset.clinic  <- x$clinic[id, ]
      #print( markerset.response)
      
      ### Draw drug response of the updated markerset
      par( mfrow=c(1,2))
      # markerset.clinic$tissue
      tissuetype <- unique(markerset.clinic$tissue)
      col.set <- c("red", "orange", "yellow", "green", "blue", "violet", "darkgray")[1:length(tissuetype)]
      names(col.set) <- tissuetype
      id.order <- order(markerset.response)
      
      plot(markerset.response[id.order], pch=16, col=col.set[markerset.clinic$tissue[id.order]], 
           cex=1.5, type="h") # sort(markerset.response))
      points(markerset.response[id.order], pch=16, col=col.set[markerset.clinic$tissue[id.order]], 
           cex=1.5) # sort(markerset.response))
      abline(h=0, lty=2)
      legend("bottomright", legend=tissuetype, col=col.set, pch=16, cex=1)
      
      ### Boxplot of drug sensitivity by cancer types...
      bb <- boxplot( markerset.response ~ markerset.clinic$tissue, las=2)
      abline(h=0, lty=2)
      mtext( text=paste("n=", bb$n, sep=""), side=3, line=0.5, at=1:length(bb$n) )
      stripchart( markerset.response ~ markerset.clinic$tissue, vertical=T, method="jitter", add=T, pch=16)
      title(main=input$drug, sub=input$MarkerSet)
    }#else{
      #plot(0, 0, type="n")
      #text("topleft", labels = "Select an anticancer drug and PDX tumor sites")
    #}
  })
  
  
  ## correlaiton biomarker
  source("./library/biomarker.discovery.correlation.R")
  ## t-test biomarker
  source("./library/biomarker.discovery.ttest.v02.R")

  
  ### Idea 1: Why don't we compare drug response profiles between different types of cancers?
  ### If no difference, we just use all rather than using one cancer type.  
  ### Biomarker Discovery 
  get.Biomarker <- reactive({
    
    a <- paste(input$drug, input$MarkerSet,  input$tissueSite, input$BiomarkerDiscoveryMethod, sep="_")
    a <- paste(a, ".rds", sep="")
    #a <- "paclitaxel_EuroPDX__correlation_pearson.rds"
    filename=paste("./RData/rds/", a, sep="")
    cat("filename is ", filename, "\n")
    
    
    #### If a biomarker was developed previously for this drug, use the saved outpu to save a time....
    if( file.exists(filename) ){
      cat("file exists\n")
      s <- readRDS(filename)
      return(s)
    }else{
      #cat( marker.drug[[input$drug]] )
      cat("file does not exist.\n")
      n.limit <- 1000; #nrow(rma)
      x <- get.MarkerSet()
      
      id <- 1:ncol(x$expr) ## all
      if( input$tissueSite!="all") id <- which(x$clinic$tissue==input$tissueSite)   
      cat( length(id), "samples from ", input$tissueSite, "tissue site are available for a biomarker discovery")  
      # x <- get.PDX.data("paclitaxel"); id <- which(x$clinic$tissue=="breast")

      
      ##########################
      
      
      ### define markerdiscovery set      
      markerset <- x$expr[,id]
      markerset.response <- x$drug.response[id]
      markerset.clinic <- x$clinic[id, ]
      method="pearson"    
      
      ### GSEA..
      cat("We get markerset of size ", ncol(markerset), "\n")
      if( input$BiomarkerDiscoveryMethod!="else"){
        
        
        if( input$BiomarkerDiscoveryMethod == "correlation-rank" ){
          method="spearman"
          s <- biomarker.discovery.correlation(markerset=markerset, markerset.response=markerset.response, method=method)
          
        }else if( input$BiomarkerDiscoveryMethod == "correlation" ){
          method="pearson"   
          s <- biomarker.discovery.correlation(markerset=markerset, markerset.response=markerset.response, method=method)
          
        }else if( input$BiomarkerDiscoveryMethod %in% c("NIRA", "NIRA-rank")) {
          
          method.type <- ifelse(input$BiomarkerDiscoveryMethod=="NIRA", "lm", "lm.rank")
          
          ####  source("./library/library_appendix_2016-08-17.R")
          ####  markers.d2 
          reference.type <- markerset.clinic$tissue[1] #"ovary"
          # args( biomarker.discovery.regression)
          #biomarker.discovery.regression.v03 <- function( markerset, markerset.response, covariate=NULL, 
          #                                                method="lm", interaction=FALSE, sig.level=0.05, plot.check=TRUE ){
          
          s <- biomarker.discovery.regression.v03( markerset=markerset, 
                                                   markerset.response=markerset.response, 
                                                   interaction=TRUE,
                                                   #is.ANOVA=TRUE, 
                                                   sig.level=0.05,
                                                   covariate=data.frame(relevel( factor(markerset.clinic$tissue), ref=reference.type)), 
                                                   method=method.type, plot.check=TRUE)
          attr(s$markerset, "name") <- paste(input$MarkerSet, ncol(s$markerset),sep="_")
          
          
        }else if(input$BiomarkerDiscoveryMethod=="t-test"){
          method="ttest"
          #s <- biomarker.discovery.ttest.v2()
          ### if 0/1 binary responses..
          n.unique <- length(unique(markerset.response))
          if(n.unique==2){
            n.res <- sum(markerset.response==1)
            n.nonres <- sum(markerset.response==0)
            s <- biomarker.discovery.ttest.v02(markerset=markerset, markerset.response=markerset.response, 
                                               n.res = n.res, n.nonres=n.nonres)
          }else{ ## more than 2 clssess
            s <- biomarker.discovery.ttest.v02(markerset=markerset, markerset.response=markerset.response)
            cat("\n Exhaustive search_based_T-test biomarker discovery step has been completed...!\n")
            
          }
          
        }else if(input$BiomarkerDiscoveryMethod=="t-test2"){
          method="ttest"
          #s <- biomarker.discovery.ttest.v2()
          ### if 0/1 binary responses..
          n.unique <- length(unique(markerset.response))
          if(n.unique==2){
            n.res <- sum(markerset.response==1)
            n.nonres <- sum(markerset.response==0)
            s <- biomarker.discovery.ttest.v02(markerset=markerset, markerset.response=markerset.response, 
                                               n.res = n.res, n.nonres=n.nonres)
          }else{ ## more than 2 clssess
            n.res <- sum(markerset.response<0)
            n.nonres <- sum(markerset.response>0)
            s <- biomarker.discovery.ttest.v02(markerset=markerset, markerset.response=markerset.response, 
                                               n.res = n.res, n.nonres=n.nonres)
          }
        }else{
          method="spearman"   
          s <- biomarker.discovery.correlation(markerset=markerset, markerset.response=markerset.response, method=method)
          
        }
        saveRDS(s,  file=filename)
        return(s)
      } ## end if (BiomarkerDiscoveryMethod)
    } ## end: if (file exists)
    # z <- biomarker.discovery.correlation(markerset=markerset, markerset.response=markerset.response, method=method)
  }) ## End of 'get.Biomarker'
  
  
  # m <- biomarker.select(x=mset, p.cutoff=0.05)
  select.BiomarkerSet <- reactive({
    
    ### update Bioamrkers
    ## x <- readRDS( file="./RData/rds/5FU_Novartis_large_intestine_t-test.rds")
    mset <- get.Biomarker()
    
    m <- biomarker.select( x=mset, q.cutoff=0.2) # 478 NCI-60 ANCOVA
    m$FDR.correction <- "controlling FDR of 0.05"
    
    if( length(m$marker) < 30 ){
      cat("Control P-value because FDR-controlled result is so conservative\n")
      m <- biomarker.select( x=mset, p.cutoff=0.05) # 478 NCI-60 ANCOVA 
      m$FDR.correction <- "controlling nominal P-value of 0.05"  
      
      if( length(m$marker) < 5){
        m <- biomarker.select( x=mset, p.cutoff=quantile(mset$marker.info$pvalue, probs=0.01) )
        m$FDR.correction <- "choosing Top 10% genes because no significant genes"  
      }
      
      cat("Control P-value because FDR-controlled result is so conservative")
      cat("# of initial biomarkers=", length(m$marker), "\n")
    }
    # m <- biomarker.select( x=s, q.cutoff=0.2) # 478 NCI-60 ANCOVA
    # m <- biomarker.select( x=s, p.cutoff=0.05) # 478 NCI-60 ANCOVA
    print( head(m$marker) )
    cat("\n# of initial biomarkers=", length(m$marker), "\n")
    
    filename <- paste(input$drug, input$MarkerSet,  input$tissueSite, input$BiomarkerDiscoveryMethod, sep="_")
    filename <- paste("./my-tmp/m_", filename, ".rds", sep="")
    cat("file name of 'm' object for the selected biomarkers is ", filename, "\n")
    #saveRDS( m, file=filename, sep="") )
    m
  }) ## End of 'get.Biomarker'
  
  
  
  ### Initial Biomarkers      
  output$markerinfo <- renderPlot({
    
    if(input$BiomarkerDiscoveryMethod!='else'){
      ## note that 'get.Biomarker() generates plots of test stat vs p-values..
      x <-  get.Biomarker()$marker.info
      # x <-  s$marker.info
      
        xlab <- "Correlation coefficient between gene expression and tumor volume change"
        if( input$BiomarkerDiscoveryMethod %in% c("t-test", "t-test2") ){
          xlab <- "Mean expression difference"
        }
      
        plot( x$estimate, -10*log10(x$pvalue), pch=".", cex=2, xlab=xlab, ylab="-10 Log (P-value)") ## q-value
        grid()
        
        gs <- select.BiomarkerSet()$marker
        
        points(x[gs, "estimate"], -10*log10(x[gs, "pvalue"]), pch=1, col="red")
        #head(x)
    }    
  })
  
  
  
  ## output: Biomarker Discovery
  output$initMarker <- renderText({
    if(input$BiomarkerDiscoveryMethod!='else'){
      x <- select.BiomarkerSet()$marker
      method <- select.BiomarkerSet()$FDR.correction
      # m.select$marker
      head(x)
      cat(length(x), "marker were discovered")
      paste(length(x), "marker were discovered using", method)
    }    
  })
  
  
  # Downloadable csv of selected dataset ----
  output$downInitMarker <- downloadHandler(
    filename = function() {
      #paste(input$dataset, ".csv", sep = "")
      paste("InitMarker", input$drug, input$MarkerSet, input$tissueSite, 
            input$BiomarkerDiscoveryMethod, ".csv", sep="_")
    },
    content = function(file) {
      write.csv( get.Biomarker()$marker.info, file, row.names = TRUE)
    }
  )
  
  
  compute.CCC <- reactive({
    ##i <- 4
    m <- select.BiomarkerSet()
    
    
    target <- m$marker
    cat("\nInput : ", length(target) )
    
    sig.level.ccc <- 0.05
    cor.cutoff <- input$sliderCCC # 0.2
    cat("CCC cor.cutoff=", cor.cutoff, "\n")
    
    method <- "pearson"
    method.1 <- "pearson"
    
    ancillary.genes = get.CCCset(origin=input$origin)$expr[target, ]
    
    ###set1 <- GSE19784[ target, ] # HOVON
    set1 <- ancillary.genes[ target, ] # 
    set2 <- as.matrix(m$markerset[ target, ])
    if( m$method=="ttest" ) set2 <- m$markerset[ target, ]
    
    a <- list()
    a[[1]] <- a.1 <- ccc(set1, set2, probeset=target, cutoff=sig.level.ccc, 
                         cor.cutoff=cor.cutoff, method=method, method.1=method.1)
    
    
    ### gather only significant
    rs.ccc <- a[[1]]$ccc
    #a.1 <- ccc(set1, set2, probeset=target, cutoff=0.05, cor.cutoff=0.2, method="pearson", method.1="pearson")
      length( a.1[[1]] )
      dim( a.1[[2]] )
      dim( a.1[[3]] )
      names(a.1)
      
    a[[1]] <- a.1$probe

    ### if no significant genes or it is too few, increase it by releasing a significance cutoff..
    if( length(a.1$probe) < 10){
      cat("No concordant markers with ccc>0.2 and q<0.2. Use the 2nd option of controlling only Type I error rate < 0.05.\n")
      a.1 <- subset(as.data.frame(a.1$ccc.cor), ccc>0 & ccc.p<0.05)
      rs.ccc <- a.1
      a[[1]] <- rownames(a.1)
    }      
    
    
    rs.ccc <- data.frame( rs.ccc, probe=rownames(rs.ccc)) 
    #cat("A range of absolute CCC is from ", round(min(rs.ccc$ccc),3), "to", round(max(rs.ccc$ccc,3), ".\n"))
    
    filename <- paste("./cccmarker/CCCMarker", input$drug, input$MarkerSet, input$tissueSite, 
                      input$BiomarkerDiscoveryMethod, input$origin, ".rds", sep="_")
    print(filename)
    #saveRDS(rs.ccc, file=filename)    
    
    tmp.result <- list()
    tmp.result$target <- target
    tmp.result$set1 <- set1
    tmp.result$set2 <- set2
    tmp.result$rs.ccc <- rs.ccc
    tmp.result$a.1 <- a.1
    filename <- paste("./cccmarker/CCEC-Summary_", input$drug, input$MarkerSet, input$tissueSite, 
                      input$BiomarkerDiscoveryMethod, input$origin, ".rds", sep="_")
    #saveRDS(tmp.result, file=filename)
    
    rs.ccc
  })    
  
  
  
  get.CCC <- reactive({    
    
    m <- select.BiomarkerSet()
    
    ### Perform CCC analysis:: compute.CCC()
    a <- list()
    a[[1]] <- rownames(compute.CCC())
    inter.all <- unlist(a) #get.intersect(a)
    cat("\n Output : ", length(inter.all) )

        
    # result <- list( marker=matrix(biomarker.select( x=m[[i]], target=inter.all, q.cutoff=1 )$marker, nrow=1), 
    #		marker.info=a.1)	
    result <- matrix(biomarker.select( x=m, target=inter.all, q.cutoff=1 )$marker, nrow=1)
    marker.list.ccc <- result
    
    
    #name.cel <- colnames(nci60.dtp)
    if( m$method=="ttest" ){
      name.cel <- colnames( m$markerset.drug )
    }else name.cel <- colnames( m$markerset )
    
    if( input$useCCC=="No") marker.list.ccc <- select.BiomarkerSet()$marker
    
    marker.list <- list()
    marker.list[[1]] <- matrix(marker.list.ccc, nrow=1) #412
    marker.list
  })  
  
  
  
  output$outputCCCMarker <- renderText({
    s <- ""
    rs.ccc <- compute.CCC()
    s <- paste(s, "A range of absolute CCC is from ", round(min(rs.ccc$ccc),3), "to", round(max(rs.ccc$ccc),3), ".<br/>")
    
    x <- unlist( get.CCC() )
      head(x)
    cat(length(x), "CCC marker were discovered. <br/>")
    ## This folliwing line is for printing a message at the Ouput panel
    paste(s, length(x), "CCC marker were discovered.<br/>") 
  })
  
  
  
  
  
  # Downloadable csv of selected dataset ----
  output$downCCCMarker <- downloadHandler(
    filename = function() {
      #paste(input$dataset, ".csv", sep = "")
      paste("CCE_Marker_", input$drug, input$MarkerSet, input$tissueSite, 
            input$BiomarkerDiscoveryMethod, input$origin, ".csv", sep="_")
    },
    content = function(file) {
      write.csv( compute.CCC(), file, row.names = TRUE)
    }
  )
  
  
  
  
  ## https://stackoverflow.com/questions/19760169/how-to-perform-random-forest-cross-validation-in-r  
  ## Random Forest ==> train and test RFs using only cancer cell lines
  evaluate.RF <- reactive({
    marker <- unlist(get.CCC()) ## Biomarkers
    #if( input$useCCC=="No") marker <- select.BiomarkerSet()$marker
    
    alternative <- "greater"
    
    ## Denote a training set as get.MarkerSet() (the reactivated markerset for the selected drug and cell lines) 
    #x <- get.MarkerSet()
    x <- select.BiomarkerSet() 
    
    trainingset <- x$markerset[marker, ] #x$expr[marker, ]
    trainingset <- t(apply(trainingset, 1, scale)) ## standardization
    trainingset.response <- x$response #x$drug.response
    
    if(input$ranked.response=='rank-transform') trainingset.response <- rank(trainingset.response)
    
    evaluationset <- data.frame(trainingset)
    evaluationset.response <- trainingset.response
    evaluationset.clinic <- data.frame( response.ordinal=trainingset.response, row.names=colnames(trainingset))
    evaluationset.clinic$is.valid <- 1
    
    set.seed(71)
    cat("fitting Random Forest...\n")
    fit.RF <- randomForest( x=t(trainingset[marker, ]), y=trainingset.response)
    # plot(fit.RF)
    
    RF.trained <- list(fit=fit.RF, x=t(trainingset[marker, ]), y=trainingset.response)
    
    filename <- paste("./my-tmp/RF.train", input$drug, input$MarkerSet, input$tissueSite, 
                      input$BiomarkerDiscoveryMethod, input$origin, input$ranked.response, ".rds", sep="_")
    # print(filename)

    # saveRDS(RF.trained, file=filename)
    return( RF.trained  )
  })
  
  
  ## https://www.rdocumentation.org/packages/randomForest/versions/4.6-12/topics/rfcv
  evaluate.RF.CV <- reactive({
    marker <- unlist( get.CCC() ) 
    
    alternative <- "greater"
    ## Denote a training set as get.MarkerSet() (the reactivated markerset for the selected drug and cell lines) 
    x <- get.MarkerSet()
    trainingset <- x$expr[marker, ]
    trainingset <- t(apply(trainingset, 1, scale)) ## standardization
    
    trainingset.response <- x$drug.response
    
    evaluationset <- trainingset
    evaluationset.response <- trainingset.response
    evaluationset.clinic <- data.frame( response.ordinal=trainingset.response, row.names=colnames(trainingset))
    evaluationset.clinic$is.valid <- 1
    
    set.seed(71)
    fit.RF <- rfcv( trainx=t(trainingset[marker, ]), trainy=trainingset.response, cv.fold=5, 
                    scale="log", step=0.5)
    fit.RF
    # plot(fit.RF)
    #return( list(fit=fit.RF, x=t(trainingset[marker, ]), y=trainingset.response)  )
  })
  
  
  ### Random Forest built on the training set
  output$step4RF <- renderPlot({
    
    x <- evaluate.RF() ##
    yhat <- predict(x$fit, newdata=x$x, type="response")
    
    ## RF without Cross-validation    
    par( mfrow=c(1,3))
    plot( x$y, yhat)
    plot( x$fit)
    title("RF without CV")
    
    ## RF with Cross-validation
    x <- evaluate.RF.CV()
    with(x, plot(n.var, error.cv, log="x", type="o", lwd=2))
    title("n.var vs error from CV on trainingset")
    
    #cat('run select.model function', '\n')        
    #cat('buildPredictor', result$p.value, '\n')
  })
  
  
  
  
  
  # Downloadable csv of selected dataset ----
  output$downRandomForest <- downloadHandler(
    filename = function() {
      #paste(input$dataset, ".csv", sep = "")
      paste("RF.train", input$drug, input$MarkerSet, input$tissueSite, 
            input$BiomarkerDiscoveryMethod, input$origin, ".csv", sep="_")
    },
    content = function(file) {
      xx <- evaluate.RF()$fit$importance
      xx <- xx[order(xx[,1], decreasing=T), ]
      gene.symbol <- u133p2.gene.symbol[names(xx), ]
      m <- get.Biomarker()$marker.info
      
      xx <- cbind(xx, gene.symbol, m[names(xx), ])
      
      write.csv(xx, file, row.names = TRUE, na="")
    }
  )
  
  
  
  output$step4RFtrainingset <- renderPlot({
    ## 1. It standardized evaluationset. 
    xx <- evaluate.RF() ##
    
    par(mfrow=c(1,2))
    fit.RF <- xx$fit
    
    ### Variable Importance
    barplot( fit.RF$importance[, 1], las=2) #IncNodePurity, las=2 )
    title("Variable Importance in RF")
    
    ### Performance of RF predictor on PDA training dataset
    yhat <- predict(object=xx$fit, newdata=xx$x)
    ## % change in tumor volume (Best Average Response)
    plot( xx$y, yhat, xlab="Observed % change in tumor volume)", ylab="Predicted Best Response(yhat)")
    grid(); abline(0,1)
    text( xx$y, yhat, labels=1:length(yhat), pos=1, cex=0.8)
    ss <- cor.test(xx$y, yhat)
    title(paste("Prediction on trainingset (PDX data)\n r=", round(ss$estimate,3), "p=", round(ss$p.value,3), sep="") )
  })  
  
  
  output$outputRandomForestVariableImportance <- renderText({
    xx <- evaluate.RF()
    number.positive <- sum(xx$fit$importance[,1]>0)
    cat("the number of genes with positive variable importance is ", number.positive, "\n")
    paste("the number of genes with positive variable importance is ", number.positive) ## for output panel.
  })
 
  
   
  validate.RF <- reactive({
    ## 1. It standardized evaluationset.
    xx <- evaluate.RF() ##
    
    ## get patient test data
    marker <- unlist( get.CCC() )
    
    x <- get.testData()
    
    ## scale a test set for each gene 
    filename.testSet <- paste("./my-tmp/my-tmp-RF", input$testData, x$name, ".RData", sep="") 
    #@ save(xx, x, marker, scaleExpressionData, file=filename.testSet )
    
    cat("\n The file name of the test set is ", filename.testSet, "\n")
    ## load("./my-tmp/my-tmp-RFPancreatic:ICGC_RNAseqPancreatic_ICGC_RNAseq_96.RData")
    ## load("./my-tmp/my-tmp-RFPancreatic:EMEXP2780Pancreatic_EMEXP2780_30.RData")
    
    x.scaled <- scaleExpressionData(x$expr, marker=marker)
    # which( is.na(x.scaled), arr.ind=T)
    evaluationset <- x.scaled
    evaluationset.clinic <- x$clinic
    evaluationset.response <- x$drug.response
    
    fit.RF <- xx$fit
    
    ### Performance of validation on a patient data set
    score <- predict(object=xx$fit, newdata=t(evaluationset), type="response" )
    cat( head(score), "\n")
    
    fileheader <- paste("./my-tmp/", input$MarkerSet,"-", input$drug, "-", input$tissueSite,"-", 
                        input$BiomarkerDiscoveryMethod, "-", input$origin, "-", sep="")
    
    #@ save(xx, x, marker, evaluationset, score, file=paste(fileheader,  "RandomForest-", x$name, ".RData",sep="") )
    
    evaluationset.clinic$score <- score
    
    file=paste(fileheader, "RandomForest", x$name, "clinic.txt", sep="")
    #@ write.table( evaluationset.clinic, file=file, sep="\t")
    
    # sum( is.na(match( marker, rownames(x$expr)) )) # number of missing markers
    # length(marker)  ##..
    
    result <- list( expr=evaluationset, clinic=evaluationset.clinic, drug.response=evaluationset.response, 
                    score=score)
      print(x$name)
    filename.testResult <- paste(fileheader, "RandomForest", x$name, "-v.rds", sep="")
      print(filename.testResult)
    saveRDS(result, file=filename.testResult )
    #saveRDS(result, file="./my-tmp/v.rds")
    cat("\nTest results are saved as ", filename.testResult, "\n")
    
    result    
  })

    
  
  output$step4RFtestset <- renderPlot({
    ##     
    v <- validate.RF()
    
    score <- v$score # v$clinic$score
    evaluationset.response <- v$drug.response
    evaluationset.clinic <- v$clinic
    
    par(mfrow=c(3,4))
    
    ## binary response is available
    if( !is.null(evaluationset.response)){
      n.class <- length( setdiff(unique(evaluationset.response),NA))
      if( n.class<=5){
        plot.strip(a=score, b=evaluationset.response, smaller=T, 
                   data.name=input$testData, 
                   group.names=c("R", "NR"), 
                   cex=1.5, cex.axis=1.5, cex.lab=1.5, col=c("black", "black"), pch=1)
        
        ## two-sided t-test
        if(n.class==2){
          t.pvalue <- t.test( score ~ evaluationset.response)$p.value 
          #mtext(3, at=1, line=-1, text=paste(t.P=round(t.pvalue,3)))
        } 
        
      }else {
        plot(score, evaluationset.response, xlab="prediction score", ylab="continuous drug reponse")
        tmp <- cor.test(score, evaluationset.response, method="spearman")
        r <- round(tmp$estimate, 3)
        r.pvalue <- round(tmp$p.value)
        title(paste("cor=",r, "(", r.pvalue, ")"))        
      }
    }
    
    if( exists("response.ordinal", v$clinic)  ){
        boxplot(v$score ~ v$clinic$response.ordinal)
        stripchart(v$score ~ v$clinic$response.ordinal, vertical=T, add=T, method="jitter")
        
        tmp <- cor.test(v$score, v$clinic$response.ordinal, method="spearman")
        r <- round(tmp$estimate, 3)
        r.pvalue <- round(tmp$p.value)
        #tmp <- cor.test(v$score, v$clinic$response.ordinal)
        #r.pvalue2 <- round(tmp$p.value)
        title(paste("cor=",r, "(", r.pvalue, ")"))        
    }    

    col=c("blue", "green", "orange", "red")
    
    if( exists("OS", evaluationset.clinic) & exists("death", evaluationset.clinic) ){
      require(survival)
      s <- list()
      s$clinic <- evaluationset.clinic
      
      ### OS   
      s$clinic$time <- s$clinic$OS
      s$clinic$censor <- s$clinic$death
      
      time.end <- 60
      s$clinic$time <- pmin(s$clinic$time, time.end)
      s$clinic$censor <- pmin(as.numeric(s$clinic$time<time.end), s$clinic$censor)
      s$clinic$time <- pmax(s$clinic$time, 0.01) ## to make non-zero survival time
      
      ## OS
      #par( mfrow=c(2,2))
      for(n.group in 2:3){
        s$clinic$group <-cut(s$clinic$score, quantile(s$clinic$score, probs=seq(0, 1, 1/n.group) ), include.lowest = T )
        sdata <- s$clinic
        plot.surv(formula="Surv(time, censor)~group", data=sdata, 
                  title=paste("OS:n=", nrow(s$clinic),sep=""), 
                  xlab="OS",
                  lwd=4, legend.pos ="bottomright", lty=1, col=col, cex.axis=2, cex.lab=2, mark.time=T )
      }

    }
    
    if( exists("PFS", evaluationset.clinic) & exists("progression", evaluationset.clinic) ){
      require(survival)
      s <- list()
      s$clinic <- evaluationset.clinic
      
      col=c("blue", "green", "orange", "red")
      
      ### PFS   
      s$clinic$time <- s$clinic$PFS
      s$clinic$censor <- s$clinic$progression
      
      time.end <- 60
      s$clinic$time <- pmin(s$clinic$time, time.end)
      s$clinic$censor <- pmin(as.numeric(s$clinic$time<time.end), s$clinic$censor)
      s$clinic$time <- pmax(s$clinic$time, 0.01) ## to make non-zero survival time
      
      ## PFS
      #par( mfrow=c(2,2))
      for(n.group in 2:3){
        s$clinic$group <-cut(s$clinic$score, quantile(s$clinic$score, probs=seq(0, 1, 1/n.group) ), include.lowest = T )
        sdata <- s$clinic
        plot.surv(formula="Surv(time, censor)~group", data=sdata, 
                  title=paste("PFS:n=", nrow(s$clinic),sep=""), lwd=4, legend.pos ="bottomright", 
                  xlab="PFS",
                  lty=1, col=col, cex.axis=2, cex.lab=2, mark.time=T )
      }     
    }
    
  }) # end of 'output$step4 '
  
  
  output$downRFvalid <- downloadHandler(
    filename = function() {
      #paste(input$dataset, ".csv", sep = "")
      paste("RF.validation", input$drug, input$MarkerSet, input$tissueSite, 
            input$BiomarkerDiscoveryMethod, input$origin, "_", input$testData, ".csv", sep="_")
    },
    content = function(file) {
      xx <- validate.RF()$clinic
      write.csv(xx, file, row.names = TRUE)
    }
  )
  
  
  # Downloadable csv of selected dataset ----
  output$downRFprediction <- downloadHandler(
    filename = function() {
      #paste(input$dataset, ".csv", sep = "")
      paste("RF.predictor-", input$drug, input$MarkerSet, input$tissueSite, 
            input$BiomarkerDiscoveryMethod, input$origin, ".csv", sep="_")
    },
    content = function(file) {
      xx <- evaluate.RF()$fit$importance
      xx <- xx[order(xx[,1], decreasing=T), ]
      gene.symbol <- u133p2.gene.symbol[names(xx), ]
      m <- get.Biomarker()$marker.info
      
      xx <- cbind(xx, gene.symbol, m[names(xx), ])
      
      write.csv(xx, file, row.names = TRUE, na="")
    }
  )
  
  
  
  
  
  
  
  ### training model on PDXs and test models on a human cancer cohort
  evaluate.Predictors <- reactive({
    
    marker <- unlist( get.CCC() ) ## 
    cat(marker, "\n")
    
    alternative <- "greater"
    p <- 100
    n.pca <- 3
    
    ## Denote a training set as get.MarkerSet() (the reactivated markerset for the selected drug and cell lines) 
    ## this object bears the marker discovery set    
    x <- select.BiomarkerSet() ## markerset, markerset.info, response
    trainingset <- x$markerset[marker, ]
    trainingset.response <- x$response
    trainingset.name <- paste("PDX", ncol(trainingset), sep="-")
    trainingset.name
    
    x <- get.testData()
    
    evaluationset <- x$expr[marker, ]
    evaluationset.response <- x$drug.response
    evaluationset.clinic <- x$clinic
    
    cat('run select.model function', '\n')        
    result <- select.model( marker=marker, 
                            evaluationset=evaluationset, evaluationset.response=evaluationset.response, 
                            evaluationset.clinic=evaluationset.clinic,
                            trainingset=trainingset, trainingset.response=trainingset.response, 
                            n.pca=n.pca,
                            model.type="regression", evaluation.type="Wilcoxon", 
                            p=p, alternative=alternative )
    cat('buildPredictor', result$p.value, '\n')
    result
  })
  
  
  ###### PCA regression and its evaluation results on a test data set
  output$step4 <- renderPlot({

    xx <- evaluate.Predictors() ## model evaluation results..
    
    rs <- data.frame(num=1:length(xx$p.value), AUC=xx$estimate, pvalue=xx$p.value, aic=xx$aic)
    plot(x=rs$num, y=rs$pvalue, type="b", ylim=c(0,1), 
         xlab="Model Complexity", ylab="AUC P-value")
    abline( h=0.025, lty=2)
    points(x=rs$num, y=rs$AUC, type="b", col="red")
    abline( h=0.05, lty=2)
    #    plot(x=rs$num, y=rs$estimate, type="b", ylim=c(0,1), 
    #         xlab="Model Complexity", ylab="AUC")
    #    abline( h=0.5, lty=2)
    title( input$testData )
    
  }) # end of 'output$step4 '
  
  
  
  ### output: use render'Print' to print out a dataframe.  
  output$step4info2 <- renderPrint({
    coord.mouse <- input$bestPCregression_click
    #paste0("x=", input$bestModel_click$x, "\ny=", input$bestModel_click$y)
    
    xx <- evaluate.Predictors() ## output of 'select.model()'
    
    rs <- data.frame(num=1:length(xx$p.value), pvalue=xx$p.value, aic=xx$aic)
    nearPoints(df=rs, coordinfo=coord.mouse, xvar = "num", yvar = "pvalue")
  })
  
  
  
  
  get.oneRF <- reactive({
    xx <- evaluate.RF()$fit ##
    #coord.mouse <- input$bestModel_click
    n.marker <- 40 #round( coord.mouse$x)
    var.importance <- names(xx$importance[ order(xx$importance, decreasing=TRUE), ])
    marker <- var.importance[1:n.marker]
    
    ## Denote a training set as get.MarkerSet() (the reactivated markerset for the selected drug and cell lines) 
    x <- get.MarkerSet()
    trainingset <- x$expr[marker, ]
    trainingset <- t(apply(trainingset, 1, scale)) ## standardization
    trainingset.response <- x$drug.response
    
    fit.RF <- randomForest( x=t(trainingset[marker, ]), y=trainingset.response)
    fit.RF
  })
  
  
  ### Output: stripchart and ROC curve of Test Set
  output$step4testData <- renderPlot({
    x <- get.onePredictor()
    par(mfrow=c(2,2))
    score <- x$coxen.score
    
    evaluationset.respone <- x$evaluationset.response
    evaluationset.clinic <- x$evaluationset.clinic
    
    ## binary response is available
    if( !is.null(evaluationset.response)){
      n.class <- length(unique(evaluationset.response))
      if( n.class<=5){
        plot.strip(a=score, b=evaluationset.response)
      }else {
        plot(score, evaluationset.response, xlab="prediction score", ylab="continuous drug reponse")
        r <- round(cor.test(score, evaluationset.response)$estimate, 3)
        title(paste("cor=",r))        
      }
    }
    
    if( exists("OS", evaluationset.clinic) & exists("death", evaluationset.clinic) ){
      require(survival)
      sfit <- survfit( Surv(OS, death) ~ I(score > median(score)), data=evaluationset.clinic)
      sdiff <- survfit( Surv(OS, death) ~ I(score > median(score)), data=evaluationset.clinic)
      plot(sfit, col=c("green", "red"))
      legend("topright", legend=c("score<median", "score>median"), col=c("green", "red"), lty=1)
      title("Overall Survival")
      
    }
    
    if( exists("PFS", evaluationset.clinic) & exists("progression", evaluationset.clinic) ){
      require(survival)
      sfit <- survfit( Surv(PFS, progression) ~ I(score > median(score)), data=evaluationset.clinic)
      sdiff <- survfit( Surv(PFS, progression) ~ I(score > median(score)), data=evaluationset.clinic)
      plot(sfit, col=c("green", "red"))
      legend("topright", legend=c("score<median", "score>median"), col=c("green", "red"), lty=1)
      title("Progression Free Survival")
    }
    
  })
  
  
  output$modelEquation <- renderPlot({
    x <- get.onePredictor()$final.coefficient
    #x$final.coefficient
    id.order <- order(x$final.reg.coef)
    barplot( x$final.reg.coef[id.order], names=as.character(x$biomarker.name[id.order]), las=2)
  }) # end of 'output$step4 '
  
  
  output$downloadPredictor <- downloadHandler(
    
    filename <- function() {
      paste("./output/", input$MarkerSet, input$drug, nrow(get.onePredictor()$final.coefficient), "-genes-predictor.csv", sep = "_")
    },
    content <-  function(file) {
      write.csv(get.onePredictor()$final.coefficient, file, row.names = FALSE)
    }
  )
  
  

  output$step5RFvalidset <- renderPlot({
    xx <- evaluate.RF() ##
    
    ## get patient test data
    x <- get.validData()
    marker <- get.3way()
    x$expr <- x$expr[marker, ] 
    
    evaluationset.clinic <- x$clinic
    ## scale a test set for each gene 
    evaluationset <- t( apply(x$expr, 1, scale))
    evaluationset.response <- x$drug.response
    
    par(mfrow=c(1,3))
    score <- predict(object=xx$fit, newdata=t(evaluationset), type="response" )
    plot.strip(a=-1*score, b=evaluationset.response)
    
    ### the best RF model    
    ##par(mfrow=c(1,3))
    fit.RF <- get.oneRF()
    marker <- rownames(fit.RF$importance)
    score <- predict(object=fit.RF, newdata=t(evaluationset[marker, ]), type="response" )
    plot.strip(a=-1*score, b=evaluationset.response)
  }) # end of 'output$step4 '
  
  
  
  
  output$useData <- renderTable({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    inFile <- input$userFile
    
    if (is.null(userFile))
      return(NULL)
    
    read.csv(inFile$datapath, header = input$header)
  })
  
  
  observeEvent(input$tissueSite, {
    x <- input$tissueSite
    if( x=="lung"){
      updateSelectInput(session, "origin", selected = "NSCLC")
    }else if( x %in% c("pancreatic", "pancreas", "PCC") ){
      updateSelectInput(session, "origin", selected = "Pancreatic")
    }else if( x %in% c("breast") ){
      updateSelectInput(session, "origin", selected = "Breast")
    }
  }) ## end of observEvent
  
  
  observeEvent(input$origin, {
    x <- input$origin
    if( input$origin=="Ovarian"){
      updateSelectInput(session, "testData", selected = "OV:TCGA-448")
      updateSelectInput(session, "validData", selected = "Hess-133")
    }else if( input$origin=="Breast"){
      updateSelectInput(session, "testData", selected = "Hess-133")
      updateSelectInput(session, "validData", selected = "OV:TCGA-448")
    }else if( input$origin %in% c('Colon', 'colorectal', 'CRC') ){
      updateSelectInput(session, "testData", selected = "CRC:Khambata")
      updateSelectInput(session, "validData", selected = "CRC:Khambata")
    }else if( input$origin %in% c("lung", "NSCLC") ){
      updateSelectInput(session, "testData", selected = "NSCLC:GSE31625CellLines")
      updateSelectInput(session, "validData", selected = "NSCLC:GSE31625CellLines")
    }else if( input$origin %in% c("pancreatic", "PCC") ){
      updateSelectInput(session, "testData", selected = "Pancreatic:EMEXP2780")
      updateSelectInput(session, "validData", selected = "Pancreatic:EMEXP2780")
    }
  }) ## end of observEvent
  
  
  
  observeEvent(input$BiomarkerDiscoveryMethod, {
    x <- input$sliderCCC
    #if( x=="lung"){
      updateSliderInput(session, "sliderCCC", value = 0.2)
    #}
    
  }) ## end of observEvent
  
  
  
  ### add UI
  observeEvent(input$add, {
    insertUI(
      selector = "#add",
      where = "afterEnd",
      ui = textInput(paste0("txt", input$add),
                     "Insert some text")
    )
  })
  
  ## remove UI
  observeEvent(input$rmv, {
    removeUI(
      selector = "div:has(> #txt)"
    )
  })
  
}) ## end of 'shinyr'


