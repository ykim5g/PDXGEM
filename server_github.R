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
source("./library/library-2016-01-25_new.R")
source("./library/library_appendix_2017-08-10.R")

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







## Illumina HT-12 v3 platform.
## It has missing data
## EuroPDX-Breast
## Illumina data was matched to Affy-u133a
##### reference: COXEN/PDX/CELL-BREAST/Drug Screening plots 
get.PDX.data.Bruna <- function( compound.name ){
  
  expression <- readRDS(file="./RData/EuroPDX_Breast/EuroPDX_Breast_Expression_Matched2AffyU133A.rds")
  
  PDX.array <- expression #[,-1]
  
  drugSensitivity <- readRDS(file="./RData/EuroPDX_Breast/EuroPDX_Breast_drugSensitivity.rds")
  drugSensitivity$BestResponse <- drugSensitivity$AUC
  drugSensitivity$tissue <- "breast"
  
  PDX.clinic <- drugSensitivity %>% dplyr::filter(Drug==compound.name)
  rownames(PDX.clinic) <- PDX.clinic$Model
  
  id.valid <- intersect( PDX.clinic$Model, colnames(PDX.array))
  probe.init <- rownames(expression) ; #rm( rma)
  #probe.init <- 1:nrow(PDX.array)
  
  PDX.expr.drug <- PDX.array[probe.init, id.valid]
  PDX.clinic.drug <- PDX.clinic[id.valid, ]
  ### % change in tumor volume (Best Average Response)
  PDX.drug.response <- PDX.clinic.drug$BestResponse
  names(PDX.drug.response) <- colnames(PDX.expr.drug)
  dim(PDX.expr.drug) # 22215 x 21 (54675 x 21)
  dim(PDX.clinic) # 117 x 71
  head( PDX.clinic)
  
  return( list( expr=PDX.expr.drug, drug.response=PDX.drug.response, clinic=PDX.clinic.drug ) )
}





build.Predictor <- function(marker, evaluation.type="Wilcoxon", n.pca=3, p=100, alternative="greater"){
  result <- list()
  #p <- 50 # length of biomarkers
  #  for( i in 1:length(marker.list) ){
  #    marker <- marker.list[[i]]
  result[[i]] <- select.model( marker=marker, 
                               evaluationset=evaluationset, evaluationset.response=evaluationset.response, 
                               evaluationset.clinic=evaluationset.clinic,
                               trainingset=trainingset, trainingset.response=trainingset.response, n.pca=n.pca,
                               model.type="regression", evaluation.type=evaluation.type, p=p, alternative=alternative )
  
} ## end of 'build.Predictor'




get.CCCset <- function(origin="Breast"){
  filename=paste("./RData/CCCset/", origin, "-CCC.rds", sep="")
  cat(filename, "\n")
  x <- readRDS( filename )
  x <- list(expr=as.matrix(x))
  return (x)
}




get.dataSet <- function(data.name="Hess-133", compound.name="paclitaxel"){
  evaluationset.clinic <- NULL
  evaluationset <- NULL
  evaluationset.response <- NULL
  
  if( data.name=="Hess-133"){
    if( tolower(compound.name) %in% c("paclitaxel", "5fU", "fivefu", "5-fluorouracil", "adriamycin", "cyclophosphamide") ) id <- 1:ncol(hess133)
    evaluationset <- hess133
    evaluationset.response <- hess133.response
    evaluationset.clinic <- hess133.clinic
  }
  
  
  if( data.name=="TFAC-100"){
    if( tolower(compound.name) %in% c("paclitaxel", "5fU", "fivefu", "5-fluorouracil", "adriamycin", "cyclophosphamide") ) id <- 1:ncol(hess100)
    evaluationset <- hess100
    evaluationset.response <- hess100.response
    evaluationset.clinic <- hess100.clinic
  }

  if( data.name=="BR:GSE20194(MAQCII)"){
    GSE20194.QC <- readRDS( file="./Combo Chemo Data/breast/GSE20194/GSE20194.QC.rds")
    GSE20194.QC.clinic <- readRDS(file="./Combo Chemo Data/breast/GSE20194/GSE20194.QC.clinic.rds")
    id <- which( GSE20194.QC.clinic$Treatment.Code!="") ## n=259
    
    evaluationset <- GSE20194.QC[,id]
    evaluationset.response <- GSE20194.QC.clinic$response[id]
    evaluationset.clinic <- GSE20194.QC.clinic[id, ]    
  }
  
  if( data.name=="BR:GSE41998:Horak(all)"){
    ### Breast
    load( file="./Combo Chemo Data/breast/GSE41998.total.RData")
    
    id <- 1:ncol(christine279)
    
    evaluationset <- christine279[,id]
    evaluationset.response <- christine279.clinic[id, ]$response
    evaluationset.clinic <- data.frame( response.binary=evaluationset.response) #
  }
  
    
  if( data.name=="BR:GSE41998:Horak(paclitaxel)"){
    ### Breast
    load( file="./Combo Chemo Data/breast/GSE41998.total.RData")
    
    id <- 1:ncol(christine279)
    #if( tolower(compound.name)=="paclitaxel")
    id <- grep( "treatment arm: Paclitaxel", as.character(christine279.clinic$treatment.arm) )
    #if( tolower(compound.name)=="ixabepilone")
    #  id <- grep( "treatment arm: Ixabepilone", as.character(christine279.clinic$treatment.arm) )
    
    evaluationset <- christine279[,id]
    evaluationset.response <- christine279.clinic[id, ]$response
    evaluationset.clinic <- data.frame( response.binary=evaluationset.response) #
  }
  
  if( data.name=="BR:GSE41998:Horak(ixabepilone)"){
    ### Breast
    load( file="./Combo Chemo Data/breast/GSE41998.total.RData")
    
    id <- 1:ncol(christine279)
    #if( tolower(compound.name)=="paclitaxel")
    #  id <- grep( "treatment arm: Paclitaxel", as.character(christine279.clinic$treatment.arm) )
    #if( tolower(compound.name)=="ixabepilone")
    id <- grep( "treatment arm: Ixabepilone", as.character(christine279.clinic$treatment.arm) )
    
    evaluationset <- christine279[,id]
    evaluationset.response <- christine279.clinic[id, ]$response
    evaluationset.clinic <- data.frame( response.binary=evaluationset.response) #
  }
  
  ### GSE20271
  if( data.name=="BR:GSE20271:Tabchy(all)"){
    ### Breast
    load("./Combo Chemo Data/breast/tabchy178.total.RData") # Lymphoma : Hummel

    id <- 1:nrow(popovici178.clinic) #grep( "TFAC", popovici178.clinic$treatment )
    
    evaluationset <- popovici178[,id]
    evaluationset.response <- popovici178.clinic$response[id]
    evaluationset.clinic <- data.frame( response.binary=evaluationset.response) #
  }    
  
  ### GSE20271
  if( data.name=="BR:GSE20271:Tabchy(TFAC)"){
    ### Breast
    load("./Combo Chemo Data/breast/tabchy178.total.RData") # Lymphoma : Hummel

    id <- grep( "TFAC", popovici178.clinic$treatment )
    
    evaluationset <- popovici178[,id]
    evaluationset.response <- popovici178.clinic$response[id]
    evaluationset.clinic <- data.frame( response.binary=evaluationset.response) #
  }    

  ### GSE20271
  if( data.name=="BR:GSE20271:Tabchy(FAC)"){
    ### Breast
    load("./Combo Chemo Data/breast/tabchy178.total.RData") # Lymphoma : Hummel
    #if( tolower(compound.name) %in% c("5fu", "fivefu", "5-fluorouracil", "adriamycin", "cyclophosphamide") ) 
    id <- which(popovici178.clinic$treatment=="FAC")
    evaluationset <- popovici178[,id]
    evaluationset.response <- popovici178.clinic$response[id]
    evaluationset.clinic <- data.frame( response.binary=evaluationset.response) #
  }    
  
  
  ### Breast Cancer Data: Miyake115: u133p2: GSE32646
  ## TFEC treatment: paclitaxel (80 mg/m2) weekly for 12 cycles followed by 5-FU(500 mg/m2), epirubicin (75 mg/m2) and cyclophosphamide(500 mg/m2) every 3 weeks for four cycles (paclitaxel followedby  5-fluorouracil/epirubicin/cyclophosphamide  [P-FEC])
  if( data.name=="BR:GSE32646:Miyake"){
    ### Breast: Miyake115 : GSE32646
    ### neoadjuvant paclitaxel followed by 5-fluorouracil/epirubicin/cyclophosphamide (P-FEC)
    evaluationset <- readRDS(file="./Combo Chemo Data/breast/GSE32646_Miyake/GSE32646.rds")
    evaluationset.clinic <- readRDS(file="./Combo Chemo Data/breast/GSE32646_Miyake/GSE32646.clinic.rds") #
    # dim( evaluationset) ## 54613 x 115
    # cbind( colnames(evaluationset), as.character(evaluationset.clinic$geo_accession) )
    evaluationset.response <- evaluationset.clinic$response
  }    
  
  ### Breast: Hatzis Christos488:GSE25055: taxane + anthracycline: u133A
  #  Discovery cohort for genomic predictor of response and survival following neoadjuvant taxane-anthracycline chemotherapy  in breast cancer
  if( data.name=="BR:GSE25055:Hatzis"){
    evaluationset <- readRDS(file="./Combo Chemo Data/breast/GSE25055_Christos488_Discovery cohort/Christos488.rds") #
    evaluationset.clinic <- readRDS(file="./Combo Chemo Data/breast/GSE25055_Christos488_Discovery cohort/Christos488.clinic.rds")
    ## 2018-11-02: the observed survial time is not overall survival but Disease-relapse free survival.
    evaluationset.clinic$PFS <- evaluationset.clinic$OS
    evaluationset.clinic$progression <- evaluationset.clinic$death
    evaluationset.clinic$OS <- NULL
    evaluationset.clinic$death <- NULL
    
    evaluationset.response <- evaluationset.clinic$response
  }    
  
  
  ## USO has two arms: TFAC +/- Transtruzumab
  if( data.name=="BR:GSE42822:USO(FAC+TX+Trastuzumab)"){
    ### Breast: USO (Herceptin)
    load( file="./Combo Chemo Data/breast/GSE42822.total.RData")
    ### Herceptin
    ##id <- grep( "FEC/TX+H", USO.clinic$treatment.type )
    #sort(unique(PDX.curve$Treatment))
    table(USO.clinic$treatment.type)
    #if( tolower(compound.name) %in% c("paclitaxel", "5fu", "fivefu", "5-fluorouracil", 
    #                                  "cyclophosphamide", "epirubicin"))  
    #id <- grep("FEC/TX", USO.clinic$treatment.type)
    #if( tolower(compound.name) %in% c("herceptin", "trastuzumab")){
    id <- which( USO.clinic$treatment.type=="FEC/TX+H")
    #}
    evaluationset <- USO[,id]
    evaluationset.response <- USO.clinic$response[id]
    evaluationset.clinic <- USO.clinic[id, ]
  }
  
  if( data.name=="BR:GSE42822:USO(FAC+TX)"){
    ### Breast: USO (Herceptin)
    load( file="./Combo Chemo Data/breast/GSE42822.total.RData")
    ### Herceptin
    ##id <- grep( "FEC/TX+H", USO.clinic$treatment.type )
    #sort(unique(PDX.curve$Treatment))
    table(USO.clinic$treatment.type)
    #if( tolower(compound.name) %in% c("paclitaxel", "5fu", "fivefu", "5-fluorouracil", 
    #                                  "cyclophosphamide", "epirubicin"))  
    
    id <- grep("FEC/TX", USO.clinic$treatment.type)
    
    #if( tolower(compound.name) %in% c("herceptin", "trastuzumab")){
    #  id <- which( USO.clinic$treatment.type=="FEC/TX+H")
    #}
    evaluationset <- USO[,id]
    evaluationset.response <- USO.clinic$response[id]
    evaluationset.clinic <- USO.clinic[id, ]
  }  
  
  
  
  if( data.name=="BR:GSE22266:I-SPY(ACT+Trastuzumab):Agilent"){
    ### Breast: USO (Herceptin)
    x <- readRDS(file="./Combo Chemo Data/breast/GSE22266_ISPY_Agilent/GSE22226.u133p2.rds")
    x.clinic <- readRDS(file="./Combo Chemo Data/breast/GSE22266_ISPY_Agilent/GSE22226.u133p2.clinic.rds")
    ## table(x.clinic$treatment)
    #ac only           ac+t ac+t+herceptin     ac+t+other 
    #4            112             11              2 
    ### Herceptin
    ## http://www.robelle.com/smugbook/regexpr.html
    if( tolower(compound.name) %in% c("paclitaxel", "cytoxan", "cyclophosphamide", "adriamycin"))  id <- grep("ac[+]t", x.clinic$treatment)
    if( tolower(compound.name) %in% c("herceptin", "trastuzumab")){
      id <- which( x.clinic$treatment=="ac+t+herceptin")
    }
    evaluationset <- x[,id]
    evaluationset.response <- x.clinic$response[id]
    evaluationset.clinic <- x.clinic[id, ]
  }  
  
  ## Neoadjuvant docetaxel-Capecitabine +/- trastuzumab
  if( data.name=="BR:GSE22358:Agilent:Trastuzumab"){
    x <- readRDS(file="./Combo Chemo Data/breast/GSE22358/GSE22358.u133p2.RDS")
    x.clinic <- readRDS(file="./Combo Chemo Data/breast/GSE22358/GSE22358.u133p2.clinic.RDS")
      table( x.clinic$NEOADJUVANT.CHEMOTHERAPY )
      #x.clinic$NEOADJUVANT.CHEMOTHERAPY==(docetaxel-Capecitabine)+ Trastuzumab
    id <- which(x.clinic$NEOADJUVANT.CHEMOTHERAPY=="(docetaxel-Capecitabine)+ Trastuzumab")
    
    evaluationset <- x[,id]
    evaluationset.response <- x.clinic$response.binary[id]
    evaluationset.clinic <- x.clinic[id, ]
  }
  
  ## Neoadjuvant docetaxel-Capecitabine +/- trastuzumab
  if( data.name=="BR:GSE22358:Agilent:all"){
    x <- readRDS(file="./Combo Chemo Data/breast/GSE22358/GSE22358.u133p2.RDS")
    x.clinic <- readRDS(file="./Combo Chemo Data/breast/GSE22358/GSE22358.u133p2.clinic.RDS")
    table( x.clinic$NEOADJUVANT.CHEMOTHERAPY )
    #x.clinic$NEOADJUVANT.CHEMOTHERAPY==(docetaxel-Capecitabine)+ Trastuzumab
    id <- 1:nrow(x.clinic) ##which(x.clinic$NEOADJUVANT.CHEMOTHERAPY=="docetaxel-Capecitabine)+ Trastuzumab")
    
    evaluationset <- x[,id]
    evaluationset.response <- x.clinic$response.binary[id]
    evaluationset.clinic <- x.clinic[id, ]
  }
  
 if( data.name=="CRC:Khambata(Cetuximab)"){
    ### colorectal: Khambata-Ford
    load( file="./Combo Chemo Data/colorectal/GSE5851_Khambata/GSE5851.total.RData")
    
    GSE5851.clinic$response.binary <- car::recode(GSE5851.clinic$Best.Clinical.Response.Assessment, 
                                             recode="c('CR','PR')=1;c('SD', 'PD')=0;else=NA")
    require(car)
    GSE5851.clinic$kras <- car::recode(GSE5851.clinic$KRAS.Mutation, 
                                  recode="c('WT')='WT';c('na','')=NA;else='Mutant'")
    GSE5851.clinic$kras2 <- car::recode(GSE5851.clinic$KRAS.Mutation, 
                                   recode="c('WT','')='WT';c('na')=NA;else='Mutant'")
    #id <- which( GSE5851.clinic$KRAS.Mutation=="WT" ) #
    id <- 1:nrow(GSE5851.clinic)
    evaluationset.clinic <- GSE5851.clinic[id, ]
    evaluationset <- GSE5851[,id]
    evaluationset.response <- GSE5851.clinic[id, ]$response.binary ##Best.Clinical.Response.Assessment
  }
  
  if( data.name=="CRC:Khambata:KRASwild(Cetuximab)"){
    ### colorectal: Khambata-Ford
    load( file="./Combo Chemo Data/tal/GSE5851_Khambata/GSE5851.total.RData")
    GSE5851.clinic$response.binary <- car::recode(GSE5851.clinic$Best.Clinical.Response.Assessment, 
                                                  recode="c('CR','PR')=1;c('SD', 'PD')=0;else=NA")
      require(car)
    GSE5851.clinic$kras <- car::recode(GSE5851.clinic$KRAS.Mutation, 
                                  recode="c('WT')='WT';c('na','')=NA;else='Mutant'")
    GSE5851.clinic$kras2 <- car::recode(GSE5851.clinic$KRAS.Mutation, 
                                   recode="c('WT','')='WT';c('na')=NA;else='Mutant'")
    #id <- which( GSE5851.clinic$KRAS.Mutation=="WT" ) #
    id <- which( GSE5851.clinic$KRAS.Mutation=="WT")
    evaluationset.clinic <- GSE5851.clinic[id, ]
    evaluationset <- GSE5851[,id]
    evaluationset.response <- GSE5851.clinic[id, ]$response.binary ##Best.Clinical.Response.Assessment
  }
  
  if( data.name=="CRC:Khambata:KRASmutation(Cetuximab)"){
    ### colorectal: Khambata-Ford
    load( file="./Combo Chemo Data/tal/GSE5851_Khambata/GSE5851.total.RData")
    GSE5851.clinic$response.binary <- car::recode(GSE5851.clinic$Best.Clinical.Response.Assessment, 
                                                  recode="c('CR','PR')=1;c('SD', 'PD')=0;else=NA")
    #GSE5851.clinic$response.ordinal <- car::recode(GSE5851.clinic$Best.Clinical.Response.Assessment, 
    #                                               recode="'CR'=1;'PR'=2;'SD'=3;'PD'=4;else=NA")
    
    #GSE5851.clinic$response.ordinal <- car::recode(GSE5851.clinic$Best.Clinical.Response.Assessment, 
    #                                               recode="'CR'=1;'PR'=2;'SD'=3;'PD'=4;'UTD'=5;else=NA")
    require(car)
    GSE5851.clinic$kras <- car::recode(GSE5851.clinic$KRAS.Mutation, 
                                  recode="c('WT')='WT';c('na','')=NA;else='Mutant'")
    GSE5851.clinic$kras2 <- car::recode(GSE5851.clinic$KRAS.Mutation, 
                                   recode="c('WT','')='WT';c('na')=NA;else='Mutant'")
    
    #id <- which( GSE5851.clinic$KRAS.Mutation=="WT" ) #
    id <- which( GSE5851.clinic$kras=="Mutant")
      cat("Selected ID for a testset are ", id, "\n")
    evaluationset.clinic <- GSE5851.clinic[id, ]
    evaluationset <- GSE5851[,id]
    ### In this case, there was no responsive patients..
    evaluationset.response <- NULL #GSE5851.clinic[id, ]$response.binary ##Best.Clinical.Response.Assessment
  }

  
  if( length(grep("CRC_GSE62322_FOLFIRI", data.name))>=1 ) {
    ### colorectal
    x <- readRDS(file="./Combo Chemo Data/colorectal/GSE62322/GSE62322.tumor.rds")
    x.clinic <-  readRDS( file="./Combo Chemo Data/colorectal/GSE62322/GSE62322.tumor.clinic.rds")
  
    ## 
    if( data.name=="CRC_GSE62322_FOLFIRI_Primary"){    
      id <-  which(x.clinic$organism == "Primary Tumor") ##grep( "treatment arm: Paclitaxel", as.character(christine279.clinic$treatment.arm) )
    }else if( data.name=="CRC_GSE62322_FOLFIRI_Metastasis"){    
      id <-  which(x.clinic$organism == "Liver Metastasis") ##grep( "treatment arm: Paclitaxel", as.character(christine279.clinic$treatment.arm) )
    }else{
      id <- 1:nrow(x.clinic)
    }
    
    #id <-  which(x.clinic$organism == "Primary Tumor") ##grep( "treatment arm: Paclitaxel", as.character(christine279.clinic$treatment.arm) )
    
    evaluationset <- x[,id]
    evaluationset.clinic <- x.clinic[id, ]
    evaluationset.response <- evaluationset.clinic$response
    evaluationset.name=paste("CRC:GSE62322-", ncol(evaluationset))
  }
  
  ### Del Rio et al: Refined dataset of gene expression profiles of 21 metastatic tmor tissues
  if( length(grep("CRC_GSE62080_FOLFIRI_Metastasis", data.name))>=1 ) {
    ### colorectal
    x <- readRDS(file="./Combo Chemo Data/colorectal/GSE62080/GSE62080.rds")
    x.clinic <-  readRDS( file="./Combo Chemo Data/colorectal/GSE62080/GSE62080.clinic.rds")
      id <- 1:nrow(x.clinic)
    evaluationset <- x[,id]
    evaluationset.clinic <- x.clinic[id, ]
    evaluationset.response <- evaluationset.clinic$response
    evaluationset.name=paste("CRC:GSE62080-", ncol(evaluationset))
  }
  
  
  
  if( length(grep("CRC_GSE54483_FOLFIRI", data.name))>=1 ) {
    ### colorectal
    x <- readRDS(file="./Combo Chemo Data/colorectal/GSE54483/GSE54483.rds")
    x.clinic <-  readRDS( file="./Combo Chemo Data/colorectal/GSE54483/GSE54483.clinic.rds")
    
    id <- 1:nrow(x.clinic)

    evaluationset <- x[,id]
    evaluationset.clinic <- x.clinic[id, ]
    evaluationset.response <- evaluationset.clinic$response.binary
    evaluationset.name=paste("CRC:GSE54483-", ncol(evaluationset))
  }

  
  #### GSE14095:Watanabe et al: leucovorin, fluorouracil and oxaliplatin(n=189) (u133p2)
  ### Response:: Metastatsis or Not..
  if( data.name=="CRC_GSE14095_5FU+Leucovorin+Oxaliplatin"){
    ### colorectal
    x <- readRDS(file="./Combo Chemo Data/colorectal/GSE14095/GSE14095.rds")
    x.clinic <-  readRDS( file="./Combo Chemo Data/colorectal/GSE14095/GSE14095.clinic.rds")
      table( x.clinic$response, x.clinic$DM )
      
    ## table( unlist( lapply( strsplit(x.clinic$description, ","), function(x) x[1] ) ) )
    id <- 1:ncol(x) 
    evaluationset <- x[,id]
    evaluationset.clinic <- x.clinic[id, ]
    evaluationset.response <- evaluationset.clinic$response
    evaluationset.name=paste("CRC:GSE14095-", ncol(evaluationset))
  }
  
 
  
  
  
  #### GSE39582:: 5FU, FOLFIRI(n=xxx) (u133p2)
  # Adjuvant: Folinic acid, fluorouracil and oxaliplatin (FOLFOX)
  # Irinotecan with fluorouracil (5FU) and folinic acid (FOLFIRI)
  # FU+FOL
  
  if( length(grep("GSE39582", data.name)) >=1 ){
    ### colorectal
    x <- readRDS(file="./Combo Chemo Data/colorectal/GSE39582/GSE39582.adj.rds")
    dim(x) ## 233
    x.clinic <- readRDS(file="./Combo Chemo Data/colorectal/GSE39582/GSE39582.adj.clinic.rds")
    table(x.clinic$chemotherapy.type)
    
    id <- 1:ncol(x) ##grep( "treatment arm: Paclitaxel", as.character(christine279.clinic$treatment.arm) )
    
    if( data.name=="CRC_GSE39582_5FU+FOLFIRI"){
      id <- which(x.clinic$chemotherapy.type %in% c("5FU", "FOLFIRI", "FOLFOX", "FUFOL"))
    }else if( data.name=="CRC_GSE39582_5FU+FOLFIRI_Metastasis"){
      id <- which(x.clinic$chemotherapy.type %in% c("5FU", "FOLFIRI", "FOLFOX", "FUFOL") & x.clinic$tnm.m=="M1")
    }else if( data.name=="CRC_GSE39582_5FU+FOLFIRI_NoMetastasis"){
      id <- which(x.clinic$chemotherapy.type %in% c("5FU", "FOLFIRI", "FOLFOX", "FUFOL") & x.clinic$tnm.m=="M0")
      
    }else if( data.name=="CRC_GSE39582_5FU"){
      id <- which(x.clinic$chemotherapy.type %in% c("5FU"))
    }else if( data.name=="CRC_GSE39582_5FU_NoMetastasis"){
      id <- which(x.clinic$chemotherapy.type %in% c("5FU") & x.clinic$tnm.m=="M0")
    }else if( data.name=="CRC_GSE39582_5FU_Metastasis"){
      id <- which(x.clinic$chemotherapy.type %in% c("5FU") & x.clinic$tnm.m=="M1")
      
    }else if ( data.name=="CRC_GSE39582_FOLFIRI"){
      id <- which(x.clinic$chemotherapy.type %in% c("FOLFIRI", "FOLFOX", "FUFOL"))
    }else if(  data.name=="CRC_GSE39582_FOLFIRI_NoMetastasis"){
      id <- which(x.clinic$chemotherapy.type %in% c("FOLFIRI", "FOLFOX", "FUFOL") & x.clinic$tnm.m=="M0")
    }else if(  data.name=="CRC_GSE39582_FOLFIRI_Metastasis"){
      id <- which(x.clinic$chemotherapy.type %in% c("FOLFIRI", "FOLFOX", "FUFOL") & x.clinic$tnm.m=="M1")
      
    }else if ( data.name=="CRC_GSE39582_NoFOLFIRI"){
      id <- which(x.clinic$chemotherapy.type %in% c("5FU", "FOLFOX", "FUFOL"))
    }else if(  data.name=="CRC_GSE39582_NoFOLFIRI_NoMetastasis"){
      id <- which(x.clinic$chemotherapy.type %in% c("5FU", "FOLFOX", "FUFOL") & x.clinic$tnm.m=="M0")
    }else if(  data.name=="CRC_GSE39582_NoFOLFIRI_Metastasis"){
      id <- which(x.clinic$chemotherapy.type %in% c("5FU", "FOLFOX", "FUFOL") & x.clinic$tnm.m=="M1")
    }else{
      id <- which(x.clinic$chemotherapy.type %in% c("5FU", "FOLFIRI", "FOLFOX", "FUFOL") & x.clinic$tnm.m=="M0")
    }      
    
    evaluationset <- x[,id]
    evaluationset.clinic <- x.clinic[id, ]
    evaluationset.response <- NULL
    evaluationset.name=paste("CRC:GSE39582-", ncol(evaluationset))
  }
  
  
  
  
  ### GSE19860: Japan: Watanabe ## No raw data
  ### Response: response to either FL or Bevactzumab
  ## FOLFOX differs from FOLFIRI
  if( data.name=="CRC_GSE19860_Metastasis_FOLFOX+Bevacizumab"){
    ### colorectal
    x <- readRDS(file="./Combo Chemo Data/colorectal/GSE19860/GSE19860.rds")
    x.clinic <-  readRDS( file="./Combo Chemo Data/colorectal/GSE19860/GSE19860.clinic.rds")
    x.clinic$response.FOLFOX <- x.clinic$response.FOLFIRI
    if( tolower(compound.name)=="5fu"){
      x.clinic$response <- x.clinic$response.FOLFOX
    }else if( tolower(compound.name)=="bevacizumab") {
      x.clinic$response <- x.clinic$response.FOLFOX
    }else{
      x.clinic$response <- x.clinic$response.FOLFOX
    }

        
    x.clinic$response <- as.numeric( x.clinic$response)
    
    id <- which( !is.na(x.clinic$response)) ###1:ncol(x) ##grep( "treatment arm: Paclitaxel", as.character(christine279.clinic$treatment.arm) )
    
    evaluationset <- x[,id]
    evaluationset.clinic <- x.clinic[id, ]
    evaluationset.response <- evaluationset.clinic$response
    evaluationset.name=paste("CRC:GSE19860-", ncol(evaluationset))
  }
  
  
  ### GSE52735
  ### GSE52847 as well as GSE52735(validation)
  if( data.name=="CRC_GSE52735_Fluoropyrimidine"){
    ### colorectal
    x <- readRDS(file="../Combo Chemo Data/colorectal/GSE52735/GSE52735.rds")
    x.clinic <-  readRDS( file="../Combo Chemo Data/colorectal/GSE52735/GSE52735.clinic.rds")

    x.clinic$response <- as.numeric( x.clinic$response.binary)
    
    id <- which( !is.na(x.clinic$response)) ###1:ncol(x) ##grep( "treatment arm: Paclitaxel", as.character(christine279.clinic$treatment.arm) )
    
    evaluationset <- x[,id]
    evaluationset.clinic <- x.clinic[id, ]
    evaluationset.response <- evaluationset.clinic$response
    evaluationset.name=paste("CRC:GSE52735-", ncol(evaluationset))
  }
  
    
  #### Its outcome is not a drug response but having metastasis or not. 
  if( data.name=="CRC_GSE14095_LiverMetastasis"){
    ### colorectal
    x <- readRDS(file="../Combo Chemo Data/colorectal/GSE14095/GSE14095.rds")
    x.clinic <-  readRDS( file="../Combo Chemo Data/colorectal/GSE14095/GSE14095.clinic.rds")
    x.clinic$response <- x.clinic$response.treatment
    id <- which( !is.na(x.clinic$response) & x.clinic$DM=="P") ###1:ncol(x) ##grep( "treatment arm: Paclitaxel", as.character(christine279.clinic$treatment.arm) )
    
    evaluationset <- x[,id]
    evaluationset.clinic <- x.clinic[id, ]
    evaluationset.response <- evaluationset.clinic$response
    evaluationset.name=paste("CRC:GSE14095-", ncol(evaluationset))
  }
  
  
  #Pancreatic:EMEXP2780: Google Goest Cancer
  if( data.name=="Pancreatic:EMEXP2780"){
    ### Pancreatic
    EMEXP2780 <- readRDS("./Combo Chemo Data/pancreatic/E-MEXP-2780/EMEXP2780.rds")
    EMEXP2780.clinic <- readRDS("./Combo Chemo Data/pancreatic/E-MEXP-2780/EMEXP2780.clinic.rds")

    EMEXP2780.clinic$valid <- 1
    EMEXP2780.clinic$OS <- EMEXP2780.clinic$OS/30

    #if( compound.name %in% c("Gemcitabine") ) id <- 1:ncol(EMEXP2780)
    evaluationset <- EMEXP2780
    evaluationset.response <- NULL
    evaluationset.clinic <- EMEXP2780.clinic
  } 
  
  if( data.name=="Pancreatic:GSE57495:Moffitt"){
    evaluationset <- readRDS(file="./Combo Chemo Data/pancreatic/GSE57495/GSE57495.u133p2.Average.rds")
    
    evaluationset.clinic <- readRDS(file="./Combo Chemo Data/pancreatic/GSE57495/GSE57495.u133p2.clinic.rds")
    evaluationset.response <- NULL
  }
  
  if( data.name=="Pancreatic:ICGC_RNAseq"){
    evaluationset <- readRDS( file="./Combo Chemo Data/pancreatic/nature2016/expr.u133p2.ICGC.rds") ## MeanMax probeset selection
    evaluationset.clinic <- readRDS(file="./Combo Chemo Data/pancreatic/nature2016/clinic.u133p2.ICGC.rds") ## Average probeset integration
    evaluationset.response <- NULL
  }    
  
  
  # Pancreatic:GSE17891CellLines (n=19)
  if( data.name=="Pancreatic:GSE17891CellLines"){
    
    ### Pancreatic: GSE17891: cancer cell lines 
    GSE17891.cell <- readRDS("./Combo Chemo Data/pancreatic/GSE17891/GSE17891.cell.rds")
    GSE17891.cell.clinic <- readRDS("./Combo Chemo Data/pancreatic/GSE17891/GSE17891.cell.clinic.rds")
    GSE17891.cell.clinic$valid <- 1
    
    #if( compound.name %in% c("Gemcitabine") ) id <- 1:ncol(EMEXP2780)
    evaluationset <- GSE17891.cell
    evaluationset.clinic <- GSE17891.cell.clinic
    evaluationset.response <- NULL
    if( length(grep("gemcitabine", tolower(compound.name)))>=1  ){
      evaluationset.response <- evaluationset.clinic$IC50.gemcitabine
      cat("We read IC50 gemcitabine\n")
    }else if( length(grep("erlotinib", tolower(compound.name)))>=1  ){
      evaluationset.response <- evaluationset.clinic$IC50.erlotinib
      cat("We read IC50 erlotinib\n")
    }
  } 
  
  
  if( data.name=="Pancreatic:GSE17891:Patient"){
    #if( compound.name %in% c("Gemcitabine") ) id <- 1:ncol(EMEXP2780)
    evaluationset <- readRDS(file="../Combo Chemo Data/pancreatic/GSE17891/GSE17891.patient.rds")
    evaluationset.clinic <- readRDS(file="../Combo Chemo Data/pancreatic/GSE17891/GSE17891.clinic.supple.rds")
    evaluationset.response <- NULL
  } 
  
  
  
  ## Pancreatic: GSE28735: HuEx-st-(gene version): OS outcome
  if( data.name=="Pancreatic_GSE28735_HuEXst1"){
    evaluationset <- readRDS(file="../Combo Chemo Data/pancreatic/GSE28735/GSE28735.tumor.u133p2.rds")
    evaluationset.clinic <- readRDS(file="../Combo Chemo Data/pancreatic/GSE28735/GSE28735.tumor.u133p2.clinic.rds")
    evaluationset.response <- NULL      
  }
  
  
  ## Pancreatic: GSE28735: HuEx-st-(gene version): OS outcome
  if( data.name=="Pancreatic_GSE28735_HuEXst1"){
    evaluationset <- readRDS(file="../Combo Chemo Data/pancreatic/GSE28735/GSE28735.tumor.u133p2.rds")
    evaluationset.clinic <- readRDS(file="../Combo Chemo Data/pancreatic/GSE28735/GSE28735.tumor.u133p2.clinic.rds")
    evaluationset.response <- NULL      
  }
  
 
  ### GSE21501: Pacnreatic cancer patients from 5 different centers
  if( data.name=="Pancreatic:GSE21501.Agilent"){
    evaluationset <- readRDS(file="../Combo Chemo Data/pancreatic/GSE21501/GSE21501.u133p2.rds")
    evaluationset.clinic <- readRDS(file="../Combo Chemo Data/pancreatic/GSE21501/GSE21501.u133p2.clinic.rds")
    cat("Finish initializing Test Data set\n")
    evaluationset.response <- NULL  
  }  
  
  
  
  ## NSCLC:GSE31625CellLines (n=49)
  # it has response to erlotinib
  if( data.name=="NSCLC:GSE31625CellLines:Erlotinib"){
    
    ### NSCLC: GSE31625 data set: Balko dataset
    ##load("./Combo Chemo Data/lung/GSE31625-cell lines and erlotinib sensitivity/GSE31625.total.RData")
    load("./Combo Chemo Data/lung/GSE31625-cell lines and erlotinib sensitivity/GSE31625.total.RData")
    # binary erlotinib sensitivity
    GSE31625.clinic$response <- car::recode(GSE31625.clinic$erlotinib.sensitivity, recode="c('Yes')=1;c('No')=0;else=NA")
    GSE31625.clinic$drug.response <- GSE31625.clinic$response
    
      colnames(GSE31625.clinic)
      ## 46 NSCLC cell line, 1 epidermoid carcinoma, 1 CML cell line, 
      table(GSE31625.clinic$Sample_source_name)
      table(GSE31625.clinic$tissue)
      ## match(marker, rownames(GSE31625))
    
    ## 46 NSCLC cell line, 1 epidermoid carcinoma, 1 CML cell line, 
      table(GSE31625.clinic$Sample_source_name)
      table(GSE31625.clinic$tissue)
    ## match(marker, rownames(GSE31625))
    ##id <- 1:ncol(GSE31625) 
    ## After QC, the 48th sample showed strange boxplot and a singleton cluster in PCA cluster analysis.
    ##id <- 1:47
    id <- which( GSE31625.clinic$CELL.TYPE=="non small cell lung cancer cell (NSCLC)")
    # id; length(id)
    # id <- which( GSE31625.clinic$CELL.TYPE!="chronic myelogenous leukemia (CML)")
    evaluationset <- GSE31625[,id]
    evaluationset.clinic <- GSE31625.clinic[id, ]
    evaluationset.response <- GSE31625.clinic$response[id] 
    attr(evaluationset, "name") <- "NSCLC:GSE31625:erlotinib"
  } 
  
  
  
  ## NSCLC:GSE4824
  # it has both NSCLC and SCLC. Response to Gefitinib is available only in 39 NSCLC cell lines
  if( data.name=="NSCLC:GSE4824:Gefitinib"){
    GSE4824 <- readRDS(file="./Combo Chemo Data/lung/GSE4824-Gefitinib/GSE4824.rds")
    GSE4824.clinic <- readRDS(file="./Combo Chemo Data/lung/GSE4824-Gefitinib/GSE4824.clinic.rds")
    # cbind(colnames(GSE4824), rownames(GSE4824.clinic))
    # table(GSE4824.clinic$Histology, GSE4824.clinic$response)
    
    ## select cases with data on gefitinib response
    id <- which( !is.na(GSE4824.clinic$response));
    evaluationset <- GSE4824[,id]
    evaluationset.clinic <- GSE4824.clinic[id, ]
    evaluationset.response <- GSE4824.clinic$response[id] 
    attr(evaluationset, "name") <- paste("NSCLC:GSE4824:Gefitinib", nrow(evaluationset.clinic)); attr(evaluationset, "name") 
  } 
 
  
  
   
  ##  NSCLC:GSE37138:Bevacizumab+Erlotinib (n=117)
  ## original platform: HuEX-st1.0
  ## Using Best-Match file, u133p2..
  if( data.name=="NSCLC:GSE37138:Bevacizumab+Erlotinib"){
    x <- readRDS(file="./Combo Chemo Data/lung/GSE37138-bevacizumab+erlotinib/GSE37138.u133p2.rds")
    x.clinic <- readRDS(file="./Combo Chemo Data/lung/GSE37138-bevacizumab+erlotinib/GSE37138.clinic.rds")
    # cbind(colnames(x), rownames(x.clinic))
    table(x.clinic$tissue) ## 75 Blood samples and 42 Bronchoscoptic NSCLCs.
    ## select cases with data on drug response
    id <- which( x.clinic$tissue=="Bronchoscopic lung biopsies from NSCLC patients" )
    
    #evaluationset[1:3, 1:3]    
    evaluationset <- x[,id] # data.frame(x)[probe.u133a, id]
    # rownames(evaluationset) <- probe.u133a
    # evaluationset[ is.na(evaluationset) ] <- 0
    evaluationset.clinic <- x.clinic[id, ]
    evaluationset.response <- x.clinic$response[id] 
    attr(evaluationset, "name") <- paste("NSCLC:GSE37138:bevacizumab+erlotinib", nrow(evaluationset.clinic)); attr(evaluationset, "name") 
  } 
  
  
  ##  NSCLC:GSE33072:Sorafenib | Erlotinib (n=117)
  ## original platform: HuEX-st1.0
  ## original platform: HuGene-st1.0 (transcript)
  ## Using Best-Match file, u133p2..
  if( data.name=="NSCLC:BATTLE(Sorafenib|Erlotinib)"){
    x <- readRDS(file="./Combo Chemo Data/lung/GSE33072-BATTLE-sorafenib+erlotinib/GSE33072.u133p2.rds")  
    x.clinic <- readRDS(file="./Combo Chemo Data/lung/GSE33072-BATTLE-sorafenib+erlotinib/GSE33072.clinic.20181015.rds") ## need to further check.

    ## select cases with data on drug response
    id <- 1:nrow(x.clinic)
    if( tolower(compound.name)=="erlotinib"){
      id <- grep( "erlotinib", x.clinic$treatment) 
    }else if( tolower(compound.name)=="sorafenib"){
      id <- grep( "sorafenib", x.clinic$treatment) 
    }
    evaluationset <- x[, id]
    evaluationset.clinic <- x.clinic[id, ]
    evaluationset.response <- as.numeric(evaluationset.clinic$response)
  } 
  
  
  if( data.name=="NSCLC:BATTLE(All)"){
    x <- readRDS(file="./Combo Chemo Data/lung/GSE33072-BATTLE-sorafenib+erlotinib/GSE33072.u133p2.rds")  
    x.clinic <- readRDS(file="./Combo Chemo Data/lung/GSE33072-BATTLE-sorafenib+erlotinib/GSE33072.clinic.20181015.rds")
    
   id <- 1:nrow(x.clinic)

    evaluationset <- x[, id]
    evaluationset.clinic <- x.clinic[id, ]
    evaluationset.response <- as.numeric(evaluationset.clinic$response)
  } 
  
  if( data.name=="NSCLC:BATTLE(Sorafenib)"){
    x <- readRDS(file="./Combo Chemo Data/lung/GSE33072-BATTLE-sorafenib+erlotinib/GSE33072.u133p2.rds")  
    x.clinic <- readRDS(file="./Combo Chemo Data/lung/GSE33072-BATTLE-sorafenib+erlotinib/GSE33072.clinic.20181015.rds")
    id <- which(x.clinic$treatment=="sorafenib")
    evaluationset <- x[, id]
    evaluationset.clinic <- x.clinic[id, ]
    evaluationset.response <- as.numeric(evaluationset.clinic$PFS)
  } 
  
  

  if( data.name=="NSCLC:GSE68793:TCGA:NoErlotinib"){
    x <- readRDS( file="./Combo Chemo Data/lung/GSE68793_TCGA/TCGA.LUSC.20.rds")
    x.clinic <- readRDS(file="./Combo Chemo Data/lung/GSE68793_TCGA/TCGA.LUSC.20.clinic.rds")
    
    id <- which( x.clinic$erlotinib==0)
    evaluationset <- x[, id]
    evaluationset.clinic <- x.clinic[id, ]
    evaluationset.response <- NULL
    
  }
  
    
  if( !is.null(evaluationset.response) ) names(evaluationset.response) <- colnames(evaluationset)
  rownames(evaluationset.clinic) <- colnames(evaluationset)
  evaluationset.clinic$is.valid <- 1
  
  attr(evaluationset, "name") <- paste(data.name, nrow(evaluationset.clinic), sep=""); attr(evaluationset, "name") 
  
  evaluationset.name=paste(gsub("[:|]", "_", data.name), ncol(evaluationset), sep="_")
  evaluationset.name
  
  x <- list(expr=evaluationset, drug.response=evaluationset.response, 
            clinic=evaluationset.clinic, name=evaluationset.name)
  x
}




## x: expression data
## marker: genes to be selected
scaleExpressionData.v02 <- function(x, marker){
  # x <- x$expr
  # marker=marker
  n.miss <- sum( is.na(match(marker, rownames(x))));
  cat(n.miss, "probsets are missing (code line=1271) \n")
  if( n.miss==0){ ## if trainingset and test are made on the same array platform
    x <- x[marker, ] 
    ##### scale a test set for each gene 
    x.scaled <- t( apply(x, 1, scale))
    
  }else if( n.miss > 0){ ## else..(different platforms were used)
    x <- data.frame(x)[marker, ]
    rownames(x) <- marker
    ## scale a test set for each gene 
    x.scaled <- t( apply(x, 1, scale))
    data.class( x.scaled)
    x.scaled[is.na(x.scaled)] <- 0
  }
  
  ## fill '0' for probsets with missing values
  id.na <- which(is.na(x.scaled), arr.ind=T)
  if( nrow(id.na)>0) x.scaled[id.na] <- 0
  cat(nrow(id.na), "probsets are missing (code line=1016) \n")
  colnames(x.scaled) <- colnames(x)
  
  return(x.scaled) ##x$expr[marker, ] 
}


scaleExpressionData <- function(x, marker){
  # x <- x$expr
  # marker=marker
  n.miss <- sum( is.na(match(marker, rownames(x))));
  cat(n.miss, "probsets are missing (code line=1271) \n")
  if( n.miss==0){ ## if trainingset and test are made on the same array platform
    x <- x[marker, ] 
    ##### scale a test set for each gene 
    x.scaled <- t( apply(x, 1, scale))
    
  }else if( n.miss > 0){ ## else..(different platforms were used)
    missing.marker <- marker[ is.na(match(marker, rownames(x))) ]
    tmp <- data.frame(x)[missing.marker, ]
    rownames(tmp) <- missing.marker
    
    present.marker <- marker[ !is.na(match(marker, rownames(x))) ]
    tmp2 <- data.frame(x[present.marker, ])
    x <- rbind(tmp, tmp2)[marker, ]
    
    #x <- data.frame(x)[marker, ]
    # rownames(x) <- marker
    ## scale a test set for each gene 
    x.scaled <- t( apply(x, 1, scale))
    colnames(x.scaled) <- colnames(x)
     data.class( x.scaled)
    x.scaled[is.na(x.scaled)] <- 0
  }
  
  ## fill '0' for probsets with missing values
  id.na <- which(is.na(x.scaled), arr.ind=T)
  if( nrow(id.na)>0) x.scaled[id.na] <- 0
  cat(nrow(id.na), "probsets are missing (code line=1016) \n")
  colnames(x.scaled) <- colnames(x)
  
  return(x.scaled) ##x$expr[marker, ] 
}


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
  
  
  
  
  
  ### 2018-01-19: 'biomarker.discovery.ttest: change 'qvalue' to 'p.adjust' because an error occurs: min(p)>1...
  ### 2018-01-20: 'biomarker.discoveyr.ttest: markerset.drug ==> markerset: to be in parallel with the same output 
  ###                 with  biomarker.discovery.correlation and other biomarker discovery functions. 
  ### 2018-01-30: Add 'how' option
  ## how=EX: exhasutive search; 
  ## how=EM; EM clustring
  
  biomarker.discovery.ttest.v02 <- function( markerset, markerset.response, qvalue.cutoff=0.2, 
                                             n.res=NULL, n.nonres=NULL, min.array=NULL, max.array=NULL, how="EX" ){
    # x <- readRDS( file="./RData/rds/5-Fluorouracil_EuroPDX-Breast_all_correlation-rank.rds")
    # markerset <- x$markerset; markerset.response <- x$markerset.response; n.res <- n.nonres <- min.array <- max.array <- NULL
    # qvalue.cutoff <- 0.2
    require( genefilter)
    
    n.array <- ncol(markerset)
    ##min.array <- 7 ## the minmum number of arrays in each samples for t-test
    if( is.null(max.array) ) max.array <- as.integer(0.75*n.array) #19
    if( is.null(min.array)  )  min.array <- max(as.integer(0.25*n.array), 2) # 25%
    cat("\n max.array=", max.array, "; min.array=", min.array, "\n")
    
    gene.count = matrix(0, nrow=max.array, ncol=max.array)
    ## why did it change the col,rownames?
    colnames(gene.count) = c(1:max.array)    # The start number of n.sesitive
    rownames(gene.count) = c(1:max.array)    # The start number of n.resistant
    
    drug.order <- order(markerset.response) ## Increasing order
    # markerset.response[order(markerset.response)] ## sort  from sensitive <-----> resistant
    par( mfrow=c(1,2))
    hist( markerset.response, nclass=15 )
    
    
    if( is.null(n.res) | is.null(n.nonres) ){
      for (i in min.array:max.array) {
        for (j in min.array:min((n.array-i),max.array) ) {
          cat(i); cat("\t"); cat(j);cat("\n")
          # i <-3; j <- 3
          n.res = i; 
          n.nonres = j
          
          drug.res = drug.order[1:n.res] ## smaller drug sensitivity values
          drug.nonres = rev(drug.order)[1:n.nonres] ## bigger drug sensitivity values
          
          # fast, but approximation
          x=markerset[,c(drug.res, drug.nonres)]
          y=factor(c( rep(1, n.res), rep(0, n.nonres)))
          
          s <- genefilter::rowttests(x=x, fac=y, tstatOnly=FALSE)
          ##temp <- rowttests( subsample, subsample$mol.biol )
          s$q.value <- p.adjust( s$p.value, method="fdr")      
          
          ## record the number of significant probesets for each # of response(sensitive) & # of non-response
          gene.sam = (1:nrow(markerset))[s$q.value < qvalue.cutoff]
          gene.count[j,i] = length(gene.sam)
        }
        cat(c(i, "."))
      }
      
      max.comb = which(gene.count == max(gene.count), arr.ind=T)
      #persp(1:max.array, 1:max.array, gene.count, xlab="n.sensitive", ylab="n.resistant", zlab="count of significant markers")
      
      # check max.comb
      n.res = max.comb[1,2]   # add to minimum number 7 for response
      n.nonres = max.comb[1,1] # add to minimum number 8 for non-response
    }
    
    
    #print( cbind(n.res,n.nonres) ) ## check result
    drug.res = drug.order[1:n.res]
    drug.nonres = rev(drug.order)[1:n.nonres]
    
    # generate binary response variables
    # That is, selecting a subset of NCI cellines which showed the sinificant difference.
    drug.o=c(drug.res, drug.nonres)
    
    markerset.drug <- markerset[, drug.o]
    markerset.drug.response <- markerset.response[drug.o] ## (11 sensitive drug response,  28 resistant drug response )
    markerset.drug.response.2 <- c(rep(1,n.res), rep(0,n.nonres))
    
    names(markerset.drug.response) <- names(markerset.drug.response.2) <- colnames(markerset.drug)
    
    ####################################################################################################
    # Biomarker method 2.1 - t.test + no coxen step
    ####################################################################################################
    #qval.cutoff = 0.2
    # fast, but approximation
    x=markerset[,c(drug.res, drug.nonres)]
    y=factor( markerset.drug.response.2)
    
    
    s <- genefilter::rowttests(x=x, fac=y, tstatOnly=FALSE)
    ##temp <- rowttests( subsample, subsample$mol.biol )
    s$q.value <- p.adjust( s$p.value, method="fdr")      
    
    ## select top 345 genes
    tstat <- s$statistic; names(tstat) <- rownames(s)
    tstat.ordered <- sort( abs(tstat), decreasing=TRUE)
    
    result <- list()
    xx <- data.frame(s) ##data.frame( mean.1=x.mean, mean.2=y.mean, teststat=t.stat, pvalue=t.pval, qvalue=t.qval )
    colnames(xx) <- c("teststat", "mean.diff", "pvalue", "qvalue")
    xx$estimate <- xx$mean.diff
    
    result$marker.info <- xx
    result$ordered.probe <- tstat.ordered
    
    result$markerset.origin <- markerset        
    
    
    ## a subset of input markerset,  used for t.test..
    result$markerset <- markerset
    result$markerset.response <- result$response <- markerset.response
    
    result$markerset.drug <- x
    result$markerset.drug.response <- markerset.drug.response ## continous
    result$markerset.drug.response.2 <- markerset.drug.response.2 ## binary
    result$method <- "ttest"
    
    return (result)
  }
  
  
  
  
  
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
    set1 <- ancillary.genes[ target, ] # COXEN-OV
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
    marker <- unlist(get.CCC()) ## Biomarkers after three-way COXEN analysis
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
    marker <- unlist( get.CCC() ) ## Biomarkers after three-way COXEN analysis
    
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
    
    marker <- unlist( get.CCC() ) ## Biomarkers after three-way COXEN analysis
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
    #xx <- dd.bCOXEN()$ff1
    #layout( matrix( c(1,1,2,2,2), nrow=1))
    #boxplot( x ~ group3, data=xx, las=2, ylim=c(-0.4, 1), col="green")
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
    ## coxen.score = -1*coxen.regression.score
    par(mfrow=c(2,2))
    #plot.strip(a=x$coxen.score, b=x$evaluationset.response)
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
    
    # roc(score=x$coxen.score, response=x$evaluationset.response, is.plot=FALSE)$power.table
    # plot(tt)
    ##wilcoxplot(x$coxen.regression.score ~ x$evaluationset.response)
    #paste0("x=",cli input$bestModel_click$x, "\ny=", input$bestModel_click$y)
    #xx <- evaluate.Predictors() ## output of 'select.model()'
    #rs <- data.frame(num=1:length(xx$p.value), pvalue=xx$p.value, aic=xx$aic)
    #nearPoints(df=rs, coordinfo=coord.mouse, xvar = "num", yvar = "pvalue")
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
  
  
  
  ### Step 5. Fit a final model & make a validation
  output$step5 <- renderPrint({
    GEM <- get.onePredictor()
    ##x <- get.validData()
    x=get.dataSet(input$validData)
    newdata <- x$expr
    
    head(predict.coxen(GEM=GEM, newdata=newdata))
  }) # end of 'output$step4 '
  
  
  output$step5validData <- renderPlot({
    GEM <- get.onePredictor() ## GEM 
    x <- get.validData() ## validation set
    newdata <- x$expr
    b <- x$drug.response;
    data.name <- input$validData
    
    score <- predict.coxen(GEM=GEM, newdata=newdata)
    ##coxen.score <- get.rankscore(-1*score)
    ## Multiply by -1, in order to make that the higer score, the more sensitive..
    plot.strip(a=-1*score, b=b, data.name=data.name)
    #title("Ovarian:TCGA-448")
  }) # end of 'output$step4 '
  
  
  output$step5LassoValidset <- renderPlot({
    #xx <- dd.bCOXEN()$ff1
    #layout( matrix( c(1,1,2,2,2), nrow=1))
    #boxplot( x ~ group3, data=xx, las=2, ylim=c(-0.4, 1), col="green")
    
    ## RF using all variables
    xx <- evaluate.Lasso() ## model evaluation results..
    ## get patient test data
    x <- get.validData()
    marker <- get.3way()
    x$expr <- x$expr[marker, ] 
    
    evaluationset.clinic <- x$clinic
    ## scale a test set for each gene 
    evaluationset <- t( apply(x$expr, 1, scale))
    evaluationset.response <- x$drug.response
    
    par(mfrow=c(1,3))
    score <- predict.glmnet(object=xx$fit, newx=t(evaluationset) )
    if( !is.null(evaluationset.response)){
      plot.strip(a=score, b=evaluationset.response)
    }
    #xx <- evaluate.RF.CV()
  }) # end of 'output$step4 '
  
  
  
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


## http://hostname:3838/APP_NAME/

# 1. can we confidenlty disprove ~
# quiescent and non-quiesecnet ..




## RX-648 t+Coxen biomarkers

