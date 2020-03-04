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
