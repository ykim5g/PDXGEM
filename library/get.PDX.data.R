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

