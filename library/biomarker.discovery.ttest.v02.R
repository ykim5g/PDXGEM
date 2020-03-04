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
  # Biomarker method 2.1 - t.test 
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
