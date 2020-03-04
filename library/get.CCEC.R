
get.CCEC <- function( set1, set2, probeset=probe.cor, cutoff=0.2, cor.cutoff=NULL, method="pearson",method.1="pearson"){
  #set1 <- rma
  #set2 <- rx
  #probeset <- target[120:130]
  
  set1.cor.cor = ( cor(t(set1[probeset, ]), method=method.1))    ## correlation among genes within PDX
  set2.cor.cor = ( cor(t(set2[probeset, ]), method=method.1) )   ## correlation among genes within patient
  
  n.cor.num <- length(probeset)
  ccc.cor<- matrix(1, nrow=n.cor.num, ncol=2)
  
  for (i in 1:n.cor.num) {
    # i <- 110;  print(i)
    ### Calculates Lin's (1989, 2000) concordance correlation coefficient for agreement on a continuous measure.
    ## epi.ccc 
    xx <- my.ccc(set1.cor.cor[i,-i], set2.cor.cor[i,-i], alternative="greater")        
    ccc.cor[i,1] = xx$rho.c$est
    ccc.cor[i,2] = xx$p.value #.asym
  }
  
  result <- list()
  qvalue <- p.adjust(ccc.cor[,2], method="fdr")
  result$probe <- probeset[ which( qvalue < cutoff ) ]
  
  if( !is.null(cor.cutoff) ){
    result$probe <- probeset[ which( (qvalue <= cutoff) & (ccc.cor[,1] >= cor.cutoff) ) ]
  }
  
  rownames(ccc.cor) <- probeset
  
  ### Put CCC 
  result$ccc.cor <- ccc.cor
  result$ccc <- result$ccc.cor[ result$probe, ]
  ### transformation of vector into a matrix
  if(length(result$probe)==1) result$ccc <- matrix(result$ccc.cor[ result$probe, ], nrow=1)
  colnames(result$ccc.cor) <- colnames(result$ccc) <- c("ccc", "ccc.p")
  
  return(result)
}
