### 2013/8/1 : modified to save 'correlation type'
### 2013/7/31 : use 'apply' for faster calculation
biomarker.discovery.correlation <- function( markerset, markerset.response, method="spearman", 	plot.check=TRUE ){
  require( qvalue)
  n.gene=nrow(markerset)
  pval=rep(0,nrow(markerset))
  cor.value = cor.stat = vector()
  
  xx <- apply( markerset, 1, cor.test, y=markerset.response, method=method, exact=FALSE )
  for (i in 1:nrow(markerset)) {
    #xx <- cor.test(markerset.response , unlist(markerset[i,]), method=method, exact=FALSE)
    pval[i]=xx[[i]]$p.value
    cor.value[i]=xx[[i]]$estimate
    cor.stat[i]=xx[[i]]$statistic
    if (i/1000 == i%/%1000) cat(i,".")
  }
  cat("\n")
  
  pval.cor <- pval
  qval.cor = p.adjust(pval.cor, method="fdr")
  
  names(qval.cor) <- rownames(markerset)
  names(cor.value) <- rownames(markerset) #correlation
  names(pval.cor) <- rownames(markerset)
  names(cor.stat) <- rownames(markerset)  # teststat
  
  par(mfrow=c(1,3));hist( pval.cor); hist(qval.cor);
  ordered.probe <- sort( abs(cor.value), decreasing=T  )
  plot( ordered.probe, pval.cor[names(ordered.probe)], xlab="ordered.probe", ylab="P-value testing Correlation" )
  
  result <- list()
  xx <- data.frame( estimate=cor.value, teststat=cor.stat, pvalue=pval.cor, qvalue=qval.cor )
  rownames(xx) <- rownames(markerset)
  
  result$marker.info <- xx
  result$ordered.probe <- ordered.probe ##names(ordered.probe)
  
  result$markerset <- markerset
  result$response <- markerset.response
  result$method <- "correlation"
  result$type <- method
  
  if( plot.check==TRUE) plot( result[[1]]$estimate, result[[1]]$pvalue, pch=16 )
  return (result)
}
