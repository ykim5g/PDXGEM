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