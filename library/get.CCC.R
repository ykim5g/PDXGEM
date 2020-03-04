#### CCC function from epiR
get.CCC <- function (x, y, ci = "z-transform", conf.level = 0.95, null.CCC=0, alternative="greater", is.blalt=FALSE) {
  ## data
  dat <- data.frame(x, y)
  ## extract complete cases
  id <- complete.cases(dat)
  nmissing <- sum(!complete.cases(dat))
  dat <- dat[id, ]
  
  ## Quantile for CI
  N. <- 1 - ((1 - conf.level)/2)
  zv <- qnorm(N., mean = 0, sd = 1)
  lower <- "lower"
  upper <- "upper"
  
  ## parameter
  k <-   length(dat$y)
  yb <-  mean(dat$y)
  sy2 <- var(dat$y) * (k - 1)/k
  sd1 <- sd(dat$y)
  
  xb <- mean(dat$x)
  sx2 <- var(dat$x) * (k - 1)/k
  sd2 <- sd(dat$x)
  
  ## correlation
  r <- cor(dat$x, dat$y)
  sl <- r * sd1/sd2
  sxy <- r * sqrt(sx2 * sy2)
  
  ### CCC value (p)
  p <- 2 * sxy/(sx2 + sy2 + (yb - xb)^2) ### CCC = 2*s12 / (s11 + s22 + (m1-m2)**2 )  = rho * C_b
  
  ## M (delta) & A (mean)
  blalt=NULL
  if( is.blalt==TRUE){
    delta <- (dat$x - dat$y)
    rmean <- apply(dat, MARGIN = 1, FUN = mean)
    blalt <- data.frame(mean = rmean, delta)
  }
  
  v <- sd1/sd2
  u <- (yb - xb)/(sqrt(sd1*sd2))
  
  C.b <- p/r
  
  
  ### Standard Deviation of CCC
  sep = sqrt(((1 - ((r)^2)) * (p)^2 * (1 - ((p)^2))/(r)^2 + 
                (4 * (p)^3 * (1 - p) * (u)^2/r) - 2 * (p)^4 * (u)^4/(r)^2)/(k -2))
  
  
  #############################
  test.ccc <- (p-null.CCC)/sep
  
  p.value <- 2*(1-pnorm(abs(test.ccc)))
  if( alternative=="greater") p.value <- 1-pnorm( test.ccc )
  if( alternative=="less") p.value <- pnorm(test.ccc)
  
  ll = p - zv * sep
  ul = p + zv * sep
  
  ### Fisher's Z transformation
  t <- log((1 + p)/(1 - p))/2 ## mean
  set = sep/(1 - ((p)^2))     ## Standard error
  
  llt = t - zv * set
  ult = t + zv * set
  
  test.ccc.asym = (t/set)
  p.value.asym <- 2*(1-pnorm(abs(test.ccc.asym)))
  if( alternative=="greater") p.value.asym <- 1-pnorm( test.ccc.asym )
  if( alternative=="less") p.value.asym <- pnorm(test.ccc.asym)
  
  
  ### inverse tangent
  llt = (exp(2 * llt) - 1)/(exp(2 * llt) + 1) ## asymtotic variance
  ult = (exp(2 * ult) - 1)/(exp(2 * ult) + 1)
  
  if (ci == "asymptotic") {
    rho.c <- as.data.frame(cbind(p, ll, ul))
    names(rho.c) <- c("est", lower, upper)
    rval <- list(rho.c = rho.c, s.shift = v, l.shift = u, 
                 C.b = C.b, blalt = blalt, nmissing = nmissing, teststat=test.ccc, p.value=p.value, p.value.asym=p.value.asym)
  }
  else if (ci == "z-transform") {
    rho.c <- as.data.frame(cbind(p, llt, ult))
    names(rho.c) <- c("est", lower, upper)
    rval <- list(rho.c = rho.c, s.shift = v, l.shift = u, 
                 C.b = C.b, blalt = blalt, nmissing = nmissing, teststat=test.ccc, p.value=p.value, p.value.asym=p.value.asym)
  }
  return(rval)
}

