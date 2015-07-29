smearing <- function(DF, theta){
  nPeriods <- max(DF$time)
  nPeriodT <- sapply(1:8, function(n) as.integer(sum(DF$time==n)))
  ncoltheta <- ncol(theta)
  wm <- .Fortran("smearing", i = nrow(DF),
                 nMCd = ncoltheta,
                 nPeriodT = nPeriodT,
                 nPeriods = as.integer(nPeriods),
                 DF = unlist(DF),
                 theta = unlist(theta),
                 wm = numeric(ncoltheta*nPeriods))$wm
  dim(wm) <- c(nPeriods,ncoltheta)
  return(wm)
}
