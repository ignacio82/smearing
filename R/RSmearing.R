fimpacts <- function(d, DF, theta, nPeriods){
  D <- ncol(theta)
  B <- nrow(DF)
  DF$epredC <- exp(DF$pred + theta[,d] * (DF$treat == 1))
  DF$epredT <- exp(DF$pred - theta[,d] * (DF$treat == 0))
  DF$diff <-
    as.vector(
      DF$eresid %*% matrix(
        DF$epredT, nrow = B, ncol = B, byrow = T
      ) -
        DF$eresid %*% matrix(
          DF$epredC, nrow = B, ncol = B, byrow = T
        )
    ) / B
  sapply(1:nPeriods, function(t)
    weighted.mean(DF$diff[DF$time == t], DF$weight[DF$time == t]))

}

ParallelRsmearing <- function(DF, theta, cores=2) {
  nPeriods <- max(DF$time)
  D <- ncol(theta)
  impacts <- parallel::mclapply(1:D, fimpacts, DF, theta, nPeriods, mc.cores = cores)
  impacts <- do.call("rbind", impacts)
  return(impacts)
}

Rsmearing <- function(DF, theta) {
  nPeriods <- max(DF$time)
  D <- ncol(theta)
  impacts <- sapply(1:D, fimpacts, DF, theta, nPeriods)
  return(impacts)
}
