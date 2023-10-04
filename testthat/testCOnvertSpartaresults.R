convertSpartaResults <- function(out){
  if(!is.character(class(out)) && class(out) == "occDet") stop("Class must be of occDet")

  #extract number of chains
  n.chains <- out$BUGSoutput$n.chains
  samplesList  <- vector('list', n.chains)

  samplesList <- lapply(1:n.chains, function(x){
    ret <- as.data.frame(out$BUGSoutput$sims.array[,1,])  
  })
  names(samplesList)  <- paste0('chain', 1:n.chains)
  
  #Return results
  retList <- list()
  retList$samples <- samplesList
  retList$n.keep <- out$BUGSoutput$n.keep
  return(retList)
}

rrt <- convertSpartaResults(out)





ret <- coda::as.mcmc(out$model$state())
ret$BUGSoutput$sims.list

retSamps <- coda::as.mcmc(ret$BUGSoutput$sims.list)

library(coda)
ttr <- as.mcmc(out)
qr <- as.data.frame(out$BUGSoutput$sims.array[,1,])
qr[,130:140]

out$BUGSoutput$sims.list
