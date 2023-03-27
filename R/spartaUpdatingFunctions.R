#library(usethis)
#usethis::edit_r_environ()
#when the tab opens up in R studio, add this to the 1st line: R_MAX_VSIZE=100Gb (or whatever memory you wish to allocate).
baselineSpartaEstimation <- function(model, #nimbleModel
                                     MCMCconfiguration = NULL, #configuration for MCMC
                                     pfType = "bootstrap",#Either 'auxiliary' or 'bootstrap'. Defaults to auxiliary
                                     pfControl = list(saveAll = TRUE,
                                                      #lookahead = "mean",
                                                      smoothing = FALSE), #list of controls for particle filter
                                     nParFiltRun = NULL, #Number of PF runs
                                     latent #the latent variable
){

  timeStart1 <- Sys.time()
  target = MCMCconfiguration[["target"]]
  additionalPars = MCMCconfiguration[["additionalPars"]] #other dependent variables you seek to monitor
  n.iter = MCMCconfiguration[["n.iter"]]
  n.chains = MCMCconfiguration[["n.chains"]]
  n.burnin = MCMCconfiguration[["n.burnin"]]
  n.thin = MCMCconfiguration[["n.thin"]]
  if(is.null(nParFiltRun)) nParFiltRun = 10000

  #create new model for weights
  estimationModel <- model$newModel(replicate = TRUE)

  message("Building particle filter for model")
  if(!is.null(pfType)){
    if(!pfType %in% c("auxiliary", "bootstrap")) stop("Function currently works for auxiliary and bootstap Particle filters")
    #particleFilter is used to run the MCMC
    #particleFilterEst is used for returning the weights at the posterior
    #values of the top-level nodes
    if(pfType == "auxiliary"){
      particleFilter <- nimbleSMC::buildAuxiliaryFilter(model,
                                                             latent,
                                                             control = pfControl)

      particleFilterEst <- nimbleSMC::buildAuxiliaryFilter(estimationModel,
                                                                latent,
                                                                control = pfControl)
    }

    if(pfType == "bootstrap"){
      particleFilter <- nimbleSMC::buildBootstrapFilter(model,
                                                             latent,
                                                             control = pfControl)

      particleFilterEst <-  nimbleSMC::buildBootstrapFilter(estimationModel,
                                                                 latent,
                                                                 control = pfControl)
    }
  }else{
    particleFilter <- nimbleSMC::buildBootstrapFilter(model,
                                                           latent,
                                                           control = pfControl)

    particleFilterEst <- nimbleSMC::buildBootstrapFilter(estimationModel,
                                                              latent,
                                                              control = pfControl)
  }

  message("Compiling the particle filter")
  #compiling the model
  compiledParticleFilter <- compileNimble(model,  particleFilter)

  #Loglikelihood of last run and the Effective sample sizes
  message("Running the particle filter")
  logLik <-   compiledParticleFilter$particleFilter$run(m = nParFiltRun)
  ESS <-   compiledParticleFilter$particleFilter$returnESS()

  message("Setting up the MCMC Configuration")
  #model <- model$newModel(replicate = TRUE)
  modelMCMCconf <- nimble::configureMCMC(model, nodes = NULL)

  if(is.null(pfType)) pfType = "bootstrap" #set bootstrap as default pfType

  modelMCMCconf$addSampler(target = target,
                           type = 'RW_PF_block',
                           control = list(latents = latent,
                                          pfControl = list(saveAll = TRUE),
                                          pfNparticles = nParFiltRun,
                                          pfType = pfType))

  modelMCMCconf$addMonitors(additionalPars)

  message("Building and compiling the PF MCMC")
  ## build and compile pMCMC sampler
  modelMCMC <- nimble::buildMCMC(modelMCMCconf)
  compiledList <- nimble::compileNimble(model,
                                        modelMCMC,
                                        resetFunctions = TRUE)

  message("Running the PF MCMC")
  #run MCMC
  timeStart2 <- Sys.time()
  mcmc.out <- nimble::runMCMC(compiledList$modelMCMC,
                              niter = n.iter,
                              nchains = n.chains,
                              nburnin = n.burnin,
                              #inits = initsList,
                              thin = n.thin,
                              setSeed = TRUE,
                              samples=TRUE,
                              samplesAsCodaMCMC = TRUE,
                              summary = TRUE,
                              WAIC = FALSE)
  timeEnd <- Sys.time()

  timetaken1 <- timeEnd - timeStart1
  timetaken2 <- timeEnd - timeStart2

#   message(paste("Estimating weights at posteror values of ", target))
#   #posteriorEstimates <- mcmc.out$summary$all.chains[target, 'Mean']
# expandTarget <- model$expandNodeNames(target)
#   for(i in 1:length(expandTarget)){
#     #estimationModel[[expandTarget[i]]] <- mcmc.out$summary$all.chains[expandTarget[i], 'Mean']
#     estimationModel[[expandTarget[i]]] <- mcmc.out$summary[expandTarget[i], 'Mean']
#   }
#
#   message("Compiling the estimation particle filter")
#   #compiling the model
#   compiledParticleFilterEst <- compileNimble(estimationModel,  particleFilterEst)
#
#   #Loglikelihood of last run and the Effective sample sizes
#   message("Running the estimation particle filter")
#   logLik <-   compiledParticleFilterEst$particleFilterEst$run(m = nParFiltRun)
#   ESS <-   compiledParticleFilterEst$particleFilterEst$returnESS()
#
#
#   #save weights and samples
#   message("Extracting the weights and samples from particle fiter")
#   weights <- as.matrix( compiledParticleFilterEst$particleFilterEst$mvWSamples, "wts")
#   unweightedSamples <- as.matrix( compiledParticleFilterEst$particleFilterEst$mvWSamples, latent)
#   weightedSamples <- as.matrix( compiledParticleFilterEst$particleFilterEst$mvEWSamples, latent)
#

  #message("Returning the results")

  message("Returning the results")
  #list to return
  retList <- list()
  retList$samples <- mcmc.out$samples
  retList$summary <- mcmc.out$summary
  retList$timeRun <-  timetaken2
  return(retList)

  # #list to return
  # returnList <- list(weights = weights,
  #                    logLike = NULL,
  #                    unweightedSamples = unweightedSamples,
  #                    weightedSamples = weightedSamples,
  #                    particleFilter = particleFilter,
  #                    mcmcSamplesAndSummary = mcmc.out,
  #                    timeTakenAll = timetaken1,
  #                    timeTakenRun = timetaken2,
  #                    ess = ESS)

  #return(returnList)
}



#This function runs MCMC using particle filter MCMC and returns weights

spartaNimWeights <- function(model, #nimbleModel
                             MCMCconfiguration = NULL, #configuration for MCMC
                             pfType = "bootstrap",#Either 'auxiliary' or 'bootstrap'. Defaults to auxiliary
                             pfControl = list(saveAll = TRUE,
                                              #lookahead = "mean",
                                              smoothing = FALSE), #list of controls for particle filter
                             nParFiltRun = NULL, #Number of PF runs
                             latent, #the latent variable
                             mcmc = TRUE # logical if MCMC was used or not
){

  timeStart1 <- Sys.time()
  target = MCMCconfiguration[["target"]]
  additionalPars = MCMCconfiguration[["additionalPars"]] #other dependent variables you seek to monitor
  n.iter = MCMCconfiguration[["n.iter"]]
  n.chains = MCMCconfiguration[["n.chains"]]
  n.burnin = MCMCconfiguration[["n.burnin"]]
  n.thin = MCMCconfiguration[["n.thin"]]
  smoothing = pfControl[["smoothing"]]

  estimationModel <- model$newModel(replicate = TRUE)

  if(mcmc == TRUE){
    nMCompile <- compileNimble(model)
    #
    cMCMC <- configureMCMC(model, monitors = c(target, latent))
    #
    bMCMC <- buildMCMC(cMCMC)
    #
    coMCMC <- compileNimble(bMCMC, project = nMCompile)
    #
    timeStart2 <- Sys.time()
    mcmc.out <- runMCMC(coMCMC,
                            niter = n.iter,
                            nchains = n.chains,
                            nburnin = n.burnin,
                            setSeed = TRUE,
                            samples=TRUE,
                            samplesAsCodaMCMC = TRUE,
                            summary = TRUE,
                            WAIC = FALSE)

    timeEnd <- Sys.time()

    timetaken1 <- timeEnd - timeStart1
    timetaken2 <- timeEnd - timeStart2

#     #save weights and samples
#     message("Extracting the weights and samples from particle fiter")
#     latentNodes <- model$expandNodeNames(latent)
#
#     m <- n.iter - n.burnin
# #
#     #save weight, weighted and unweighted samples
#     weights <- matrix(1, nrow = m, ncol = length(latentNodes))
#     unweightedSamples <- mcmc.out$samples$chain1[, latentNodes]
#     weightedSamples <- mcmc.out$samples$chain1[, latentNodes]
#
#     message("Saving unsampled and sampled values in model values for updating")
#     namesEst <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
#     names <- c(target, latent)
#     type <- rep("double", length(names))
#     size <- as.list(sapply(names, function(x)length(estimationModel$expandNodeNames(x))))
#
#     mvSamplesEst <- modelValues(modelValuesConf(vars = names,
#                                                 types = type,
#                                                 sizes = size))
#     resize(mvSamplesEst, m)
# #mvSamplesEst <- list()
# logLike <- matrix(NA, nrow = m, ncol = length(latentNodes))
# for(i in 1: m){
# for(j in 1:length(names)){
#       mvSamplesEst[[names[j]]][[i]] <- mcmc.out$samples$chain1[i, estimationModel$expandNodeNames(names[j])]
#       nimCopy(from = mvSamplesEst, to = estimationModel, nodes = names[j],row = i)
#     }
#
#     for(k in 1:length(latentNodes)){
#       logLike[i,k] <- 0
#     }
# }

#Model values for the weighted samples
mvEWS <- mvWS <- mvSamplesEst <- NULL
logLike <- NULL
ESS <- particleFilter <- NULL
unweightedSamples = NULL
weightedSamples = NULL
weights = NULL

  }else{
  if(is.null(nParFiltRun)) nParFiltRun = 10000

  #create new model for weights
  estimationModel <- model$newModel(replicate = TRUE)

  message("Building particle filter for model")
  if(!is.null(pfType)){
    if(!pfType %in% c("auxiliary", "bootstrap")) stop("Function currently works for auxiliary and bootstap Particle filters")
#particleFilter is used to run the MCMC
#particleFilterEst is used for returning the weights at the posterior
    #values of the top-level nodes
    if(pfType == "auxiliary"){
      particleFilter <- myphdthesis::buildAuxiliaryFilterNew(model,
                                                             latent,
                                                             control = pfControl)

        particleFilterEst <- myphdthesis::buildAuxiliaryFilterNew(estimationModel,
                                                             latent,
                                                             control = pfControl)
    }

    if(pfType == "bootstrap"){
      particleFilter <- myphdthesis::buildBootstrapFilterNew(model,
                                                        latent,
                                                        control = pfControl)

        particleFilterEst <-  myphdthesis::buildBootstrapFilterNew(estimationModel,
                                                              latent,
                                                              control = pfControl)
    }
  }else{
    particleFilter <- myphdthesis::buildBootstrapFilterNew(model,
                                                           latent,
                                                           control = pfControl)

      particleFilterEst <- myphdthesis::buildBootstrapFilterNew(estimationModel,
                                                           latent,
                                                           control = pfControl)
  }

  message("Compiling the particle filter")
  #compiling the model
  compiledParticleFilter <- compileNimble(model,  particleFilter)

  #Loglikelihood of last run and the Effective sample sizes
  #message("Running the particle filter")
  #logLik <-   compiledParticleFilter$particleFilter$run(m = nParFiltRun)
  #ESS <-   compiledParticleFilter$particleFilter$returnESS()

  message("Setting up the MCMC Configuration")
  #model <- model$newModel(replicate = TRUE)
  modelMCMCconf <- nimble::configureMCMC(model, nodes = NULL)

  if(is.null(pfType)) pfType = "bootstrap" #set auxiliary as default pfType

  modelMCMCconf$addSampler(target = target,
                           type = 'RW_PF_block',
                           control = list(latents = latent,
                                          pfControl = list(saveAll = TRUE),
                                          pfNparticles = nParFiltRun,
                                          pfType = pfType))

  modelMCMCconf$addMonitors(additionalPars)

  message("Building and compiling the PF MCMC")
  ## build and compile pMCMC sampler
  modelMCMC <- nimble::buildMCMC(modelMCMCconf)
  compiledList <- nimble::compileNimble(model,
                                        modelMCMC,
                                        resetFunctions = TRUE)

  message("Running the PF MCMC")
  #run MCMC
  timeStart2 <- Sys.time()
  mcmc.out <- nimble::runMCMC(compiledList$modelMCMC,
                              niter = n.iter,
                              nchains = n.chains,
                              nburnin = n.burnin,
                              #inits = initsList,
                              thin = n.thin,
                              setSeed = TRUE,
                              samples=TRUE,
                              samplesAsCodaMCMC = TRUE,
                              summary = TRUE,
                              WAIC = FALSE)
timeEnd <- Sys.time()

timetaken1 <- timeEnd - timeStart1
timetaken2 <- timeEnd - timeStart2

 message(paste("Estimating effective sample size"))
 if(n.chains > 1){
 posteriorEstimates <- mcmc.out$summary$all.chains[target, 'Mean']

 expandTarget <- model$expandNodeNames(target)
 for(i in 1:length(expandTarget)){
   estimationModel[[expandTarget[i]]] <- mcmc.out$summary$all.chains[expandTarget[i], 'Mean']
 }
 }else{
   posteriorEstimates <- mcmc.out$summary[target, 'Mean']

   expandTarget <- model$expandNodeNames(target)
   for(i in 1:length(expandTarget)){
     estimationModel[[expandTarget[i]]] <- mcmc.out$summary[expandTarget[i], 'Mean']
   }
}
# message("Compiling the estimation particle filter")
 #compiling the model
 compiledParticleFilterEst <- compileNimble(estimationModel,  particleFilterEst)
#
# #Loglikelihood of last run and the Effective sample sizes
 message("Running the estimation particle filter")
 logLik <-   compiledParticleFilterEst$particleFilterEst$run(m = nParFiltRun)
 ESS <-   compiledParticleFilterEst$particleFilterEst$returnESS()
#
#
# #save weights and samples
# message("Extracting the weights and samples from particle fiter")
# weights <- as.matrix( compiledParticleFilterEst$particleFilterEst$mvWSamples, "wts")
# unweightedSamples <- as.matrix( compiledParticleFilterEst$particleFilterEst$mvWSamples, latent)
# weightedSamples <- as.matrix( compiledParticleFilterEst$particleFilterEst$mvEWSamples, latent)
# if(pfType == "auxiliary"){
#   logLike <- as.matrix( compiledParticleFilterEst$particleFilterEst$mvWSamples, "auxlog")
# }else{
#   logLike <- as.matrix( compiledParticleFilterEst$particleFilterEst$mvWSamples, "bootLL")
# }
#
#
# message("Saving unsampled and sampled values in model values for updating")
# m = getsize(compiledParticleFilterEst$particleFilterEst$mvWSamples)
# names <- compiledParticleFilterEst$particleFilterEst$mvWSamples$varNames
# size <- compiledParticleFilterEst$particleFilterEst$mvWSamples$sizes
# type <- c("double", "double", "double")
#
# if(pfType == "bootstrap"){
# namesEWS <- names[!names %in% c("bootLL", "wts")]
# }else{
#   namesEWS <- names[!names %in% c("bootLL", "auxlog")]
# }
# index <- which(names == namesEWS)
# mvEWS <- modelValues(modelValuesConf(vars = namesEWS,
#                                      types = type[index],
#                                      sizes = size[index]))
# if(smoothing){
#   size$wts <- 1
#   size$bootLL <- 1
# }
# mvWS  <- modelValues(modelValuesConf(vars = names,
#                                      types = type,
#                                      sizes = size))
#
# resize(mvWS, m)
# resize(mvEWS, m)
#
# for(i in 1: m){
#   mvEWS[[namesEWS]][[i]] <- weightedSamples[i,]
#   mvWS[[namesEWS]][[i]] <- unweightedSamples[i,]
# }
#
# mvSamplesEst <- mvEWS

 mvEWS <- mvWS <- mvSamplesEst <- NULL
}

  message("Returning the results")
  # #list to return
  # returnList <- list(#weights = weights,
  #                    #logLike = logLike,
  #                    #unweightedSamples = unweightedSamples,
  #                    #weightedSamples = weightedSamples,
  #                    particleFilter = particleFilter,
  #                    mcmcSamplesAndSummary = mcmc.out,
  #                    timeTakenAll = timetaken1,
  #                    timeTakenRun = timetaken2,
  #                    ess = ESS,
  #                    #compiledParticleFilterEst = compiledParticleFilterEst,
  #                    mvWS = mvWS,
  #                    mvEWS = mvEWS,
  #                    mvSamplesEst = mvSamplesEst)
  #
  # return(returnList)
  message("Returning the results")
  #list to return
  retList <- list()
  retList$samples <- mcmc.out$samples
  retList$summary <- mcmc.out$summary
  retList$timeRun <-  timetaken2
  return(retList)
}




##################################
#save weights and samples
#message("Extracting the weights and samples from particle fiter")

updateUtils <- function(model, #reduced model
                        mcmcOut,
                        latent, target, n.iter, m, timeIndex){

latentNodes <- model$expandNodeNames(latent)
nodes <- findLatentNodes(model, latent, timeIndex)

dims <- lapply(nodes, function(n) nimDim(model[[n]]))

if(length(unique(dims)) > 1)
  stop('sizes or dimensions of latent states varies')
vars <- c(model$getVarNames(nodes =  nodes), target)

modelSymbolObjects <- model$getSymbolTable()$getSymbolObjects()[vars]

names <- sapply(modelSymbolObjects, function(x)return(x$name))
type <- sapply(modelSymbolObjects, function(x)return(x$type))
size <- lapply(modelSymbolObjects, function(x){
  y <- x$size
  t <- length(y)
  rr <- c()
  if(t > 1){
rr <- y
  }else{
    if(length(y)>0){
   rr <- y
    }else{
 rr <- 1}
  }
  return(rr)
  }
)


# #names of targets
# modelSymbolObjects1 <- model$getSymbolTable()$getSymbolObjects()[target[2]]
# names1 <- sapply(modelSymbolObjects1, function(x)return(x$name))
# type1 <- sapply(modelSymbolObjects1, function(x)return(x$type))
# size1 <- lapply(modelSymbolObjects1, function(x)return(x$size))
#
#
#
# #expand target names
# targetAsScaler <- model$expandNodeNames(target, returnScalarComponents = TRUE)
# names <- c(names, targetAsScaler)
# type <- c(type, rep("double", length(targetAsScaler)))
#
# #create additional sizes for the target vars
# addSize <- list()
# for(i in 1:length(targetAsScaler)){
# addSize[[i]] <- length(dims)
# #names(size)[i + 1] <- target[i]
# }
# names(addSize) <- targetAsScaler
#
# #include them in the sizes
# size <- c(size, addSize)

mvSamplesEst <- modelValues(modelValuesConf(vars = names,
                                            types = type,
                                            sizes = size))
#my_initializeModel <- initializeModel(model, silent = silent)


# Create mv variables for x state and sampled x states.  If saveAll=TRUE,
# the sampled x states will be recorded at each time point.

#if(mcmc == TRUE){
# #save weight, weighted and unweighted samples
# if(timeIndex = 1){
# weights <- matrix(1, nrow = n.iter, ncol = length(latentNodes))
# unweightedSamples <- mcmcOut[, latentNodes]
# weightedSamples <- mcmcOut[, latentNodes]
# names <- c(target, latent)
# type <- rep("double", length(names))
# size <- as.list(sapply(names, function(x)length(model$expandNodeNames(x))))
#
# mvSamplesEst <- modelValues(modelValuesConf(vars = names,
#                                             types = type,
#                                             sizes = size))
#
# logLike <- matrix(NA, nrow = n.iter, ncol = length(latentNodes))
# }else{
#   weights <- matrix(1, nrow = n.iter, ncol = dims[[1]])
#   unweightedSamples <- NA
#   weightedSamples <- NA
#   names <- c(target, latent)
#   type <- rep("double", length(names))
#   size <- as.list(sapply(names, function(x)length(model$expandNodeNames(x))))
#
#   size[[length(names)]] <- dims[[1]]
#   mvSamplesEst <- modelValues(modelValuesConf(vars = names,
#                                               types = type,
#                                               sizes = size))
# }

message("Saving unsampled and sampled values in model values for updating")
#namesEst <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)

 resize(mvSamplesEst, n.iter)
#mvSamplesEst <- nimbleList(type = double())
#nimbleListTypes <- list(nimbleType(name = 'est', type = 'character'))

## this nimbleList definition is identical to the one created above
#mvSamplesEst <- nimbleList(nimbleListTypes)


  # mvSamplesEst[[iter]] <- modelValues(modelValuesConf(vars = names,
  #                                             types = type,
  #                                             sizes = size))


for(iter in 1:n.iter){
  #for(i in 1: m){
  for(j in 1:length(names)){
    if(names[j] == latent & length(size[[1]]) > 1){
    mvSamplesEst[[names[j]]][[iter]] <- matrix(mcmcOut[iter, model$expandNodeNames(names[j])], nrow = size[[1]][1], ncol = size[[1]][2])
    }else{
      mvSamplesEst[[names[j]]][[iter]] <-  mcmcOut[iter, model$expandNodeNames(names[j])]
    }
    #nimCopy(from = mvSamplesEst, to = model, nodes = names[j],row = i)
  #}
  }

  # for(k in 1:length(latentNodes)){
  #   logLike[,k] <- 0
  # }
}


#Model values for the weighted samples
mvEWS <- mvWS <- mvSamplesEst
#}else{
#  mvEWS <- mvWS <- mvSamplesEst < NA
#}

returnlist = list(weights = NA,
  logLike = NA,
  unweightedSamples = NA,
  weightedSamples = NA,
  #compiledParticleFilterEst = compiledParticleFilterEst,
  mvWS = mvWS,
  mvEWS = mvEWS,
  mvSamplesEst = mvSamplesEst
)
return(returnlist)
}

###########################################

spartaNimUpdates <- function(model, #nimbleModel
                             reducedModel,
                             MCMCconfiguration = NULL, #configuration for MCMC
                             pfType = NULL,#Either 'auxiliary' or 'bootstrap'. Defaults to auxiliary
                             pfControl = NULL, #list of controls for particle filter
                             nParFiltRun = NULL, #Number of PF runs
                             #updatePFControl = list(iNodePrev = NULL, M = NULL),
                             latent, #the latent variable
                             #newData,
                             postReducedMCMC,
                             target
){

  target = MCMCconfiguration[["target"]]
  additionalPars = MCMCconfiguration[["additionalPars"]] #other dependent variables you seek to monitor
  n.iter = MCMCconfiguration[["n.iter"]]
  n.chains = MCMCconfiguration[["n.chains"]]
  n.burnin = MCMCconfiguration[["n.burnin"]]
  n.thin = MCMCconfiguration[["n.thin"]]
  iNodePrev = pfControl[["iNodePrev"]]
  M = pfControl[["M"]]
  timeIndex <- pfControl[["timeIndex"]]
  #iNodePrev = updatePFControl[["iNodePrev"]]
  #M = updatePFControl[["M"]]
  if(is.null(nParFiltRun)) nParFiltRun = 10000
  if(is.null(timeIndex)) timeIndex = 1

  samplesList  <- vector('list', n.chains)

  samplesList <- lapply(as.list(1:n.chains), function(chain.iter){
    timeStart1 <- Sys.time()

    updateVars <- updateUtils(model = reducedModel, #reduced model
                #mcmcOut = postReducedMCMC$samples$chain1,
                mcmcOut = postReducedMCMC$samples[[chain.iter]],
                          latent = latent,
              target = target,
              n.iter = n.iter ,
              m = nParFiltRun,
              timeIndex = timeIndex)



  weights = updateVars$weights #weights from reduced model
  unweightedLatentSamples = updateVars$mvWS #saved model values for unweighted samples from reduced model
  weightedLatentSamples = updateVars$mvEWS  #saved model values for weighted samples from reduced model
  loglike = updateVars$logLike  #saved loglikelihoods from reduced model
  mvSamplesEst =  updateVars$mvSamplesEst

   inits <- as.list(target)
   names(inits) <- target
   for(i in 1:length(target)){
     expandTarget <- model$expandNodeNames(target[i])
     inits[[target[i]]] <- c(postReducedMCMC$summary[[chain.iter]][expandTarget, 'Mean'])
   }

   model$setInits(inits)

  #create new model for weights
  estimationModel <- model$newModel(replicate = TRUE)

  message("Building particle filter for model")

  if(is.null(pfType)){
    particleFilter <- myphdthesis::buildBootstrapFilterUpdate(model,
                                                              latent,
                                                              mvWSamplesWTSaved = weights,
                                                              mvWSamplesXSaved = unweightedLatentSamples,
                                                              mvEWSamplesXSaved = weightedLatentSamples,
                                                              logLikeVals = loglike,
                                                              mvSamplesEst = mvSamplesEst,
                                                              target = target,
                                                              control = pfControl)

    particleFilterEst <- myphdthesis::buildBootstrapFilterUpdate(estimationModel,
                                                                 latent,
                                                                 mvWSamplesWTSaved = weights,
                                                                 mvWSamplesXSaved = unweightedLatentSamples,
                                                                 mvEWSamplesXSaved = weightedLatentSamples,
                                                                 logLikeVals = loglike,
                                                                 mvSamplesEst = mvSamplesEst,
                                                                 target = target,
                                                                 control = pfControl)
  }


  if(!is.null(pfType)){
   if(!pfType %in% c("auxiliary", "bootstrap")) stop("Function currently works for auxiliary and bootstap Particle filters")
   if(pfType == "bootstrap"){
     particleFilter <- myphdthesis::buildBootstrapFilterUpdate(model,
                                                              latent,
                                                              target = target,
                                                             mvWSamplesWTSaved = weights,
                                                              mvWSamplesXSaved = unweightedLatentSamples,
                                                              mvEWSamplesXSaved = weightedLatentSamples,
                                                              logLikeVals = loglike,
                                                             mvSamplesEst = mvSamplesEst,
                                                              control = pfControl)

    particleFilterEst <- myphdthesis::buildBootstrapFilterUpdate(estimationModel,
                                                                    latent,
                                                                 target = target,
                                                                    mvWSamplesWTSaved = weights,
                                                                    mvWSamplesXSaved = unweightedLatentSamples,
                                                                    mvEWSamplesXSaved = weightedLatentSamples,
                                                                    logLikeVals = loglike,
                                                                 mvSamplesEst = mvSamplesEst,
                                                                    control = pfControl)
   }else{
     particleFilter <- myphdthesis::buildAuxiliaryFilterUpdate(model,
                                                              latent,
                                                               target = target,
                                                               mvWSamplesWTSaved = weights,
                                                               mvWSamplesXSaved = unweightedLatentSamples,
                                                               mvEWSamplesXSaved = weightedLatentSamples,
                                                               logLikeVals = loglike,
                                                               mvSamplesEst = mvSamplesEst,
                                                               control = pfControl)

     particleFilterEst <- myphdthesis::buildAuxiliaryFilterUpdate(estimationModel,
                                                                  latent,
                                                                  target = target,
                                                                  mvWSamplesWTSaved = weights,
                                                                  mvWSamplesXSaved = unweightedLatentSamples,
                                                                  mvEWSamplesXSaved = weightedLatentSamples,
                                                                  logLikeVals = loglike,
                                                                  mvSamplesEst = mvSamplesEst,
                                                                  control = pfControl)
   }
   }

  # message("Compiling the particle filter")
  # #compiling the model
  # compiledParticleFilter <- nimble::compileNimble(model,  particleFilter, showCompilerOutput = FALSE, resetFunctions = FALSE)
  #
  # #Log likelihood of last run and the Effective sample sizes
  # message("Running the particle filter")
  # vals <- c(rep(0, length(model$expandNodeNames(target, returnScalarComponents = TRUE))))
  # logLik <-   compiledParticleFilter$particleFilter$run(m = nParFiltRun, iterRun = 1, storeModelValues =  vals)
  # ESS <-   compiledParticleFilter$particleFilter$returnESS()

  message("Setting up the MCMC Configuration")
  #newModel <- model$newModel(replicate = TRUE)
  modelMCMCconf <- nimble::configureMCMC(model, nodes = NULL)

  if(is.null(pfType)){
    pfTypeUpdate = 'bootstrapUpdate'
  }else{
if(pfType == "bootstrap"){
  pfTypeUpdate = 'bootstrapUpdate'
  }else{
  if(pfType == "auxiliary") {
    pfTypeUpdate = 'auxiliaryUpdate'
    }
  }
  }


  modelMCMCconf$addSampler(target = target,
                           type = 'RW_PF_blockUpdate',
                           control = list(latents = latent,
                                          target = target,
                                          adaptInterval = n.iter,
                                          pfControl = list(saveAll = TRUE, M = M, iNodePrev = iNodePrev),
                                          pfNparticles = nParFiltRun,
                                          pfType = pfTypeUpdate,
                                          weights = weights,
                                          unweightedSamples= unweightedLatentSamples,
                                          weightedSamples = weightedLatentSamples,
                                          mvSamplesEst = mvSamplesEst,
                                          logLikeVals = loglike))

  modelMCMCconf$addMonitors(additionalPars)

  message("Building and compiling the PF MCMC")
  ## build and compile pMCMC sampler
  modelMCMC <- buildMCMC(modelMCMCconf)
  compiledList <- nimble::compileNimble(model,
                                        modelMCMC,
                                        resetFunctions = TRUE)
  timeStart2 <- Sys.time()
  message("Running the PF MCMC")
  #run MCMC
  mcmc.out <- nimble::runMCMC(compiledList$modelMCMC,
                              niter = n.iter,
                              nchains = 1,
                              nburnin = n.burnin,
                              #inits = initsList,
                              thin = n.thin,
                              setSeed = TRUE,
                              samples=TRUE,
                              samplesAsCodaMCMC = TRUE,
                              summary = TRUE,
                              WAIC = FALSE)
  timeEnd <- Sys.time()

  timetaken1 <- timeEnd - timeStart1
  timetaken2 <- timeEnd - timeStart2

  retList <- list()
  retList$timeRun <- timetaken2
  retList$samples <- mcmc.out$samples
  return(retList)
  })

  #set names for samples frpm various chains
  names(samplesList)  <- paste0('chain', 1:n.chains)

  # Time taken for each chain to run
  timetakenRun <- lapply(samplesList, function(x){x[[1]]})
  timetakenRun$all.chains <- sum(do.call('c', timetakenRun))

  #as samplesAsCodaMCMC
  samplesList1 <- coda::as.mcmc.list(lapply(samplesList, function(x){coda::as.mcmc(x[[2]])}))

  #Estimating summary
  summaryObject <- lapply(samplesList1, samplesSummary)
  names(summaryObject) <- paste0('chain', 1:n.chains)
  summaryObject$all.chains <- samplesSummary(do.call('rbind', samplesList1))



  #posteriorEstimates <- mcmc.out$summary$all.chains[target, 'Mean']

#   expandTarget <- model$expandNodeNames(target)
#   for(i in 1:length(expandTarget)){
#     estimationModel[[expandTarget[i]]] <- mcmc.out$summary[expandTarget[i], 'Mean']
#       #mcmc.out$summary$all.chains[expandTarget[i], 'Mean']
#   }
#
#   message("Compiling the estimation particle filter")
#   #compiling the model
#   compiledParticleFilterEst <- compileNimble(estimationModel,  particleFilterEst)
#
#   #Loglikelihood of last run and the Effective sample sizes
#   message("Running the estimation particle filter")
#   vals <- nimble::values(estimationModel, estimationModel$expandNodeNames(target, returnScalarComponents = TRUE))
#   logLik <-   compiledParticleFilterEst$particleFilterEst$run(m = nParFiltRun, storeModelValues = vals)
#   ESS <-   compiledParticleFilterEst$particleFilterEst$returnESS()
#
#
#   #save weights and samples
#   message("Extracting the weights and samples from posterior particle fiter")
#   if(is.null(pfType)) pfType = "bootstrap"
#   weights <- as.matrix(compiledParticleFilterEst$particleFilterEst$mvWSamples, "wts")
#   unweightedSamples <- as.matrix(compiledParticleFilterEst$particleFilterEst$mvWSamples, latent)
#   weightedSamples <- as.matrix( compiledParticleFilterEst$particleFilterEst$mvEWSamples, latent)
#   if(pfType == "auxiliary"){
#   logLike <- as.matrix(compiledParticleFilterEst$particleFilterEst$mvWSamples, "auxlog")
#   }else{
#     logLike <- as.matrix(compiledParticleFilterEst$particleFilterEst$mvWSamples, "bootLL")
# }


  message("Returning the results")

  #truth
  #list to return
  retList <- list()
retList$samples <- samplesList1
retList$summary <- summaryObject
retList$timeRun <- timetakenRun
  return(retList)
}


