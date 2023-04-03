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

#Model values for the weighted samples
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

#  message(paste("Estimating effective sample size"))
#  if(n.chains > 1){
#  posteriorEstimates <- mcmc.out$summary$all.chains[target, 'Mean']
#
#  expandTarget <- model$expandNodeNames(target)
#  for(i in 1:length(expandTarget)){
#    estimationModel[[expandTarget[i]]] <- mcmc.out$summary$all.chains[expandTarget[i], 'Mean']
#  }
#  }else{
#    posteriorEstimates <- mcmc.out$summary[target, 'Mean']
#
#    expandTarget <- model$expandNodeNames(target)
#    for(i in 1:length(expandTarget)){
#      estimationModel[[expandTarget[i]]] <- mcmc.out$summary[expandTarget[i], 'Mean']
#    }
# }
# # message("Compiling the estimation particle filter")
#  #compiling the model
#  compiledParticleFilterEst <- compileNimble(estimationModel,  particleFilterEst)
# #
# # #Loglikelihood of last run and the Effective sample sizes
#  message("Running the estimation particle filter")
#  logLik <-   compiledParticleFilterEst$particleFilterEst$run(m = nParFiltRun)
#  ESS <-   compiledParticleFilterEst$particleFilterEst$returnESS()
}

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

mvSamplesEst <- modelValues(modelValuesConf(vars = names,
                                            types = type,
                                            sizes = size))

#resize
 resize(mvSamplesEst, n.iter)

 message("Saving unsampled and sampled values in model values for updating")
for(iter in 1:n.iter){
  for(j in 1:length(names)){
    if(names[j] == latent & length(size[[1]]) > 1){
    mvSamplesEst[[names[j]]][[iter]] <- matrix(mcmcOut[iter, model$expandNodeNames(names[j])], nrow = size[[1]][1], ncol = size[[1]][2])
    }else{
      mvSamplesEst[[names[j]]][[iter]] <-  mcmcOut[iter, model$expandNodeNames(names[j])]
    }

  }
}


returnlist = list(
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

  mvSamplesEst =  updateVars$mvSamplesEst #saved weights

  #set initial values to the posterior mean of the saved targets and latent states
  if(!is.null(postReducedMCMC$summary)){
   inits <- as.list(target)
   names(inits) <- target
   for(i in 1:length(target)){
     expandTarget <- model$expandNodeNames(target[i])
     inits[[target[i]]] <- c(postReducedMCMC$summary[[chain.iter]][expandTarget, 'Mean'])
   }

   model$setInits(inits)
  }

  #create new model for weights
  estimationModel <- model$newModel(replicate = TRUE)

  message("Building particle filter for model")

  if(is.null(pfType)){
    particleFilter <- myphdthesis::buildBootstrapFilterUpdate(model,
                                                              latent,
                                                              mvSamplesEst = mvSamplesEst,
                                                              target = target,
                                                              control = pfControl)

    particleFilterEst <- myphdthesis::buildBootstrapFilterUpdate(estimationModel,
                                                                 latent,
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
                                                             mvSamplesEst = mvSamplesEst,
                                                              control = pfControl)

    particleFilterEst <- myphdthesis::buildBootstrapFilterUpdate(estimationModel,
                                                                    latent,
                                                                 target = target,
                                                                 mvSamplesEst = mvSamplesEst,
                                                                    control = pfControl)
   }else{
     particleFilter <- myphdthesis::buildAuxiliaryFilterUpdate(model,
                                                              latent,
                                                               target = target,
                                                               mvSamplesEst = mvSamplesEst,
                                                               control = pfControl)

     particleFilterEst <- myphdthesis::buildAuxiliaryFilterUpdate(estimationModel,
                                                                  latent,
                                                                  target = target,
                                                                  mvSamplesEst = mvSamplesEst,
                                                                  control = pfControl)
   }
   }

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

  message("Returning the results")

  #truth
  #list to return
  retList <- list()
retList$samples <- samplesList1
retList$summary <- summaryObject
retList$timeRun <- timetakenRun
  return(retList)
}


