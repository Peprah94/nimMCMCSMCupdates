# Functions to run the reduced and updated models

baselineSpartaEstimation <- function(model, #nimbleModel
                                     MCMCconfiguration = NULL, #configuration for MCMC
                                     pfType = "bootstrap",#Either 'auxiliary' or 'bootstrap'. Defaults to auxiliary
                                     pfControl = list(saveAll = TRUE,
                                                      #lookahead = "mean",
                                                      smoothing = FALSE), #list of controls for particle filter
                                     nParFiltRun = NULL, #Number of PF runs
                                     propCov = 'identity',
                                     latent #the latent variable
){

  timeStart1 <- Sys.time()
  target = MCMCconfiguration[["target"]]
  additionalPars = MCMCconfiguration[["additionalPars"]] #other dependent variables you seek to monitor
  n.iter = MCMCconfiguration[["n.iter"]]
  n.chains = MCMCconfiguration[["n.chains"]]
  n.burnin = MCMCconfiguration[["n.burnin"]]
  n.thin = MCMCconfiguration[["n.thin"]]
  if(is.null(nParFiltRun)) nParFiltRun = 5000

  #create new model for weights
  estimationModel <- model$newModel(replicate = TRUE)

  message("Building particle filter for model")
  if(is.null(pfType)) pfType <- "bootstrap" #set bootstrap as default
 # if(!is.null(pfType)){
    if(!pfType %in% c("auxiliary", "bootstrap")) stop("Function currently works for auxiliary and bootstap Particle filters")
    #particleFilter is used to run the MCMC
    #particleFilterEst is used for returning the weights at the posterior
    #values of the top-level nodes
    if(pfType == "auxiliary"){
      particleFilter <- nimbleSMC::buildAuxiliaryFilter(model,
                                                            nodes =  latent,
                                                             control = pfControl)

      particleFilterEst <- nimbleSMC::buildAuxiliaryFilter(estimationModel,
                                                               nodes =  latent,
                                                                control = pfControl)
    }

    if(pfType == "bootstrap"){
      particleFilter <- nimbleSMC::buildBootstrapFilter(model,
                                                           nodes =  latent,
                                                             control = pfControl)

      particleFilterEst <-  nimbleSMC::buildBootstrapFilter(estimationModel,
                                                               nodes =   latent,
                                                                 control = pfControl
                                                            )
    }
  #}else{
    # set default as bootstrap
    # particleFilter <- nimbleSMC::buildBootstrapFilter(model,
    #                                                       nodes =  latent,
    #                                                        control = pfControl)
    #
    # particleFilterEst <- nimbleSMC::buildBootstrapFilter(estimationModel,
    #                                                           nodes = latent,
    #                                                           control = pfControl)
  #}

  message("Compiling the particle filter")
  #compiling the model
  compiledParticleFilter <- compileNimble(model,
                                          particleFilter)

  #Loglikelihood of last run and the Effective sample sizes
  #message("Running the particle filter")
  compiledParticleFilterEst <- compileNimble(estimationModel,
                                             particleFilterEst)
  logLik <-   compiledParticleFilterEst$particleFilterEst$run(m = nParFiltRun)
  ESS <-   compiledParticleFilterEst$particleFilterEst$returnESS()

  message("Setting up the MCMC Configuration")
  #model <- model$newModel(replicate = TRUE)
  modelMCMCconf <- nimble::configureMCMC(model, nodes = NULL)

  if(is.null(pfType)) pfType = "bootstrap" #set bootstrap as default pfType

  modelMCMCconf$addSampler(target = target,
                           type = 'RW_PF_block',
                           control = list(latents = latent,
                                          pfControl = pfControl,
                                          pfNparticles = nParFiltRun,
                                          propCov = propCov,
                                          pfType = pfType))

  modelMCMCconf$addMonitors(additionalPars)

  message("Building and compiling the PF MCMC")
  ## build and compile pMCMC sampler
  modelMCMC <- nimble::buildMCMC(modelMCMCconf)
  compiledList <- nimble::compileNimble(model,
                                        modelMCMC,
                                        resetFunctions = FALSE)

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



# This function runs fits model with either particle filtering algorithms
# or mcmc. Can be used to fit both the reduced and baseline (full) models.

#' Fit a full or reduced state-space model using either MCMC or SMC methods
#'
#' @description Fit a given NIMBLE state space model using either with particle MCMC using nimbleSMC or mcmc using NIMBLE.
#'
#' @param model A NIMBLE model object, typically representing a state space model or a hidden Markov model.
#' @param MCMCconfiguration  A list specifying the parameterisation for nimble's runMCMC function.
#' @param pfType  The type of particle filter algorithm to fit. Defaults to 'bootstrap'.
#' @param pfControl  A list specifying different control options for the particle filter. Options are described in the details section below.
#' @param nParFiltRun  Number of particles for the particle filtering algorithm.
#' @param latent  A character specifying which variable in the state-space model is the latent variable.
#' @param block  A logical value indicating whether to block the model parameters and assign a random-walk block sampler.
#' @param mcmc  A logical value indicating whether to fit the model with 'mcmc' or 'particle filter'. Defaults to TRUE.
#' @author  Kwaku Peprah Adjei
#' @details
#' Each of the \code{pfControl()} list options are described in detail here:
#' \describe{
#' \item{lookahead}{The lookahead function used to calculate auxiliary weights.  Can choose between \code{'mean'} and \code{'simulate'}.
#'  Defaults to \code{'simulate'}.}
#'  \item{resamplingMethod}{The type of resampling algorithm to be used within the particle filter.  Can choose between \code{'default'} (which uses NIMBLE's \code{rankSample()} function), \code{'systematic'}, \code{'stratified'}, \code{'residual'}, and \code{'multinomial'}.  Defaults to \code{'default'}. Resampling methods other than \code{'default'} are currently experimental.}
#'  \item{saveAll}{Indicates whether to save state samples for all time points (\code{TRUE}), or only for the most recent time point (\code{FALSE})}
#'  \item{smoothing}{Decides whether to save smoothed estimates of latent states, i.e., samples from f(x[1:t]|y[1:t]) if \code{smoothing = TRUE}, or instead to save filtered samples from f(x[t]|y[1:t]) if \code{smoothing = FALSE}. \code{smoothing = TRUE} only works if \code{saveAll = TRUE}.}
#'  \item{timeIndex}{An integer used to manually specify which dimension of the latent state variable indexes time. This need only be set if the number of time points is less than or equal to the size of the latent state at each time point.}
#'  \item{initModel}{A logical value indicating whether to initialize the model before running the filtering algorithm.  Defaults to TRUE.}
#'  }
#'  \item{iNodePrev}{An integer specifying the number of years used to fit the reduced model.}
#' }
#'
#'  This function provides an implementation to fit reduced models that can be used to
#'  update the model parameters and latent state distribution in the future i.e. fit the reduced
#'  model. It can also be used to fit the full model that can be compared to the
#'  reduced and updated models
#'
#' @export
#'
#' @family smc update methods
#'
#' @examples
#' ## For illustration only.
#' stateSpaceCode <- nimbleCode({
#'   x0 ~ dnorm(0, var = 1)
#'   x[1] ~ dnorm(a* x0, var = 1)
#'   y[1] ~ dnorm(x[1]*c, var = .5)
#'   for(t in 2:10){
#'     x[t] ~ dnorm(a * x[t-1], var = 1)
#'     y[t] ~ dnorm(x[t]*c, var = .5)
#'   }
#' })
#'
#' model <- nimbleModel(code = exampleCode, data = list(y = rnorm(10)),
#'                      inits = list(x0 = 0, x = rnorm(10)))
#'
#'  spartaNimWeights(model,
#'  MCMCconfiguration = list(target = c('a', 'c'),
#' n.iter = 1000,
#' n.chains = 2,
#' n.burnin = 100,
#' n.thin = 2),
#' latent = "x",
#' mcmc = TRUE)
#'
#'

spartaNimWeights <- function(model, #nimbleModel
                             MCMCconfiguration = NULL, #configuration for MCMC
                             pfType = "bootstrap",#Either 'auxiliary' or 'bootstrap'. Defaults to bootstrap
                             pfControl = list(saveAll = TRUE,
                                              smoothing = FALSE), #list of controls for particle filter
                             nParFiltRun = NULL, #Number of PF runs
                             latent, #the latent variable
                             block = FALSE,
                             mcmc = TRUE # logical: if model should be fitted using MCMC or not. Default is TRUE
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
    #compile model
    nMCompile <- compileNimble(model)

    # configure MCMC
    cMCMC <- configureMCMC(model,
                           monitors = c(target, latent, additionalPars))#,
                           #useConjugacy = TRUE)

    # If we want to assign a block sampler to the target variables
    if(block == TRUE){
      cMCMC$removeSampler(target)
      cMCMC$addSampler(target, type = "RW_block")
    }else{
      cMCMC = cMCMC
    }

    # build, compile and run MCMC
    bMCMC <- buildMCMC(cMCMC)#,
                       #useConjugacy = FALSE)

    coMCMC <- compileNimble(bMCMC,
                            project = nMCompile)

    timeStart2 <- Sys.time()



    postburninTimeResult <- system.time({
    mcmc.out <- runMCMC(coMCMC,
                            niter = n.iter,
                            nchains = n.chains,
                            nburnin = n.burnin,
                            thin = n.thin,
                            setSeed = TRUE,
                            samples=TRUE,
                            samplesAsCodaMCMC = TRUE,
                            summary = TRUE,
                            WAIC = FALSE)
    })

    postburninTime <- postburninTimeResult[3]

    timeEnd <- Sys.time()

    timetaken1 <- timeEnd - timeStart1
    timetaken2 <- timeEnd - timeStart2

  }else{
    #fit model using particle filter
  if(is.null(nParFiltRun)) nParFiltRun = 5000

  #create new model for weights
  estimationModel <- model$newModel(replicate = TRUE)

  message("Building particle filter for model")
  #if(!is.null(pfType)){
  if(is.null(pfType)) pfType <- "bootstrap" # defaults to boostrap PF
    if(!pfType %in% c("auxiliary", "bootstrap")) stop("Function currently works for auxiliary and bootstap Particle filters")
#particleFilter is used to run the MCMC
#particleFilterEst is used for testing
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
  postburninTimeResult <- system.time({
  mcmc.out <- nimble::runMCMC(compiledList$modelMCMC,
                              niter = n.iter,
                              nchains = n.chains,
                              nburnin = n.burnin,
                              inits = initsList,
                              thin = n.thin,
                              setSeed = TRUE,
                              samples=TRUE,
                              samplesAsCodaMCMC = TRUE,
                              summary = TRUE,
                              WAIC = FALSE)
  })

  postburninTime <- postburninTimeResult[3]


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
  retList$timeRun <-  postburninTime
  return(retList)
}




##################################
#save weights and samples
#message("Extracting the weights and samples from particle fiter")

# This function creates a modelValue object for the MCMC samples.
# This is needed to fit the updated models

#' Create a modelValue object for the MCMC samples
#'
#' @description Create a modelValue object for MCMC object.
#'
#' @param model A NIMBLE model object, typically representing a state space model or a hidden Markov model for the updated model we want to fit.
#' @param reducedModel  A NIMBLE model object, typically representing a state space model or a hidden Markov model for the reduced model we fitted.
#' @param mcmcOut  The MCMC object of class mcmc.
#' @param iNodeAll  A logical value indicating whether we want to copy all the years.
#' @param target  A character specifying which variable is the target variable.
#' @param latent  A character specifying which variable in the state-space model is the latent variable.
#' @param n.iter  A integer specifying the number of iterations we want to use for the updated model.
#' @param m  An integer specifying the number of particles for the particle filter.
#' @param timeIndex  An integer specifying which dimension of the data is the time component.
#' @author  Kwaku Peprah Adjei
#'
#'  This function provides an implementation to fit reduced models that can be used to
#'  update the model parameters and latent state distribution in the future i.e. fit the reduced
#'  model. It can also be used to fit the full model that can be compared to the
#'  reduced and updated models
#'
#' @export
#'
#' @family smc update methods
#'
#' @examples
#' ## For illustration only.
#' stateSpaceCode <- nimbleCode({
#'   x0 ~ dnorm(0, var = 1)
#'   x[1] ~ dnorm(a* x0, var = 1)
#'   y[1] ~ dnorm(x[1]*c, var = .5)
#'   for(t in 2:10){
#'     x[t] ~ dnorm(a * x[t-1], var = 1)
#'     y[t] ~ dnorm(x[t]*c, var = .5)
#'   }
#' })
#'
#' model <- nimbleModel(code = exampleCode, data = list(y = rnorm(10)),
#'                      inits = list(x0 = 0, x = rnorm(10)))
#'
#'  spartaNimWeights(model,
#'  MCMCconfiguration = list(target = c('a', 'c'),
#' n.iter = 1000,
#' n.chains = 2,
#' n.burnin = 100,
#' n.thin = 2),
#' latent = "x",
#' mcmc = TRUE)
#'
#'
updateUtils <- function(model, #updated model
                        reducedModel, #reduced model
                        mcmcOut, # mcmc object
                        iNodeAll, # all nodes
                        latent,  #latent variable
                        target, # target parameters
                        n.iter, # number of iteration
                        m, # number of particles
                        timeIndex){ #time index for particle filter

# iNodeAll specifies whether the we want to just use time t to update the model or all the data point.
# if iNodeAll = TRUE, we use all the time points.
if(iNodeAll == FALSE){

# get the latent nodes
latentNodes <- model$expandNodeNames(latent)
nodes <- findLatentNodes(model, latent, timeIndex)

#get the last node
lastNodes <- findLatentNodes(reducedModel, latent, timeIndex )
lastNode <- lastNodes[length(lastNodes)]

#get the dimensions of the parameters
dims <- lapply(nodes[1:2], function(n) nimDim(model[[n]]))

# get the target variable
target <- target[!target %in% latent]
if(length(unique(dims)) > 1)
  stop('sizes or dimensions of latent states varies')
vars <- c(model$getVarNames(nodes =  nodes), target)

#create model object
modelSymbolObjects <- model$getSymbolTable()$getSymbolObjects()[vars]

names <- sapply(modelSymbolObjects, function(x)return(x$name))
type <- sapply(modelSymbolObjects, function(x)return(x$type))

# Get the size of the modelValue object
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

# create the modelValue object
mvSamplesEst <- modelValues(modelValuesConf(vars = names,
                                            types = type,
                                            sizes = size))

# resize the modelValue object
 resize(mvSamplesEst, n.iter)

 message("Saving unsampled and sampled values in model values for updating")
for(iter in 1:n.iter){
  for(j in 1:length(names)){
    if(names[j] == latent & length(size[[1]]) > 1){
      # Obtain the length of extra target nodes that needs to be simulated
      extraVars <- length(model$expandNodeNames(names[j])) == length(reducedModel$expandNodeNames(names[j]))
      if(!extraVars){
      extraValues <- values(model, model$expandNodeNames(nodes[-1]))
        #rep(0, length(model$expandNodeNames(nodes[-1])))
      #names(extraValues) <- extraVars
      allVals <- c(mcmcOut[iter, reducedModel$expandNodeNames(lastNode)], extraValues)
      names(allVals) <- model$expandNodeNames(nodes)
    mvSamplesEst[[names[j]]][[iter]] <- matrix(allVals, nrow = size[[1]][1], ncol = size[[1]][2])
     }else{

        allVals <- mcmcOut[iter, reducedModel$expandNodeNames(names[j])]
        mvSamplesEst[[names[j]]][[iter]] <- matrix(allVals, nrow = size[[1]][1], ncol = size[[1]][2])
      }
    }else{
      extraVars <- length(model$expandNodeNames(names[j])) == length(reducedModel$expandNodeNames(names[j]))
      #print(extraVars)
     if(!extraVars){
      if(names[j] == latent){
      estValues <- c(mcmcOut[iter, reducedModel$expandNodeNames(lastNode)], rep(0, length(model$expandNodeNames(names[j]))-1))
      names(estValues) <- model$expandNodeNames(names[j])
      }else{
        namesExpanded <- reducedModel$expandNodeNames(names[j])
        lastName <- namesExpanded[length(namesExpanded)]
        extraVals <- values(model, (model$expandNodeNames(names[j]))[-1])
       estValues <- c(mcmcOut[iter, lastName ], extraVals)#rep(0, length(model$expandNodeNames(names[j]))-1))
       names(estValues) <- model$expandNodeNames(names[j])
      }
       }else{
       estValues <- c(mcmcOut[iter, model$expandNodeNames(names[j])])
     }
      mvSamplesEst[[names[j]]][[iter]] <-  estValues
    }

  }
}


returnlist = list(
  mvSamplesEst = mvSamplesEst
)
}else{
# get nodes
  latentNodes <- model$expandNodeNames(latent)
  nodes <- findLatentNodes(model, latent, timeIndex)

# Get latent nodes and extra nodes
  lastNodes <- findLatentNodes(reducedModel, latent, timeIndex )
  lastNode <- lastNodes[length(lastNodes)]
  extraNodes <- nodes[!nodes %in% lastNodes]
  extraNodes <- model$expandNodeNames(extraNodes)

  # Get the dimensions of latent nodes
dims <- lapply(nodes[1:2], function(n) nimDim(model[[n]]))

# target variables
  target <- target[!target %in% latent]

  # Dimensions of variables
  if(length(unique(dims)) > 1)
    stop('sizes or dimensions of latent states varies')
  vars <- c(model$getVarNames(nodes =  nodes), model$getVarNames(nodes = target))

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

  # Create modelValue object
  mvSamplesEst <- modelValues(modelValuesConf(vars = names,
                                              types = type,
                                              sizes = size))
  resize(mvSamplesEst, n.iter) #resize modelvalue object

  message("Saving unsampled and sampled values in model values for updating")
  for(iter in 1:n.iter){
    for(j in 1:length(names)){
      if(names[j] == latent & length(size[[1]]) > 1){
        # Extra target variables to simulate in sampler
        extraVars <- length(model$expandNodeNames(names[j])) == length(reducedModel$expandNodeNames(names[j]))
        # extraVars <- model$expandNodeNames(names[j])[!model$expandNodeNames(names[j]) %in% reducedModel$expandNodeNames(names[j])]
        if(!extraVars){
          #if there are extra variables, initialise the model with the initial values
          extraValues <- nimble::values(model, extraNodes)#rep(0, length(extraNodes))
            #values(model, model$expandNodeNames(nodes[-1]))
          #rep(0, length(model$expandNodeNames(nodes[-1])))
          #names(extraValues) <- extraVars
          allVals <- c(mcmcOut[iter, model$expandNodeNames(lastNodes)], extraValues)
          names(allVals) <- model$expandNodeNames(nodes)
          mvSamplesEst[[names[j]]][[iter]] <- matrix(allVals, nrow = size[[1]][1], ncol = size[[1]][2])
        }else{
#if not, then save all the target values
          allVals <- mcmcOut[iter, reducedModel$expandNodeNames(names[j])]
          mvSamplesEst[[names[j]]][[iter]] <- matrix(allVals, nrow = size[[1]][1], ncol = size[[1]][2])
        }
      }else{
        extraVars <- length(model$expandNodeNames(names[j])) == length(reducedModel$expandNodeNames(names[j]))
        #print(extraVars)
        if(!extraVars){
          if(names[j] == latent){
            if(names[j] == "N"){
              meanN <- mean(mcmcOut[iter, lastNodes])
              estValues <- c(mcmcOut[iter, lastNodes], rep(meanN, length(extraNodes)))
            }else{
              estValues <- c(mcmcOut[iter, lastNodes], rep(0, length(extraNodes)))
            }

            names(estValues) <- model$expandNodeNames(names[j])
          }else{
            namesExpanded <- reducedModel$expandNodeNames(names[j])
            extraValNodes <- model$expandNodeNames(names[j])[!model$expandNodeNames(names[j]) %in% namesExpanded]
            #lastName <- namesExpanded[length(namesExpanded)]
            extraVals <- nimble::values(model, extraValNodes)#rep(0, (length(model$expandNodeNames(names[j])) - length(namesExpanded)))#values(model, (model$expandNodeNames(names[j]))[-1])
            estValues <- c(mcmcOut[iter, namesExpanded ], extraVals)#rep(0, length(model$expandNodeNames(names[j]))-1))
            names(estValues) <- model$expandNodeNames(names[j])
          }
        }else{
          estValues <- c(mcmcOut[iter, model$expandNodeNames(names[j])])
        }
        mvSamplesEst[[names[j]]][[iter]] <-  estValues
      }

    }
  }

#save results
  returnlist = list(
    mvSamplesEst = mvSamplesEst
  )

}
return(returnlist)
}

###########################################
# This function runs fits model with either particle filtering algorithms
# or mcmc. Can be used to fit both the reduced and baseline (full) models.

#' Fit a full or reduced state-space model using either MCMC or SMC methods
#'
#' @description Fit a given updated NIMBLE state space model.
#'
#' @param model A NIMBLE model object, typically representing a state space model or a hidden Markov model.
#' @param reducedModel A NIMBLE model object, typically representing a state space model or a hidden Markov model for the reduced method.
#' @param MCMCconfiguration  A list specifying the parameterisation for nimble's runMCMC function.
#' @param pfType  The type of particle filter algorithm to fit. Defualts to 'bootstrap'.
#' @param pfControl  A list specifying different control options for the particle filter. Options are described in the details section below.
#' @param nParFiltRun  Number of particles for the particle filtering algorithm.
#' @param latent  A character specifying which variable in the state-space model is the latent variable.
#' @param postReducedMCMC MCMC samples from the reduced model
#' @param iNodeAll  A logical value indicating whether we want to copy all the years.
#' @param target  A character specifying which variable is the target variable.
#' @param extraVars  A character specifying which variable is the target variable.
#' @param mcmcScale  A logical value indicating whether to fit the model with 'mcmc' or 'particle filter'. Defaults to TRUE.
#' @param propCov  A character specifying which variable is the target variable.
#' @param propCov1  A character specifying which variable is the target variable.
#' @author  Kwaku Peprah Adjei
#' @details
#' Each of the \code{pfControl()} list options are described in detail here:
#' \describe{
#' \item{lookahead}{The lookahead function used to calculate auxiliary weights.  Can choose between \code{'mean'} and \code{'simulate'}.
#'  Defaults to \code{'simulate'}.}
#'  \item{resamplingMethod}{The type of resampling algorithm to be used within the particle filter.  Can choose between \code{'default'} (which uses NIMBLE's \code{rankSample()} function), \code{'systematic'}, \code{'stratified'}, \code{'residual'}, and \code{'multinomial'}.  Defaults to \code{'default'}. Resampling methods other than \code{'default'} are currently experimental.}
#'  \item{saveAll}{Indicates whether to save state samples for all time points (\code{TRUE}), or only for the most recent time point (\code{FALSE})}
#'  \item{smoothing}{Decides whether to save smoothed estimates of latent states, i.e., samples from f(x[1:t]|y[1:t]) if \code{smoothing = TRUE}, or instead to save filtered samples from f(x[t]|y[1:t]) if \code{smoothing = FALSE}. \code{smoothing = TRUE} only works if \code{saveAll = TRUE}.}
#'  \item{timeIndex}{An integer used to manually specify which dimension of the latent state variable indexes time. This need only be set if the number of time points is less than or equal to the size of the latent state at each time point.}
#'  \item{initModel}{A logical value indicating whether to initialize the model before running the filtering algorithm.  Defaults to TRUE.}
#'  }
#'  \item{iNodePrev}{An integer specifying the number of years used to fit the reduced model.}
#' }
#'
#'  This function provides an implementation to fit reduced models that can be used to
#'  update the model parameters and latent state distribution in the future i.e. fit the reduced
#'  model. It can also be used to fit the full model that can be compared to the
#'  reduced and updated models
#'
#' @export
#'
#' @family smc update methods
#'
#' @examples
#' ## For illustration only.
#' stateSpaceCode <- nimbleCode({
#'   x0 ~ dnorm(0, var = 1)
#'   x[1] ~ dnorm(a* x0, var = 1)
#'   y[1] ~ dnorm(x[1]*c, var = .5)
#'   for(t in 2:10){
#'     x[t] ~ dnorm(a * x[t-1], var = 1)
#'     y[t] ~ dnorm(x[t]*c, var = .5)
#'   }
#' })
#'
#' model <- nimbleModel(code = exampleCode, data = list(y = rnorm(10)),
#'                      inits = list(x0 = 0, x = rnorm(10)))
#'
#'  spartaNimUpdates(model,
#'  reducedModel,
#'  MCMCconfiguration = list(target = c('a', 'c'),
#' n.iter = 1000,
#' n.chains = 2,
#' n.burnin = 100,
#' n.thin = 2),
#' latent = "x",
#' mcmc = TRUE)
#'
#'
spartaNimUpdates <- function(model, #nimbleModel
                             reducedModel,
                             MCMCconfiguration = NULL, #configuration for MCMC
                             pfType = NULL,#Either 'auxiliary' or 'bootstrap'. Defaults to auxiliary
                             pfControl = NULL, #list of controls for particle filter
                             nParFiltRun = NULL, #Number of PF runs
                             #updatePFControl = list(iNodePrev = NULL, M = NULL),
                             latent, #the latent variable
                             #newData,
                             postReducedMCMC, #MCMC results from reduced model
                             iNodeAll = TRUE, # whether we want to use all the years from the reduced model in our updated model
                             target = NULL, # target variable
                             extraVars = NULL, #indicator of extra target variabøes
                             mcmcScale = 1,
                             propCov = 'identity',
                             propCov1 = 'identity',
                             adaptiveSampler = FALSE,
                             adaptScaleOnly = FALSE,
                             samplerType = "type1",
                             nCores = NULL
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

  # set default values for the number of particles and timeIndex
  if(is.null(nParFiltRun)) nParFiltRun = 5000
  if(is.null(timeIndex)) timeIndex = 1
  if(is.null(nCores)) nCores = n.chains
   if(is.null(samplerType)) samplerType <- "type1"

  #samplesList  <- vector('list', n.chains)

  if(nCores > 1){
    message(paste("Parallelising updated method using", nCores, "cores"))
  }else{
    message("Fitting updated method without parallelising")
  }

  runUpdatingFunction <- function(i,
                                  model, #nimbleModel
                                  reducedModel,
                                  MCMCconfiguration, #configuration for MCMC
                                  pfType,#Either 'auxiliary' or 'bootstrap'. Defaults to auxiliary
                                  pfControl, #list of controls for particle filter
                                  nParFiltRun, #Number of PF runs
                                  #updatePFControl = list(iNodePrev = NULL, M = NULL),
                                  latent, #the latent variable
                                  #newData,
                                  postReducedMCMC, #MCMC results from reduced model
                                  iNodeAll, # whether we want to use all the years from the reduced model in our updated model
                                  target, # target variable
                                  extraVars, #indicator of extra target variabøes
                                  mcmcScale,
                                  propCov,
                                  propCov1,
                                  adaptiveSampler,
                                  adaptScaleOnly,
                                  samplerType,
                                  nCores){

    target = MCMCconfiguration[["target"]]
    additionalPars = MCMCconfiguration[["additionalPars"]] #other dependent variables you seek to monitor
    n.iter = MCMCconfiguration[["n.iter"]]
    n.chains = MCMCconfiguration[["n.chains"]]
    n.burnin = MCMCconfiguration[["n.burnin"]]
    n.thin = MCMCconfiguration[["n.thin"]]
    iNodePrev = pfControl[["iNodePrev"]]
    M = pfControl[["M"]]
    timeIndex <- pfControl[["timeIndex"]]

    # set default values for the number of particles and timeIndex
    if(is.null(nParFiltRun)) nParFiltRun = 5000
    if(is.null(timeIndex)) timeIndex = 1
    if(is.null(nCores)) nCores = n.chains
    if(is.null(propCov)) propCov = "identity"
    if(is.null(propCov1)) propCov1 = "identity"

    if(nCores > 1) {
      dirName <- file.path(tempdir(), 'nimble_generatedCode', paste0("worker_", i))
    } else dirName = NULL

print(dirName)
      timeStart1 <- Sys.time()

      # create modelValues for MCMC output
      if(!is.null(extraVars)){
      updateVars <- updateUtils(model = model, #reduced model
                                reducedModel = reducedModel,
                                mcmcOut = postReducedMCMC$samples[[i]],
                                latent = latent,
                                target = c(target, additionalPars, extraVars),
                                n.iter = n.iter ,
                                m = nParFiltRun,
                                iNodeAll = TRUE,
                                timeIndex = timeIndex)
      } else {
        updateVars <- updateUtils(model = model, #reduced model
                                  reducedModel = reducedModel,
                                  mcmcOut = postReducedMCMC$samples[[i]],
                                  latent = latent,
                                  target = c(target, additionalPars),
                                  n.iter = n.iter ,
                                  m = nParFiltRun,
                                  iNodeAll = TRUE,
                                  timeIndex = timeIndex)
      }

      mvSamplesEst =  updateVars$mvSamplesEst #saved weights

      #set initial values to the posterior mean of the saved targets and latent states
      #if(!is.null(postReducedMCMC$summary)){
      # inits <- as.list(target)
      # names(inits) <- target
      # for(i in 1:length(target)){
      #   expandTarget <- model$expandNodeNames(target[i])
      #   inits[[target[i]]] <- c(postReducedMCMC$summary[[chain.iter]][expandTarget, 'Mean'])
      # }

      # model$setInits(inits)
      #}
      ## set options to make history accessible
      #nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
      #nimbleOptions(MCMCsaveHistory = TRUE)
      ## Next, set up and run your MCMC.
      ## Now access the history information:
      ## only for RW_block
      #create new model for weights
      estimationModel <- model$newModel(replicate = TRUE)

      message("Building particle filter for model")



      if(is.null(pfType)) pfType <- "bootstrap" #defaults to bootstrap PF

      if(!pfType %in% c("auxiliary", "bootstrap")) stop("Function currently works for auxiliary and bootstap Particle filters")
      if(pfType == "bootstrap"){
        particleFilter <- nimMCMCSMCupdates::buildBootstrapFilterUpdate(model,
                                                                        latent,
                                                                        target = target,
                                                                        mvSamplesEst = mvSamplesEst,
                                                                        control = pfControl)

        particleFilterEst <- nimMCMCSMCupdates::buildBootstrapFilterUpdate(estimationModel,
                                                                           latent,
                                                                           target = target,
                                                                           mvSamplesEst = mvSamplesEst,
                                                                           control = pfControl)
      }else{
        particleFilter <- nimMCMCSMCupdates::buildAuxiliaryFilterUpdate(model,
                                                                        nodes = latent,
                                                                        target = target,
                                                                        mvSamplesEst = mvSamplesEst,
                                                                        control = pfControl)

        particleFilterEst <- nimMCMCSMCupdates::buildAuxiliaryFilterUpdate(estimationModel,
                                                                           latent,
                                                                           target = target,
                                                                           mvSamplesEst = mvSamplesEst,
                                                                           control = pfControl)
      }


      #checking if the updated pF works very well
      #targetAsScalar <- estimationModel$expandNodeNames(target, returnScalarComponents = TRUE)
       #compiledParticleFilter <- compileNimble(estimationModel,  particleFilterEst)
#
      #compiledParticleFilter$particleFilterEst$run(m = 10, iterRun = 3, storeModelValues = values(estimationModel, targetAsScalar))
 #colMeans(compiledParticleFilter$particleFilterEst$mvWSamples[["z"]][[5]])

      message("Setting up the MCMC Configuration")
      #newModel <- model$newModel(replicate = TRUE)

      if("samplerTimes" %in% c(additionalPars)){
       # additionalPars <- additionalPars[!additionalPars %in% "samplerTimes"]
      modelMCMCconf <- nimble::configureMCMC(model,
                                             monitors = c(target, latent,
                                                          additionalPars[!additionalPars %in% "samplerTimes"]),
                                             monitors2 = c("samplerTimes"),
                                             nodes = NULL)#, monitors = c(target, latent, additionalPars))

      } else {
        modelMCMCconf <- nimble::configureMCMC(model,
                                               monitors = c(target, latent, additionalPars),
                                               nodes = NULL)
      }
    if(samplerType == "type2"){

      model <-  modelMCMCconf$model
      # set sampler to use
      if(pfType == "bootstrap"){
        pfTypeUpdate = 'bootstrapUpdate'
      }else{
        if(pfType == "auxiliary") {
          pfTypeUpdate = 'auxiliaryUpdate'
        }
      }
      targetSamplesRun <- lapply(model$getVarNames(nodes = target), function(x){
        ret <- model$getDependencies(x, stochOnly = TRUE,  downstream = TRUE, self = FALSE)
        nodesToReturn <- model$getVarNames(nodes = ret)
        return(nodesToReturn)
      })


      #################################
      # Last run: Those that depend on the observation
      #######################################

      lastRun <- sapply(targetSamplesRun, function(x){
        ret <- all(model$isData(x))
        return(ret)
        #any(x %in% c(model$expandNodeNames(nodes = latents), model$getVarNames(nodes = target) ))
      })

      lastRunPars <- model$expandNodeNames(model$getVarNames(nodes = target)[lastRun])
      print(lastRunPars)
      #####################################
      # Third run: Hyperparameters of the last run parameters
      # I will call it third run
      #####################################
      thirdRun <- sapply(targetSamplesRun, function(x){
        any(x %in% model$getVarNames(nodes = lastRunPars))
      })
      thirdRunPars <-model$expandNodeNames(model$getVarNames(nodes = target)[thirdRun])
      print(thirdRunPars)

      ############################################
      # Second run: The parameters that directly affect the latent state
      ####################################
      suppressWarnings(secondRun <- sapply(targetSamplesRun, function(x){
        ret <- x == c(model$getVarNames(nodes = latent), model$getVarNames(nodes = model$getDependencies(lastRunPars, self = FALSE, stochOnly = TRUE)))
        ret <- all(ret == TRUE)
        return(ret)
      }))
      secondRunPars <- model$expandNodeNames(model$getVarNames(nodes = target)[secondRun])
      print(secondRunPars)

      ############################################
      # First run: The parameters that affect second run
      ####################################
      firstRun <- sapply(targetSamplesRun, function(x){
        any(x %in% model$getVarNames(nodes = secondRunPars))
      })
      firstRunPars <- model$expandNodeNames(model$getVarNames(nodes = target)[firstRun])
      print(firstRunPars)


      #####################################
      # Set the conditions fot the MCMC run
      #####################################
      simLastPars <- TRUE
      simThirdPars <- TRUE

      if(length(lastRunPars) == 0) simLastPars <- FALSE
      if(length(thirdRunPars) == 0) simThirdPars <- FALSE


      if(length(firstRunPars) == 0){
        simFirstPars <- FALSE
      } else{
        simFirstPars <- TRUE
        topParams <- firstRunPars

      }
      if(length(secondRunPars) == 0){
        simSecondPars <- FALSE
      }else{
        simSecondPars <- TRUE
        #if(length(firstRunPars)==0){
         # topParamsInter <- firstRunPars
       # }else{
          topParamsInter <- secondRunPars
        #}
      }

      print(c(simFirstPars, simSecondPars,
              simThirdPars, simLastPars))

      #set extraVars
      targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
      extraVars <- targetAsScalar[!targetAsScalar %in% topParamsInter]

      if(simThirdPars){
        isDiscrete <- any(model$isDiscrete(thirdRunPars) == TRUE)
        if(isDiscrete){
          if(length(thirdRunPars) > 1){
            modelMCMCconf$addSampler(target = thirdRunPars,
                                     "AF_slice")
          } else {
            modelMCMCconf$addSampler(target = thirdRunPars,
                                     "slice")
          }
        } else {
          if(length(thirdRunPars) > 1){
            modelMCMCconf$addSampler(target = thirdRunPars,
                                     "RW_block")
          } else {
            modelMCMCconf$addSampler(target = thirdRunPars,
                                     "RW")
          }
        }
      }

      if(simLastPars){
        isDiscrete <- any(model$isDiscrete(lastRunPars) == TRUE)
        if(isDiscrete){
          if(length(lastRunPars) > 1){
            modelMCMCconf$addSampler(target = lastRunPars,
                                     "AF_slice")
          } else {
            modelMCMCconf$addSampler(target = lastRunPars,
                                     "slice")
          }
        } else {
          if(length(lastRunPars) > 1){
            modelMCMCconf$addSampler(target = lastRunPars,
                                     "RW_block",
                                     control = list(scale = mcmcScale,
                                                    adaptive = adaptiveSampler))
          } else {
            modelMCMCconf$addSampler(target = lastRunPars,
                                     "RW",
                                     control = list(scale = mcmcScale))
          }
        }
      }


      modelMCMCconf$addMonitors(additionalPars)

      # if dimfirstPars
      if(simFirstPars){
        modelMCMCconf$addSampler(target = topParams,
                                 type = 'RW_PF_blockUpdateV2',
                                 control = list(latents = latent,
                                                #target = target,
                                                adaptive = adaptiveSampler,
                                                propCov = propCov,
                                                propCov1 = propCov1,
                                                #adaptInterval = 100,
                                                scale = mcmcScale,
                                                # pf = particleFilter,
                                                pfControl = pfControl, #list( M = M, iNodePrev = iNodePrev),
                                                pfNparticles = nParFiltRun,
                                                pfType = pfTypeUpdate,
                                                postSamples = postReducedMCMC$samples[[i]],
                                                mvSamplesEst = mvSamplesEst,
                                                extraVars = extraVars,
                                                iNodePrev = iNodePrev,
                                                additionalPars = additionalPars,
                                                reducedModel = reducedModel)
        )
      }#else{

      #set up MCMC configuration
      modelMCMCconf$addSampler(target = topParamsInter,
                               type = 'RW_PF_blockUpdateV2',
                               control = list(latents = latent,
                                              #target = target,
                                              adaptive = adaptiveSampler,
                                              propCov = propCov,
                                              propCov1 = propCov1,
                                              #adaptInterval = 100,
                                              scale = mcmcScale,
                                              # pf = particleFilter,
                                              pfControl = pfControl, #list( M = M, iNodePrev = iNodePrev),
                                              pfNparticles = nParFiltRun,
                                              pfType = pfTypeUpdate,
                                              postSamples = postReducedMCMC$samples[[i]],
                                              mvSamplesEst = mvSamplesEst,
                                              extraVars = extraVars,
                                              iNodePrev = iNodePrev,
                                              additionalPars = additionalPars,
                                              reducedModel = reducedModel)
      )
#}



      modelMCMCconf$printSamplers(executionOrder = TRUE)

    } else if(samplerType == "type1") {
      # set sampler to use
      if(pfType == "bootstrap"){
        pfTypeUpdate = 'bootstrapUpdate'
      }else{
        if(pfType == "auxiliary") {
          pfTypeUpdate = 'auxiliaryUpdate'
        }
      }
      #set up MCMC configuration
      modelMCMCconf$addSampler(target = target,
                               type = 'RW_PF_blockUpdate',
                               control = list(latents = latent,
                                              #target = target,
                                              adaptive = adaptiveSampler,
                                              adaptScaleOnly = adaptScaleOnly,
                                              propCov = propCov,
                                              propCov1 = propCov1,
                                              #adaptInterval = 100,
                                              scale = mcmcScale,
                                              # pf = particleFilter,
                                              pfControl = pfControl, #list( M = M, iNodePrev = iNodePrev),
                                              pfNparticles = nParFiltRun,
                                              pfType = pfTypeUpdate,
                                              postSamples = postReducedMCMC$samples[[i]],
                                              mvSamplesEst = mvSamplesEst,
                                              extraVars = extraVars,
                                              iNodePrev = iNodePrev,
                                              additionalPars = additionalPars,
                                              reducedModel = reducedModel)
      )

}

      message("Building and compiling the PF MCMC")
      ## build and compile pMCMC sampler
      modelMCMC <- buildMCMC(modelMCMCconf)

      message("Building and compiling the PF MCMC")
      compiledList <- nimble::compileNimble(model,
                                            modelMCMC,
                                            resetFunctions = TRUE#,
                                            #dirName = NULL
                                            )
      timeStart2 <- Sys.time()
      message("Running the PF MCMC")
      # #run MCMC
      # compiledList$modelMCMC$run(niter = n.iter,
      #                            #nburnin = n.burnin,
      #                            time = TRUE)
     #
     # system.time(mcmc.out <- nimble::runMCMC(compiledList$modelMCMC,
     #                              niter = n.iter,
     #                              nchains = 1,
     #                              nburnin = n.burnin,
     #                              #inits = initsList,
     #                              thin = n.thin,
     #                              setSeed = TRUE,
     #                              samples= 123,
     #                              samplesAsCodaMCMC = TRUE,
     #                              summary = TRUE,
     #                              WAIC = FALSE
     #                              )
     # )

     postburninTimeResult <- system.time({
       mcmc.out <- nimble::runMCMC(compiledList$modelMCMC,
                                   niter = n.iter,
                                   nchains = 1,
                                   nburnin = n.burnin,
                                   #inits = initsList,
                                   thin = n.thin,
                                   setSeed = TRUE,
                                   samples= 123,
                                   samplesAsCodaMCMC = TRUE,
                                   summary = TRUE,
                                   WAIC = FALSE
       )
     })
     postburninTime <- postburninTimeResult[3]
      timeEnd <- Sys.time()

      timetaken1 <- timeEnd - timeStart1
      timetaken2 <- timeEnd - timeStart2
      if(getNimbleOption('MCMCsaveHistory')){
      scaleHistory <- compiledList$modelMCMC$samplerFunctions[[1]]$getScaleHistory()
      acceptanceHistory <- compiledList$modelMCMC$samplerFunctions[[1]]$getAcceptanceHistory()
      covarianceHistory <- compiledList$modelMCMC$samplerFunctions[[1]]$getPropCovHistory()
      } else {
        scaleHistory <- NA
        acceptanceHistory <- NA
        covarianceHistory <- NA
      }
#compiledList$modelMCMC$mvSamples
#compiledList$modelMCMC$mvSaved
      retList <- list()
      retList$timeRun <- timetaken2
      retList$samples <- mcmc.out$samples
      retList$elapsedTime <- postburninTime

if("samplerTimes" %in% c(additionalPars)) retList$samples2 <- mcmc.out$samples2
if(getNimbleOption('MCMCsaveHistory')){
retList$scaleHistory <- scaleHistory
retList$acceptanceHistory <- acceptanceHistory
retList$covarianceHistory <- covarianceHistory
}
      return(retList)
    }


  if(nCores > 1){
    samplesList <- parallel::mclapply(1:n.chains,
                                      runUpdatingFunction,
                                      model, #nimbleModel
                                      reducedModel,
                                      MCMCconfiguration,
                                      pfType,
                                      pfControl,
                                      nParFiltRun,
                                      #updatePFControl = list(iNodePrev = NULL, M = NULL),
                                      latent, #the latent variable
                                      #newData,
                                      postReducedMCMC, #MCMC results from reduced model
                                      iNodeAll, # whether we want to use all the years from the reduced model in our updated model
                                      target, # target variable
                                      extraVars, #indicator of extra target variabøes
                                      mcmcScale,
                                      propCov,
                                      propCov1,
                                      adaptiveSampler,
                                      adaptScaleOnly,
                                      samplerType,
                                      nCores = n.chains,
                                      mc.cores = nCores)
  } else{
    samplesList <- lapply(1:n.chains,
                          runUpdatingFunction,
                          model, #nimbleModel
                          reducedModel,
                          MCMCconfiguration,
                          pfType,
                          pfControl,
                          nParFiltRun,
                          #updatePFControl = list(iNodePrev = NULL, M = NULL),
                          latent, #the latent variable
                          #newData,
                          postReducedMCMC, #MCMC results from reduced model
                          iNodeAll, # whether we want to use all the years from the reduced model in our updated model
                          target, # target variable
                          extraVars, #indicator of extra target variabøes
                          mcmcScale,
                          propCov,
                          propCov1,
                          adaptiveSampler,
                          adaptScaleOnly,
                          samplerType,
                          nCores = 1)
  }

#########################################
# POST-PROCESS SAMPLES GENERATED
#########################################

  #set names for samples frpm various chains
  names(samplesList)  <- paste0('chain', 1:n.chains)

  # Time taken for each chain to run
  timetakenRun <- lapply(samplesList, function(x){x[[3]]})
  timetakenRun$all.chains <- sum(do.call('c', timetakenRun))

  # get sampler times for all
  if("samplerTimes" %in% c(additionalPars)){
    samplesList2 <- coda::as.mcmc.list(lapply(samplesList, function(x){coda::as.mcmc(x[[4]])}))
    #Estimating summary
    summaryObjectTimes <- lapply(samplesList2, samplesSummary)
    names(summaryObjectTimes) <- paste0('chain', 1:n.chains)
    summaryObjectTimes$all.chains <- samplesSummary(do.call('rbind', samplesList2))
  } else {
    samplesList2 <- NA
    summaryObjectTimes <- NA
  }


  # Post process the covariance and acceptance history
  if(getNimbleOption('MCMCsaveHistory')){
scaleHistory <- lapply(samplesList, function(x){x[[5]]})
scaleHistory$all.chains <- do.call("c", scaleHistory)

acceptanceHistory <- lapply(samplesList, function(x){x[[6]]})
acceptanceHistory$all.chains <- do.call("c", acceptanceHistory)

propCovHistory <- lapply(samplesList, function(x){x[[7]]})
}

  # Target samples of interest
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
retList$sampleTimes <- samplesList2
retList$summarySampleTimes <- summaryObjectTimes
if(getNimbleOption('MCMCsaveHistory')){
retList$scaleHistory <- scaleHistory
retList$acceptanceHistory <- acceptanceHistory
retList$propCovHistory <- propCovHistory
}
  return(retList)
}






#' Fit a full or reduced state-space model using either MCMC or SMC methods
#'
#' @description Replicate MCMC samples from an already fitted model.
#'
#' @param mcmcout MCMC samples.
#' @param N Number of samples needed now.
#' @param nChains  Number of MCMC chains.
#' @author  Kwaku Peprah Adjei
#'
#'  This function provides an implementation to fit reduced models that can be used to
#'  update the model parameters and latent state distribution in the future i.e. fit the reduced
#'  model. It can also be used to fit the full model that can be compared to the
#'  reduced and updated models
#'
#' @export
#'
#' @family smc update methods
#'
#' @examples
#' ## For illustration only.

replicateMCMCchains <- function(mcmcout, N, nChains){

  for(j in 1:nChains){
    nRows <- nrow(mcmcout$samples[[j]])
    nCols <- ncol(mcmcout$samples[[j]])
    samplingSamples <- sample(1:nRows, N, replace = TRUE)
    newMatrix <- matrix(NA, nrow = N, ncol = nCols)
    colnames(newMatrix) <- colnames(mcmcout$samples[[j]])
    for(i in 1:N){
      newMatrix[i,] <-  mcmcout$samples[[j]][samplingSamples[i],]
    }
    mcmcout$samples[[j]] <- newMatrix
  }
  return(mcmcout)

}

# Convert Sparta results to mcmc results
convertSpartaResults <- function(out, samplerTimes = TRUE){
  if(!is.character(class(out)) && class(out) == "occDet") stop("Class must be of occDet")

  #extract number of chains
  n.chains <- out$BUGSoutput$n.chains
  samplesList  <- vector('list', n.chains)

  samplesList <- lapply(1:n.chains, function(x){
    ret <- out$BUGSoutput$sims.array[,x,]
    #get the results in an appropriate format to match the node names in NIMBLE
    #this works for nodes with dimension greater than one
    colnames(ret) <- stringr::str_replace_all(colnames(ret), c("," = ", "))

if(samplerTimes){
  samplerTimes <- matrix(1000, nrow = nrow(ret), ncol = 4)
  colnames(samplerTimes) <- paste0("samplerTimes[", 1:4, "]")

  # Put both ret and samplerTimes together
  ret <- cbind(ret, samplerTimes)
}

    return(ret)
  }
  )
  names(samplesList)  <- paste0('chain', 1:n.chains)


  # Target samples of interest
  samplesList1 <- coda::as.mcmc.list(lapply(samplesList, function(x){coda::as.mcmc(x)}))

  #Estimating summary
  summaryObject <- lapply(samplesList1, samplesSummary)
  names(summaryObject) <- paste0('chain', 1:n.chains)
  summaryObject$all.chains <- samplesSummary(do.call('rbind', samplesList1))

  #Return results
  retList <- list()
  retList$samples <- samplesList
  retList$summaryObject <- summaryObject
  retList$n.keep <- out$BUGSoutput$n.keep
  return(retList)
}

