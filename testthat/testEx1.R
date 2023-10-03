## load the nimble library and set seed
library('nimble')
library(nimbleSMC)
#devtools::load_all(".")
set.seed(1)
mcmcRun <- TRUE #use mcmc or nimbleSMC for reduced Model
pfTypeRun = "bootstrap"
#source("spartaUpdatingFunctions.R")
#source("~/Documents/GitHub/myphdthesis/R/MCMCSamplersUpdate.R")
#load("case2SimData.RData")
#load("reducedCase2.RData")

# Setting up MCMC configuration values and variation of parameters
nIterations = 30000
nBurnin = 20000
nChains = 2
nThin = 2
nyears = 50
aVars <- c(0.1, 0.5) # changing the intercept
#High and small values of a
iNodePrev <- 49  # The number of years for reduced model

aVarstag = 2
#for(aVarstag in 1:length(aVars)){

#Function to simulate data
sim2 <- function(a, c, t, seed){
  set.seed(seed)
  x <- y <- numeric(t)
  x[1] <- rnorm(1, a, sd = 1 )
  y[1] <- rnorm(1, x[1], sd=  1)

  for(k in 2:t){
    x[k] <- rnorm(1, a*x[k -1], sd = 1)
    y[k] <- rnorm(1, x[k]*c, sd= 1)
  }
  return(list(x=x, y=y))
}

message("simulating data for a = ", aVars[aVarstag])
simData <- sim2(a = aVars[aVarstag],
                c = 1,
                t = nyears,
                seed = 1)

#save data
#save(simData, file = paste0("example1SimData",aVarstag,".RData"))

####################
# define the model
# ###########################
# NIMBLE CODE
stateSpaceCode <- nimbleCode({
  x[1] ~ dnorm(a, sd = 1)
  #x[1] <- 0
  y[1] ~ dnorm(c*x[1], sd = 1)
  for(i in 2:t){
    x[i] ~ dnorm(a*x[i-1], sd = 1)
    y[i] ~ dnorm(x[i] * c, 1)
  }
  a ~ dunif(0, 0.999)
  #b ~ dnorm(0, 1)
  c ~ dnorm(1,1)
  #d ~ T(dnorm(0, 1), 0.001, 4)


  #estimate biases
  ahat <- a - aTrue
  #bhat <- b - bTrue
  chat <- c - cTrue
})

# ## define data
data <- list(
  y = simData$y
)

# define constants
constants <- list(
  t = nyears,
  aTrue = 0.5,
  #bTrue = 0.5,
  cTrue = 1,
  dTrue = 1
)

#define initial values
inits <- list(
  a = 0.2,
  #b = 0.8,
  # mu0= 0.8,
  c = 1
)

# ## build the model
stateSpaceModel <- nimbleModel(stateSpaceCode,
                               data = data,
                               constants = constants,
                               inits = inits,
                               check = FALSE)



  data <- list(
    y = simData$y[-c((iNodePrev+1):50)]
  )
  constants <- list(
    t = iNodePrev,
    aTrue = 0.5,
    bTrue = 0.5,
    cTrue = 1,
    dTrue = 1
  )

  #Define model for reduced model
  newModelReduced <- nimbleModel(stateSpaceCode,
                                 data = data,
                                 constants = constants,
                                 inits = inits,
                                 check = FALSE)


  message(paste("Fitting reduced model for iNodePrev = ", iNodePrev, "and a = ", aVars[aVarstag]))
  example1ReducedModel <- nimMCMCSMCupdates::spartaNimWeights(model = newModelReduced,
                                                              latent = "x",
                                                              nParFiltRun = numParticles,
                                                              mcmc = TRUE,
                                                              block = TRUE,
                                                              pfType = pfTypeRun,
                                                              MCMCconfiguration = list(target = c('a',  'c'),
                                                                                       additionalPars = c("x", "ahat", "chat"),
                                                                                       n.iter = nIterations,
                                                                                       n.chains = nChains,
                                                                                       n.burnin = nBurnin,
                                                                                       n.thin = nThin)
  )

  #save results
  #save(example1ReducedModel, file= paste0("example1ReducedIn",iNodePrev[iNodetag],"A",aVarstag,".RData"))

  ################
  # Updated Model
  ################
  data <- list(
    y = simData$y#[c((iNodePrev[iNodetag]):50)]
  )
  #define initial values
  inits <- list(
    a = example1ReducedModel$summary$all.chains["a",1 ],
    #b = 0.8,
    # mu0= 0.8,
    c = example1ReducedModel$summary$all.chains["c",1 ]
  )

  # constants for updated model
  constants <- list(
    t = 50 ,#- iNodePrev[iNodetag] + 1,
    aTrue = 0.5,
    #bTrue = 0.5,
    cTrue = 1
  )


  newModelUpdated <- nimbleModel(stateSpaceCode,
                                 data = data,
                                 constants = constants,
                                 inits = inits,
                                 check = FALSE)

  example1UpdatedModel <- spartaNimUpdates(model = newModelUpdated, #nimble model
                                                              reducedModel = newModelReduced,
                                                              latent = "x", #latent variable
                                                              nParFiltRun = 1000,
                                                              propCov = diag(2)* c(0.1,1),
                                                              pfType = "bootstrap",
                                                              extraVars = NULL,
                                                              MCMCconfiguration = list(target = c('a', 'c'),
                                                                                       additionalPars = c("x", "ahat", "chat"),
                                                                                       n.iter = (nIterations - nBurnin)/nThin,
                                                                                       n.chains = nChains,
                                                                                       n.burnin = 0,
                                                                                       n.thin = 1),  #saved loglikelihoods from reduced model
                                                              postReducedMCMC = example1ReducedModel,# MCMC summary to use as initial values
                                                              pfControl = list(saveAll = TRUE,
                                                                               smoothing = TRUE,
                                                                               mcmc = TRUE,
                                                                               M = nyears - iNodePrev,
                                                                               iNodePrev = iNodePrev)
  )


  #check convergence of MCMC samples

 ggmcmc::ggs(example1UpdatedModel$samples)%>%
   ggmcmc::ggs_traceplot(., family = c("c"))

 ggmcmc::ggs(example1ReducedModel$samples)%>%
   ggmcmc::ggs_traceplot(., family = c("c"))

  #save(example1UpdatedModel, file= paste0("example1UpdatedIn",iNodePrev[iNodetag],"A",aVarstag,".RData"))

#}

#}

 #### Unit testing for particle filters
 model = newModelUpdated #nimble model
 reducedModel = newModelReduced
 latent = "x" #latent variable
 nParFiltRun = 1000
 propCov = diag(2)* c(0.1,1)
 pfType = "bootstrap"
 extraVars = NULL
 MCMCconfiguration = list(target = c('a', 'c'),
                          additionalPars = c("x", "ahat", "chat"),
                          n.iter = (nIterations - nBurnin)/nThin,
                          n.chains = nChains,
                          n.burnin = 0,
                          n.thin = 1)  #saved loglikelihoods from reduced model
 postReducedMCMC = example1ReducedModel# MCMC summary to use as initial values
 pfControl = list(saveAll = TRUE,
                  smoothing = TRUE,
                  mcmc = TRUE,
                  M = nyears - iNodePrev,
                  iNodePrev = iNodePrev)

 target = MCMCconfiguration[["target"]]
 additionalPars = MCMCconfiguration[["additionalPars"]] #other dependent variables you seek to monitor
 n.iter = MCMCconfiguration[["n.iter"]]
 n.chains = MCMCconfiguration[["n.chains"]]
 n.burnin = MCMCconfiguration[["n.burnin"]]
 n.thin = MCMCconfiguration[["n.thin"]]
 iNodePrev = pfControl[["iNodePrev"]]
 M = pfControl[["M"]]
 timeIndex <- pfControl[["timeIndex"]]
 timeIndex = 1
 chain.iter = 1
 updateVars <- updateUtils(model = model, #reduced model
                           reducedModel = reducedModel,
                           mcmcOut = postReducedMCMC$samples[[chain.iter]],
                           latent = latent,
                           target = c(target, additionalPars),
                           n.iter = n.iter ,
                           m = nParFiltRun,
                           iNodeAll = TRUE,
                           timeIndex = timeIndex)

 mvSamplesEst =  updateVars$mvSamplesEst #saved weights
 all(unlist(as.matrix(mvSamplesEst["a"])) == as.matrix(postReducedMCMC$samples[[chain.iter]])[,"a"])
 all(unlist(as.matrix(mvSamplesEst["c"])) == as.matrix(postReducedMCMC$samples[[chain.iter]])[,"c"])

 rr <- lapply(1:1000, function(x){
   all(as.matrix(mvSamplesEst["x"])[[x]][1:10] == as.matrix(postReducedMCMC$samples[[chain.iter]])[x,5:14])
 }
 )%>%
   do.call("c",.)%>%
   sum

 testthat::expect_equal(rr, 1000)

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
 targetAsScalar <- estimationModel$expandNodeNames(target, returnScalarComponents = TRUE)
 compiledParticleFilter <- compileNimble(estimationModel,  particleFilterEst)

 compiledParticleFilter$particleFilterEst$run(m = 1000, iterRun = 10, storeModelValues = values(estimationModel, targetAsScalar))

 rr <- lapply(1:1000, function(x){
   all(compiledParticleFilter$particleFilterEst$mvEWSamples[["x"]][[x]][1:10] == as.matrix(postReducedMCMC$samples[[chain.iter]])[10,5:14])
 }
 )%>%
   do.call("c",.)%>%
   sum()

 testthat::expect_equal(rr, 1000)


 message("Setting up the MCMC Configuration")
 #newModel <- model$newModel(replicate = TRUE)
 modelMCMCconf <- nimble::configureMCMC(model,
                                        monitors = c(target, latent, additionalPars),
                                        nodes = NULL)#, monitors = c(target, latent, additionalPars))

 # set sampler to use
 if(pfType == "bootstrap"){
   pfTypeUpdate = 'bootstrapUpdate'
 }else{
   if(pfType == "auxiliary") {
     pfTypeUpdate = 'auxiliaryUpdate'
   }
 }

 propCov1 <- "identity"
 mcmcScale <- 1
 #set up MCMC configuration


 #######
 #Unit testing for MCMC Sampler
 control = list(latents = latent,
                #target = target,
                adaptive = FALSE,
                propCov = propCov,
                propCov1 = propCov1,
                #adaptInterval = 100,
                scale = mcmcScale,
                # pf = particleFilter,
                pfControl = pfControl, #list( M = M, iNodePrev = iNodePrev),
                pfNparticles = nParFiltRun,
                pfType = pfTypeUpdate,
                postSamples = postReducedMCMC$samples[[chain.iter]],
                mvSamplesEst = mvSamplesEst,
                extraVars = extraVars,
                iNodePrev = iNodePrev,
                additionalPars = additionalPars,
                reducedModel = reducedModel)








 modelMCMCconf$addSampler(target = target,
                          type = 'RW_PF_blockUpdate',
                          control = list(latents = latent,
                                         #target = target,
                                         adaptive = FALSE,
                                         propCov = propCov,
                                         propCov1 = propCov1,
                                         #adaptInterval = 100,
                                         scale = mcmcScale,
                                         # pf = particleFilter,
                                         pfControl = pfControl, #list( M = M, iNodePrev = iNodePrev),
                                         pfNparticles = nParFiltRun,
                                         pfType = pfTypeUpdate,
                                         postSamples = postReducedMCMC$samples[[chain.iter]],
                                         mvSamplesEst = mvSamplesEst,
                                         extraVars = extraVars,
                                         iNodePrev = iNodePrev,
                                         additionalPars = additionalPars,
                                         reducedModel = reducedModel)
 )

 modelMCMCconf$addMonitors(additionalPars)

 message("Building and compiling the PF MCMC")
 ## build and compile pMCMC sampler
 modelMCMC <- buildMCMC(modelMCMCconf)
 compiledList <- nimble::compileNimble(model,
                                       modelMCMC)

 mcmc.out <- nimble::runMCMC(compiledList$modelMCMC,
                             niter = 150,
                             nchains = 1,
                             nburnin = n.burnin,
                             #inits = initsList,
                             thin = n.thin,
                             setSeed = TRUE,
                             samples= 123,
                             samplesAsCodaMCMC = TRUE,
                             summary = TRUE,
                             WAIC = FALSE)


 rr <- lapply(1:1000, function(x){
   all(as.matrix(mcmc.out$samples) == as.matrix(postReducedMCMC$samples[[chain.iter]])[10,5:14])
 }
 )%>%
   do.call("c",.)%>%
   sum()

 testthat::expect_equal(rr, 1000)
