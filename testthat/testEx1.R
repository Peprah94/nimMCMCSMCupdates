## load the nimble library and set seed
library('nimble')
library(nimbleSMC)
devtools::load_all(".")
set.seed(1)
mcmcRun <- FALSE #use mcmc or nimbleSMC for reduced Model
pfTypeRun = "bootstrap"
#source("spartaUpdatingFunctions.R")
#source("~/Documents/GitHub/myphdthesis/R/MCMCSamplersUpdate.R")
#load("case2SimData.RData")
#load("reducedCase2.RData")

# Setting up MCMC configuration values and variation of parameters
nIterations = 50000
nBurnin = 20000
nChains = 2
nThin = 5
nyears = 50
aVars <- c(0.1, 0.8) # changing the intercept
#High and small values of a
iNodePrev <- 25  # The number of years for reduced model

aVarstag = 2
#for(aVarstag in 1:length(aVars)){

#Function to simulate data
sim2 <- function(a, b, c, t, mu0){
  x <- y <- numeric(t)
  x[1] <- rnorm(1, mu0, 1 )
  y[1] <- rnorm(1, x[1], 1)

  for(k in 2:t){
    x[k] <- rnorm(1, a*x[k -1] + b, 1)
    y[k] <- rnorm(1, x[k]*c, 1)# + (sigOE * (sqrt(df -2)/df) * rt(1, df))
  }
  return(list(x=x, y=y))
}

message("simulating data for a = ", aVars[aVarstag])
simData <- sim2(a = aVars[aVarstag],
                b = 1,
                c = 1.5,
                t = nyears,
                mu0 = 0.2)

#save data
#save(simData, file = paste0("example1SimData",aVarstag,".RData"))

####################
# define the model
# ###########################

stateSpaceCode <- nimbleCode({
  x[1] ~ dnorm(mu0, 1)
  y[1] ~ dnorm(x[1], 1)
  for(i in 2:t){
    x[i] ~ dnorm(x[i-1] * a + b, 1)
    y[i] ~ dnorm(x[i] * c, 1)
  }
  a ~ dnorm(0, 1)
  b ~ dnorm(0, 1)
  c ~ dnorm(0,1)
  mu0 ~ dnorm(0, 1)
})
#
# ## define data, constants, and initial values
data <- list(
  #   #y = c(0.213, 1.025, 0.314, 0.521, 0.895, 1.74, 0.078, 0.474, 0.656, 0.802)
  y = simData$y
)
constants <- list(
  t = nyears
)
inits <- list(
  a = 0.1,
  b = 0,
  mu0= 0.2,
  c = 1
)
#
#
# ## build the model
stateSpaceModel <- nimbleModel(stateSpaceCode,
                               data = data,
                               constants = constants,
                               inits = inits,
                               check = FALSE)
#

message("Running baseline model (particle MCMC from nimbleSMC) for a = ", aVars[aVarstag])

#Define a new model
newModel <- stateSpaceModel$newModel(replicate = TRUE)

# Function to run the baseline model
# baselineModel <- myphdthesis::baselineSpartaEstimation(model = newModel, latent = "x",
#                                                        MCMCconfiguration = list(target = c('a', 'b', 'c', 'mu0'),
#                                                                                 additionalPars = "x",
#                                                                                 n.iter = nIterations,
#                                                                                 n.chains = nChains,
#                                                                                 n.burnin = nBurnin,
#                                                                                 n.thin = nThin))

#Save results
#save(baselineModel, file = paste0("example1BaselineSMC",aVarstag,".RData"))


############################
# Fitting reduced and updated model
###########################

#for(iNodetag in length(iNodePrev)){
  ####################
  # Reduced model
  #####################

  message(paste("Running reduced model for iNodePrev = ", iNodePrev, "and a = ", aVars[aVarstag]))
  message(paste("Reduced model configuration for iNodePrev = ", iNodePrev, "and a = ", aVars[aVarstag]))
  data <- list(
    #y = c(0.213, 1.025, 0.314, 0.521, 0.895, 1.74, 0.078, 0.474, 0.656, 0.802)
    y = simData$y[-c((iNodePrev+1):50)]
  )
  constants <- list(
    t = iNodePrev
  )
  newModelReduced <- nimbleModel(stateSpaceCode,
                                 data = data,
                                 constants = constants,
                                 inits = inits,
                                 check = FALSE)

  message(paste("Fitting reduced model for iNodePrev = ", iNodePrev, "and a = ", aVars[aVarstag]))
  example1ReducedModel <- spartaNimWeights(model = newModelReduced,
                                           latent = "x",
                                           mcmc = mcmcRun,
                                           pfType = pfTypeRun,
                                           MCMCconfiguration = list(target = c('a', 'b', 'c', 'mu0'),
                                                                    additionalPars = "x",
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
  message(paste("Running updated model for iNodePrev = ", iNodePrev, "and a = ", aVars[aVarstag]))
  message(paste("Running reduced model for iNodePrev = ", iNodePrev, "and a = ", aVars[aVarstag]))
  message(paste("Reduced model configuration for iNodePrev = ", iNodePrev, "and a = ", aVars[aVarstag]))
  data <- list(
    #y = c(0.213, 1.025, 0.314, 0.521, 0.895, 1.74, 0.078, 0.474, 0.656, 0.802)
    y = simData$y[c((iNodePrev):50)]
  )
  constants <- list(
    t = 50 - iNodePrev + 1
  )
  newModelUpdated <- nimbleModel(stateSpaceCode,
                                 data = data,
                                 constants = constants,
                                 inits = inits,
                                 check = FALSE)

  example1UpdatedModel <- spartaNimUpdates(model = newModelUpdated, #nimble model
                                           reducedModel = newModelReduced,
                                           latent = "x", #latent variable
                                           #nParFiltRun = 1000,
                                           pfType = pfTypeRun,
                                           MCMCconfiguration = list(target = c('a', 'b', 'c', 'mu0'),
                                                                    additionalPars = "x",
                                                                    n.iter = (nIterations - nBurnin)/nThin,
                                                                    n.chains = nChains,
                                                                    n.burnin = 0,
                                                                    n.thin = 1),  #saved loglikelihoods from reduced model
                                           postReducedMCMC = example1ReducedModel,# MCMC summary to use as initial values
                                           pfControl = list(saveAll = TRUE,
                                                            lookahead = "mean",
                                                            smoothing = TRUE,
                                                            mcmc = mcmcRun,
                                                            M = nyears - iNodePrev,
                                                            iNodePrev = 1)
  )

  #save(example1UpdatedModel, file= paste0("example1UpdatedIn",iNodePrev[iNodetag],"A",aVarstag,".RData"))

#}

#}
