# Occupancy models
library('nimble')
library(nimbleSMC)
#library(nimMCMCSMCupdates)
#devtools::install_github("Peprah94/nimMCMCSMCupdates")
devtools::load_all(".")
nimbleOptions(MCMCusePredictiveDependenciesInCalculations = TRUE)
nimbleOptions(MCMCorderPosteriorPredictiveSamplersLast = FALSE)
#library(myphdthesis)
nyears = 10
nsites = 5
nvisits = 3
iNodePrev = 7
pfTypeRun = "auxiliary"

nIterations = 5000
nBurnin = 2000
nChains = 2
nThin = 1
numParticles = 10

dynOccupancyModels <- function(nyears,
                               nsites,
                               nvisits,
                               fixedPars = list(alphaPhi = 2,
                                                betaPhi = 1.5),
                               hyperParsSig = list(alphaP = 2,
                                                   betaP = 0.5,
                                                   alphaGamma = 1,
                                                   betaGamma = 1.5,
                                                   alphaPsi = 0.5,
                                                   betaPsi = 1)){

  # Simulating covariates
  windSpeed <- array(rnorm(nyears*nsites*nvisits),
                     dim = c(nsites, nvisits, nyears))
  elevation <- rnorm(nsites)
  springPrecipitation <- matrix(rnorm(nyears*nsites),
                                nrow = nsites,
                                ncol = nyears)
  sizeOfBeak <- matrix(rnorm(nyears*nsites),
                       nrow = nsites,
                       ncol = nyears)

  # Simulating parameters
  alphaP <- rnorm(nyears, mean = 0, sd = hyperParsSig$alphaP)
  betaP <- rnorm(nyears, mean = 0, sd = hyperParsSig$betaP)
  alphaPsi <- rnorm(nyears, mean = 0, sd = hyperParsSig$alphaPsi)
  betaPsi <- rnorm(nyears, mean = 0, sd = hyperParsSig$betaPsi)
  alphaGamma <- rnorm(nyears, mean = 0, sd = hyperParsSig$alphaGamma)
  betaGamma <- rnorm(nyears, mean = 0, sd = hyperParsSig$betaGamma)

  # Detection Probability
  detectProb <- array(NA, dim = c(nsites, nvisits, nyears))

  for(year.tag in 1:nyears){
    for(site.tag in 1:nsites){
      for(visit.tag in 1:nvisits){
        detectProb[site.tag, visit.tag, year.tag] <- plogis(alphaP[year.tag] + betaP[year.tag]* windSpeed[site.tag, visit.tag, year.tag])
      }
    }
  }

  # Initial occupancy probability psi
  initOccuProb <- plogis(fixedPars$alphaPhi + fixedPars$betaPhi*elevation)

  # Persistence and colonisation probability
  persistenceProb <- colonisationProb <- matrix(NA,  nrow = nsites, ncol = nyears)
  for(year.tag in 1:nyears){
    for(site.tag in 1:nsites){
      persistenceProb[site.tag,  year.tag] <- plogis(alphaPsi[year.tag] + betaPsi[year.tag]* springPrecipitation[site.tag, year.tag])
      colonisationProb[site.tag,  year.tag] <- plogis(alphaGamma[year.tag] + betaGamma[year.tag]* sizeOfBeak[site.tag, year.tag])
    }
  }

  # Latent state and observations
  y <- array(NA, dim = c(nsites, nvisits, nyears))
  z <- matrix(NA,nrow = nsites,ncol = nyears)

  # Initial Presence/Absence
  for(site.tag in 1:nsites){
    z[ site.tag, 1] <- rbinom(1, 1, initOccuProb[site.tag])
  }

  # True presence/absence
  for(year.tag in 2:nyears){
    for(site.tag in 1:nsites){
      z[ site.tag, year.tag] <- rbinom(1, 1, z[site.tag, (year.tag -1)] * persistenceProb[site.tag,  year.tag] + (1 - z[site.tag, (year.tag -1)])*colonisationProb[site.tag,  year.tag])
    }
  }


  #observations
  for(year.tag in 1:nyears){
    for(site.tag in 1:nsites){
      for(visit.tag in 1:nvisits){
        y[site.tag, visit.tag, year.tag] <- rbinom(1, 1,  z[site.tag, year.tag] * detectProb[site.tag, visit.tag, year.tag])
      }
    }
  }

  # Proportion of occupied sites

  psi.fs <- colMeans(z)


  # Return list
  retList <- list()
  retList$y = y
  retList$z = z
  retList$covariates = list(windSpeed = windSpeed,
                            elevation = elevation,
                            springPrecipitation = springPrecipitation,
                            sizeOfBeak = sizeOfBeak)
  retList$trueSigma = hyperParsSig
  retList$truePars = fixedPars
  retList$covariateEffects = list(alphaP = alphaP,
                                  betaP = betaP,
                                  alphaPsi = alphaPsi,
                                  betaPsi = betaPsi,
                                  alphaGamma = alphaGamma,
                                  betaGamma = betaGamma)

  retList$occSites = psi.fs

  return(retList)
}


# Simulating data

simData <- dynOccupancyModels(nyears = 10,
                               nsites = 5,
                              nvisits = 3)
 #save(simData, file = "simDataDynamicOccupancy.RData")

#load("simDataDynamicOccupancy.RData")
# NIMBLE MODEL

dynOccupancyCode <- nimbleCode({

  # Prior distributions of hyperparameters
  alphaPSig ~ dgamma(1,1)
  betaPSig ~ dgamma(1,1)
  alphaPsiSig ~ dgamma(1,1)
  betaPsiSig ~ dgamma(1,1)
  alphaGammaSig ~ dgamma(1,1)
  betaGammaSig ~ dgamma(1,1)

  alphaPhi ~ dnorm(0, sd = 10)
  betaPhi ~ dnorm(0, sd = 10)

  # Prior distributions
  for(year.tag in 1: nyears){
    alphaP[year.tag] ~ dnorm(mean = 0, tau = alphaPSig)
    betaP[year.tag] ~ dnorm(mean = 0, tau =  betaPSig)
    alphaPsi[year.tag] ~ dnorm(mean = 0, tau =  alphaPsiSig)
    betaPsi[year.tag] ~ dnorm(mean = 0, tau =  betaPsiSig)
    alphaGamma[year.tag] ~ dnorm(mean = 0, tau =  alphaGammaSig)
    betaGamma[year.tag] ~ dnorm( mean = 0, tau =  betaGammaSig)
  }

  # Detection Probability
  for(year.tag in 1:nyears){
    for(site.tag in 1:nsites){
      for(visit.tag in 1:nvisits){
        logit(detectProb[site.tag, visit.tag, year.tag]) <- alphaP[year.tag] + betaP[year.tag]* windSpeed[site.tag, visit.tag, year.tag]
      }
    }
  }

  # Initial occupancy probability psi
  for(site.tag in 1:nsites){
    logit(initOccuProb[site.tag]) <- alphaPhi + betaPhi*elevation[site.tag]
  }

  # Persistence and colonisation probability
  for(year.tag in 1:nyears){
    for(site.tag in 1:nsites){
      logit(persistenceProb[site.tag,  year.tag]) <- alphaPsi[year.tag] + betaPsi[year.tag]* springPrecipitation[site.tag, year.tag]
      # logit(colonisationProb[site.tag,  year.tag]) <- alphaGamma[year.tag] + betaGamma[year.tag]* sizeOfBeak[site.tag, year.tag]
    }
  }

  for(year.tag in 1:nyears){
    for(site.tag in 1:nsites){
      #   logit(persistenceProb[site.tag,  year.tag]) <- alphaPsi[year.tag] + betaPsi[year.tag]* springPrecipitation[site.tag, year.tag]
      logit(colonisationProb[site.tag,  year.tag]) <- alphaGamma[year.tag] + betaGamma[year.tag]* sizeOfBeak[site.tag, year.tag]
    }
  }

  # Initial Presence/Absence
  for(site.tag in 1:nsites){
    z[site.tag, 1] ~ dbin(prob = initOccuProb[site.tag], 1)
  }

  # True presence/absence
  for(year.tag in 2:nyears){
    for(site.tag in 1:nsites){
      z[site.tag, year.tag] ~ dbin(prob = (z[site.tag, (year.tag -1)] * persistenceProb[site.tag,  year.tag] + (1 - z[site.tag, (year.tag -1)])*colonisationProb[site.tag,  year.tag]), 1)
    }
  }

  #observations
  for(year.tag in 1:nyears){
    for(site.tag in 1:nsites){
      for(visit.tag in 1:nvisits){
        y[site.tag, visit.tag, year.tag] ~ dbin(prob = z[site.tag, year.tag] * detectProb[site.tag, visit.tag, year.tag], 1)
      }
    }
  }

  # Derived quantities
  for(year.tag in 1:nyears){
    psi.fs[year.tag] <- sum(z[1:nsites, year.tag])/nsites
  }

})
#
# ## define data, constants, and initial values
data <- list(
  y = simData$y,
  windSpeed = simData$covariates$windSpeed,
  elevation = simData$covariates$elevation,
  springPrecipitation = simData$covariates$springPrecipitation,
  sizeOfBeak = simData$covariates$sizeOfBeak
)

constants <- list(
  nyears = nyears,
  nsites = nsites,
  nvisits = nvisits
)
inits <- list(
  alphaPSig = 1,
  betaPSig  = 1,
  alphaPsiSig = 1,
  betaPsiSig  = 1,
  alphaGammaSig = 1,
  betaGammaSig = 1,
  alphaPhi  = 0,
  betaPhi  = 0,
  alphaP =  rnorm(constants$nyears, mean = 0, sd = 1),
  betaP =  rnorm(constants$nyears, mean = 0, sd = 1),
  alphaPsi=  rnorm(constants$nyears, mean = 0, sd = 1),
  betaPsi= rnorm(constants$nyears, mean = 0, sd = 1),
  alphaGamma=  rnorm(constants$nyears, mean = 0, sd = 1),
  betaGamma=  rnorm(constants$nyears, mean = 0, sd = 1),
  alphaPhi=  rnorm(1),
  betaPhi= rnorm(1),
  z = matrix(1, nrow = constants$nsites, ncol = constants$nyears)
)
#
#
# ## build the model
dynOccModel <- nimbleModel(dynOccupancyCode,
                           data = data,
                           constants = constants,
                           inits = inits,
                           check = FALSE)


# Fitting Baseline model
newModel <- dynOccModel$newModel(replicate = TRUE)

# Function to run the baseline model
# baselineModel <-nimMCMCSMCupdates::baselineSpartaEstimation(model = newModel,
#                                                             latent = "z",
#                                                             nParFiltRun = numParticles,
#                                                             pfControl = list(saveAll = TRUE,
#                                                                              smoothing = FALSE,
#                                                                              timeIndex = 2),
#                                                             MCMCconfiguration = list(target = c('alphaPSig', 'betaPSig',
#                                                                                                 'alphaPsiSig', 'betaPsiSig',
#                                                                                                 'alphaGammaSig', 'betaGammaSig',
#                                                                                                 'alphaPhi', 'betaPhi'),
#                                                                                      additionalPars = c("z", "psi.fs",
#                                                                                                         'alphaP', 'betaP',
#                                                                                                         'alphaPsi', 'betaPsi',
#                                                                                                         'alphaGamma', 'betaGamma'),
#                                                                                      n.iter = nIterations,
#                                                                                      n.chains = nChains,
#                                                                                      n.burnin = nBurnin,
#                                                                                      n.thin = nThin))
# #save results
# save(baselineModel, file = "example2BaselineSMCboostrap.RData")

#####################
#   Reduced Model
########################

data <- list(
  y = simData$y[,,-c((iNodePrev +1):nyears)],
  windSpeed = simData$covariates$windSpeed,
  elevation = simData$covariates$elevation,
  springPrecipitation = simData$covariates$springPrecipitation,
  sizeOfBeak = simData$covariates$sizeOfBeak
)

constants <- list(
  nyears = iNodePrev,
  nsites = nsites,
  nvisits = nvisits
)
inits <- list(
  alphaPSig = 1,
  betaPSig  = 1,
  alphaPsiSig = 1,
  betaPsiSig  = 1,
  alphaGammaSig = 1,
  betaGammaSig = 1,
  alphaPhi  = 0,
  betaPhi  = 0,
  alphaP =  rnorm(constants$nyears, mean = 0, sd = 1),
  betaP =  rnorm(constants$nyears, mean = 0, sd = 1),
  alphaPsi=  rnorm(constants$nyears, mean = 0, sd = 1),
  betaPsi= rnorm(constants$nyears, mean = 0, sd = 1),
  alphaGamma=  rnorm(constants$nyears, mean = 0, sd = 1),
  betaGamma=  rnorm(constants$nyears, mean = 0, sd = 1),
  alphaPhi=  rnorm(1),
  betaPhi= rnorm(1),
  z = matrix(1, nrow = constants$nsites, ncol = constants$nyears)
)

# nimbleModel reduced
newModelReduced <- nimbleModel(dynOccupancyCode,
                               data = data,
                               constants = constants,
                               inits = inits,
                               check = FALSE)

example2ReducedModelTrue <- spartaNimWeights(model = newModelReduced,
                                             latent = "z",
                                             nParFiltRun = numParticles,
                                             mcmc = TRUE,
                                             pfType = pfTypeRun,
                                             #pfType = pfTypeRun,
                                             pfControl = list(saveAll = TRUE,
                                                              smoothing = TRUE,
                                                              timeIndex = 2),
                                             MCMCconfiguration = list(target = c('alphaPSig', 'betaPSig',
                                                                                 'alphaPsiSig', 'betaPsiSig',
                                                                                 'alphaGammaSig', 'betaGammaSig',
                                                                                 'alphaPhi', 'betaPhi',
                                                                                 'alphaP', 'betaP',
                                                                                 'alphaPsi', 'betaPsi',
                                                                                 'alphaGamma', 'betaGamma'
                                                                                ),
                                                                      additionalPars = c("z", "psi.fs",
                                                                                         'alphaP', 'betaP',
                                                                                         'alphaPsi', 'betaPsi',
                                                                                         'alphaGamma', 'betaGamma'
                                                                                         ),
                                                                      n.iter = nIterations,
                                                                      n.chains = nChains,
                                                                      n.burnin = nBurnin,
                                                                      n.thin = nThin)
)

#save results
#save(example2ReducedModelTrue, file= "example2ReducedAuxiliaryTrue.RData")


# example2ReducedModelFalse <- spartaNimWeights(model = newModelReduced,
#                                               latent = "z",
#                                               mcmc = FALSE,
#                                               nParFiltRun = numParticles,
#                                               pfType = pfTypeRun,
#                                               pfControl = list(saveAll = TRUE,
#                                                                smoothing = FALSE,
#                                                                timeIndex = 2),
#                                               MCMCconfiguration = list(target = c('alphaPSig', 'betaPSig',
#                                                                                   'alphaPsiSig', 'betaPsiSig',
#                                                                                   'alphaGammaSig', 'betaGammaSig',
#                                                                                   'alphaPhi', 'betaPhi',
#                                                                                   'alphaPSig', 'betaPSig'),
#                                                                        additionalPars = c("z", "psi.fs",
#                                                                                           'alphaP', 'betaP',
#                                                                                           'alphaPsi', 'betaPsi',
#                                                                                           'alphaGamma', 'betaGamma'),
#                                                                        n.iter = nIterations,
#                                                                        n.chains = nChains,
#                                                                        n.burnin = nBurnin,
#                                                                        n.thin = nThin)
# )
#
# #save results
# save(example2ReducedModelFalse, file= "example2ReducedAuxiliaryFalse.RData")

################
# Updated Model
################

data <- list(
  y = simData$y[,,c((iNodePrev ):nyears)],
  windSpeed = simData$covariates$windSpeed,
  elevation = simData$covariates$elevation,
  springPrecipitation = simData$covariates$springPrecipitation,
  sizeOfBeak = simData$covariates$sizeOfBeak
)

constants <- list(
  nyears = nyears - iNodePrev +1,
  nsites = nsites,
  nvisits = nvisits
)
inits <- list(
  alphaPSig = 1,
  betaPSig  = 1,
  alphaPsiSig = 1,
  betaPsiSig  = 1,
  alphaGammaSig = 1,
  betaGammaSig = 1,
  alphaPhi  = 0,
  betaPhi  = 0,
  alphaP =  rnorm(constants$nyears, mean = 0, sd = 1),
  betaP =  rnorm(constants$nyears, mean = 0, sd = 1),
  alphaPsi=  rnorm(constants$nyears, mean = 0, sd = 1),
  betaPsi= rnorm(constants$nyears, mean = 0, sd = 1),
  alphaGamma=  rnorm(constants$nyears, mean = 0, sd = 1),
  betaGamma=  rnorm(constants$nyears, mean = 0, sd = 1),
  alphaPhi=  rnorm(1),
  betaPhi= rnorm(1),
  z = matrix(1, nrow = constants$nsites, ncol = constants$nyears)
)

# nimbleModel reduced
newModelUpdated <- nimbleModel(dynOccupancyCode,
                               data = data,
                               constants = constants,
                               inits = inits,
                               check = FALSE)
example2UpdatedModelTrue <- spartaNimUpdates(model = newModelUpdated, #nimble model
                                             reducedModel = newModelReduced,
                                             latent = "z", #latent variable
                                             nParFiltRun = numParticles,
                                             #nParFiltRun = 1000,
                                             pfType = pfTypeRun,
                                             MCMCconfiguration = list(target = c('alphaPSig', 'betaPSig',
                                                                                 'alphaPsiSig', 'betaPsiSig',
                                                                                 'alphaGammaSig', 'betaGammaSig',
                                                                                 'alphaPhi', 'betaPhi',
                                                                                 'alphaP', 'betaP',
                                                                                 'alphaPsi', 'betaPsi',
                                                                                 'alphaGamma', 'betaGamma'),
                                                                      additionalPars = c("z", "psi.fs"),
                                                                      n.iter = (nIterations - nBurnin),
                                                                      n.chains = nChains,
                                                                      #n.burnin = 0,
                                                                      n.thin = 1
                                                                      ),  #saved loglikelihoods from reduced model
                                             postReducedMCMC = example2ReducedModelTrue,# MCMC summary to use as initial values
                                             pfControl = list(saveAll = TRUE,
                                                              timeIndex = 2,
                                                              smoothing = TRUE,
                                                              mcmc = TRUE,
                                                              M = nyears - iNodePrev,
                                                              iNodePrev = 1)
)

#save(example2UpdatedModelTrue, file= "example2UpdatedAuxiliaryTrue.RData")


# example2UpdatedModelFalse <- spartaNimUpdates(model = dynOccModel, #nimble model
#                                               reducedModel = newModelReduced,
#                                               latent = "z", #latent variable
#                                               nParFiltRun = numParticles,
#                                               #nParFiltRun = 1000,
#                                               pfType = pfTypeRun,
#                                               MCMCconfiguration = list(target = c('alphaPSig', 'betaPSig',
#                                                                                   'alphaPsiSig', 'betaPsiSig',
#                                                                                   'alphaGammaSig', 'betaGammaSig',
#                                                                                   'alphaPhi', 'betaPhi'),
#                                                                        additionalPars = c("z", "psi.fs",
#                                                                                           'alphaP', 'betaP',
#                                                                                           'alphaPsi', 'betaPsi',
#                                                                                           'alphaGamma', 'betaGamma'),
#                                                                        n.iter = (nIterations - nBurnin)/nThin,
#                                                                        n.chains = nChains,
#                                                                        n.burnin = 0,
#                                                                        n.thin = 1),  #saved loglikelihoods from reduced model
#                                               postReducedMCMC = example2ReducedModelFalse,# MCMC summary to use as initial values
#                                               pfControl = list(saveAll = TRUE,
#                                                                timeIndex = 2,
#                                                                smoothing = FALSE,
#                                                                mcmc = FALSE,
#                                                                M = nyears - iNodePrev,
#                                                                iNodePrev = iNodePrev)
# )
#
#
# save(example2UpdatedModelFalse, file= "example2UpdatedAuxiliaryFalse.RData")
