

# Occupancy models
library(nimble)
library(nimbleSMC)
library(nimMCMCSMCupdates)

# set the configurations for fitting the model
nyears = 10
nsites = 30
nvisits = 6
iNodePrev = 7
pfTypeRun = "auxiliary"

# Set the configurations for MCMC
nIterations = 1000
nBurnin = 50
nChains = 2
nThin = 2
#nThin = 1
numParticles = 10


# Function to simulate dynamic occupancy models
dynOccupancyModels <- function(nyears,
                               nsites,
                               nvisits,
                               fixedPars = list(alphaPhi = 2,
                                                betaPhi = 1.5),
                               hyperParsSig = list(alphaP = 2,
                                                   betaP = 3,
                                                   alphaGamma = 1,
                                                   betaGamma = 1.5,
                                                   alphaPsi = 2,
                                                   betaPsi = 2)){

  # Simulating covariates
  set.seed(1994)

  windSpeed <- array(runif(nyears*nsites*nvisits, -1, 1),
                     dim = c(nsites, nvisits, nyears))
  elevation <- runif(nsites, -1,1)
  springPrecipitation <- matrix(runif(nyears*nsites, -1, 1),
                                nrow = nsites,
                                ncol = nyears)
  sizeOfBeak <- matrix(runif(nyears*nsites, -1,1),
                       nrow = nsites,
                       ncol = nyears)

  # Simulating parameters
  alphaP <- rnorm(1, mean = 0, sd = hyperParsSig$alphaP)
  betaP <- rnorm(1, mean = 0, sd = hyperParsSig$betaP)
  alphaPsi <- rnorm(1, mean = 0, sd = hyperParsSig$alphaPsi)
  betaPsi <- rnorm(1, mean = 0, sd = hyperParsSig$betaPsi)
  alphaGamma <- rnorm(1, mean = 0, sd = hyperParsSig$alphaGamma)
  betaGamma <- rnorm(1, mean = 0, sd = hyperParsSig$betaGamma)

  # Detection Probability
  detectProb <- array(NA, dim = c(nsites, nvisits, nyears))

  for(year.tag in 1:nyears){
    for(site.tag in 1:nsites){
      for(visit.tag in 1:nvisits){
        detectProb[site.tag, visit.tag, year.tag] <- plogis(alphaP + betaP* windSpeed[site.tag, visit.tag, year.tag])
      }
    }
  }

  # Initial occupancy probability psi
  initOccuProb <- plogis(fixedPars$alphaPhi + fixedPars$betaPhi*elevation)

  # Persistence and colonisation probability
  persistenceProb  <- matrix(NA,  nrow = nsites, ncol = nyears)
  for(year.tag in 1:nyears){
    for(site.tag in 1:nsites){
      persistenceProb[site.tag,  year.tag] <- plogis(alphaPsi + betaPsi* springPrecipitation[site.tag, year.tag])
      # colonisationProb[site.tag,  year.tag] <- plogis(alphaGamma[year.tag] + betaGamma[year.tag]* sizeOfBeak[site.tag, year.tag])
    }
  }

  colonisationProb <- 0.05

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
      z[ site.tag, year.tag] <- rbinom(1, 1, z[site.tag, (year.tag -1)] * persistenceProb[site.tag,  year.tag] + (1 - z[site.tag, (year.tag -1)])*colonisationProb)#[site.tag,  year.tag])
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

# Simulating and saving data
simData <- dynOccupancyModels(nyears = 10,
                              nsites = 30,
                              nvisits = 6)

# save(simData, file = "simDataDynamicOccupancy.RData")

load("simDataDynamicOccupancy.RData")

# NIMBLE code
dynOccupancyCode <- nimbleCode({

  # Prior distributions of hyperparameters
  alphaPSig ~ T(dnorm(0, 0.1), 0.001, )
  betaPSig ~ T(dnorm(0, 0.1), 0.001, )
  alphaPsiSig ~ T(dnorm(0, 0.1), 0.001, )
  betaPsiSig ~ T(dnorm(0, 0.1), 0.001, )


  alphaPhi ~ dnorm(0, sd = 10)
  betaPhi ~ dnorm(0, sd = 10)

  # Prior distributions
  alphaP ~ dnorm(mean = 0, sd = alphaPSig)
  betaP ~ dnorm(mean = 0, sd =  betaPSig)
  alphaPsi ~ dnorm(mean = 0, sd =  alphaPsiSig)
  betaPsi ~ dnorm(mean = 0, sd =  betaPsiSig)

  # Detection Probability
  for(year.tag in 1:nyears){
    for(site.tag in 1:nsites){
      for(visit.tag in 1:nvisits){
        logit(detectProb[site.tag, visit.tag, year.tag]) <- alphaP + betaP* windSpeed[site.tag, visit.tag, year.tag]
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
      logit(persistenceProb[site.tag,  year.tag]) <- alphaPsi + betaPsi* springPrecipitation[site.tag, year.tag]
      # logit(colonisationProb[site.tag,  year.tag]) <- alphaGamma[year.tag] + betaGamma[year.tag]* sizeOfBeak[site.tag, year.tag]
    }
  }

  colonisationProb ~ dunif(0.01, 0.5)


  # Initial Presence/Absence
  for(site.tag in 1:nsites){
    z[site.tag, 1] ~ dbin(prob = initOccuProb[site.tag], 1)
  }

  # True presence/absence
  for(year.tag in 2:nyears){
    for(site.tag in 1:nsites){
      z[site.tag, year.tag] ~ dbin(prob = (z[site.tag, (year.tag -1)] * persistenceProb[site.tag,  year.tag] + (1 - z[site.tag, (year.tag -1)])*colonisationProb), 1)
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


# ## define data, constants, and initial values
z <- apply(simData$y, c(1, 3), function(x){
  ifelse(sum(x) > 0, 1, NA)
})


#####################
#   Reduced Model
########################

data <- list(
  y = simData$y[,,-c((iNodePrev +1):nyears)],
  windSpeed = simData$covariates$windSpeed[,,-c((iNodePrev +1):nyears)],
  elevation = simData$covariates$elevation,
  springPrecipitation = simData$covariates$springPrecipitation[,-c((iNodePrev +1):nyears)],
  sizeOfBeak = simData$covariates$sizeOfBeak,
  z = z[,-c((iNodePrev +1):nyears)]
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
  alphaP =  rnorm(1, mean = 0, sd = 1),
  betaP =  rnorm(1, mean = 0, sd = 1),
  alphaPsi=  rnorm(1, mean = 0, sd = 1),
  betaPsi= rnorm(1, mean = 0, sd = 1),
  colonisationProb = 0.01,
  alphaPhi=  rnorm(1),
  betaPhi= rnorm(1)
)

# NIMBLE model for reduced model
newModelReduced <- nimbleModel(dynOccupancyCode,
                               data = data,
                               constants = constants,
                               inits = inits,
                               check = FALSE)

# Fitting the model
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
                                                                                 'alphaPhi', 'betaPhi',
                                                                                 'alphaP', 'betaP',
                                                                                 'alphaPsi', 'betaPsi',
                                                                                 "colonisationProb"
                                             ),
                                             additionalPars = c( "psi.fs"
                                             ),
                                             n.iter = nIterations,
                                             n.chains = nChains,
                                             n.burnin = nBurnin,
                                             n.thin = nThin)
)
#save results
#save(example2ReducedModelTrue, file= "example4ReducedBootstrapTrue.RData")
load("example4ReducedBootstrapTrue.RData")
################
# Updated Model
################

data <- list(
  y = simData$y,#[,,c((iNodePrev ):nyears)],
  windSpeed = simData$covariates$windSpeed,#[,,c((iNodePrev ):nyears)],
  elevation = simData$covariates$elevation,
  springPrecipitation = simData$covariates$springPrecipitation,#[,c((iNodePrev ):nyears)],
  sizeOfBeak = simData$covariates$sizeOfBeak,
  z = z#[,c((iNodePrev ):nyears)]
)

constants <- list(
  nyears = nyears,# - iNodePrev +1,
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
  alphaP =  rnorm(1, mean = 0, sd = 1),
  betaP =  rnorm(1, mean = 0, sd = 1),
  alphaPsi=  rnorm(1, mean = 0, sd = 1),
  betaPsi= rnorm(1, mean = 0, sd = 1),
  colonisationProb = 0.01,
  # alphaGamma=  rnorm(1, mean = 0, sd = 1),
  # betaGamma=  rnorm(1, mean = 0, sd = 1),
  alphaPhi=  rnorm(1),
  betaPhi= rnorm(1)
)

# NIMBLE model for updated model
newModelUpdated <- nimbleModel(dynOccupancyCode,
                               data = data,
                               constants = constants,
                               inits = inits,
                               check = FALSE)

# Fitting the model
example2UpdatedModelTrue <- spartaNimUpdates(model = newModelUpdated, #nimble model
                                             reducedModel = newModelReduced,
                                             latent = "z", #latent variable
                                             nParFiltRun = numParticles,
                                             #nParFiltRun = 1000,
                                             pfType = pfTypeRun,
                                             propCov = c(1,1,1,1,1,1,0.01)*diag(7),
                                             mcmcScale = 1,
                                             extraVars = NULL,
                                             MCMCconfiguration = list(target = c('alphaPSig', 'betaPSig',
                                                                                 'alphaPsiSig', 'betaPsiSig',
                                                                                 'alphaPhi', 'betaPhi',
                                                                                 'alphaP', 'betaP',
                                                                                 'alphaPsi', 'betaPsi',
                                                                                 "colonisationProb"),
                                                                      additionalPars = c("z", "psi.fs"),
                                                                      n.iter = (nIterations - nBurnin)/nThin,
                                                                      n.chains = nChains,
                                                                      n.burnin = 10,
                                                                      n.thin = 1
                                             ),  #saved loglikelihoods from reduced model
                                             postReducedMCMC = example2ReducedModelTrue,# MCMC summary to use as initial values
                                             pfControl = list(saveAll = TRUE,
                                                              timeIndex = 2,
                                                              smoothing = TRUE,
                                                              mcmc = TRUE,
                                                              M = nyears - iNodePrev,
                                                              iNodePrev = iNodePrev)
)

#save(example2UpdatedModelTrue, file= "example4UpdatedBootstrapTrue1.RData")



save(example2UpdatedModelTrue, file= "example2UpdatedAuxiliaryTrue1.RData")


#############
# Baseline Model
###############

# data <- list(
#   y = simData$y,
#   windSpeed = simData$covariates$windSpeed,
#   elevation = simData$covariates$elevation,
#   springPrecipitation = simData$covariates$springPrecipitation,
#   sizeOfBeak = simData$covariates$sizeOfBeak,
#   z = z
# )
#
# constants <- list(
#   nyears = nyears,
#   nsites = nsites,
#   nvisits = nvisits
# )
# inits <- list(
#   alphaPSig = 1,
#   betaPSig  = 1,
#   alphaPsiSig = 1,
#   betaPsiSig  = 1,
#   alphaGammaSig = 1,
#   betaGammaSig = 1,
#   alphaPhi  = 0,
#   betaPhi  = 0,
#   alphaP =  rnorm(1, mean = 0, sd = 1),
#   betaP =  rnorm(1, mean = 0, sd = 1),
#   alphaPsi=  rnorm(1, mean = 0, sd = 1),
#   betaPsi= rnorm(1, mean = 0, sd = 1),
#   #alphaGamma=  rnorm(constants$nyears, mean = 0, sd = 1),
#   #betaGamma=  rnorm(constants$nyears, mean = 0, sd = 1),
#   alphaPhi=  rnorm(1),
#   betaPhi= rnorm(1)
# )
# #
# #
# # ## build the model
# dynOccModel <- nimbleModel(dynOccupancyCode,
#                            data = data,
#                            constants = constants,
#                            inits = inits,
#                            check = FALSE)
#
#
# # Fitting Baseline model
# newModel <- dynOccModel$newModel(replicate = TRUE)
#
# # Function to run the baseline model
# baselineModel <- nimMCMCSMCupdates::baselineSpartaEstimation(model = newModel,
#                                                              latent = "z",
#                                                              nParFiltRun = numParticles,
#                                                              pfControl = list(saveAll = TRUE,
#                                                                               smoothing = FALSE,
#                                                                               timeIndex = 2),
#                                                              MCMCconfiguration = list(target = c('alphaPSig', 'betaPSig',
#                                                                                                  'alphaPsiSig', 'betaPsiSig',
#                                                                                                  'alphaPhi', 'betaPhi',
#                                                                                                  'alphaP', 'betaP',
#                                                                                                  'alphaPsi', 'betaPsi',
#                                                                                                  "colonisationProb"),
#                                                                                       additionalPars = c("z", "psi.fs",
#                                                                                                          'alphaP', 'betaP',
#                                                                                                          'alphaPsi', 'betaPsi'),
#                                                                                       n.iter = nIterations,
#                                                                                       n.chains = nChains,
#                                                                                       n.burnin = nBurnin,
#                                                                                       n.thin = nThin))
# #save results
# save(baselineModel, file = "example2BaselineSMCAuxiliary.RData")
#
#
#
