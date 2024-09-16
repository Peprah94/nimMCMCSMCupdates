# Occupancy models
library(nimble)
library(nimbleSMC)
library(nimMCMCSMCupdates)
nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
nimbleOptions(MCMCsaveHistory = TRUE)


# set the configurations for fitting the model
nyears = 20
nsites = 100
nvisits = 6
iNodePrevAll = c(18, 15)
pfTypeRun = "bootstrap"

# Set the configurations for MCMC
nIterations = 50000
nBurnin = 20000
nChains = 2
nThin = 5
#nThin = 1
numParticles = 50


# Function to simulate dynamic occupancy models
dynOccupancyModels <- function(nyears,
                               nsites,
                               nvisits,
                               fixedPars = list(alphaPhi = 2,
                                                betaPhi = - 1.5),
                               hyperParsSig = list(alphaP = 0.8,
                                                   betaP = 1,
                                                   alphaGamma = 0.3,
                                                   betaGamma = 0.5,
                                                   alphaPsi = 0.6,
                                                   betaPsi = 1)){

  # Simulating covariates
  set.seed(2020)

  windSpeed <- array(rnorm(nyears*nsites*nvisits, 0, 1),
                     dim = c(nsites, nvisits, nyears))
  elevation <- rnorm(nsites, 0, 1)
  springPrecipitation <- matrix(rnorm(nyears*nsites, 0, 1),
                                nrow = nsites,
                                ncol = nyears)
  sizeOfBeak <- matrix(rnorm(nyears*nsites, 0,1),
                       nrow = nsites,
                       ncol = nyears)

  # Simulating parameters
  alphaP <- 1
  betaP <- -2 #rnorm(1, mean = 0, sd = hyperParsSig$betaP)
  alphaPsi <- 1  #rnorm(1, mean = 0, sd = hyperParsSig$alphaPsi)
  betaPsi <- -1 #rnorm(1, mean = 0, sd = hyperParsSig$betaPsi)
  alphaGamma <- -2 #rnorm(1, mean = 0, sd = hyperParsSig$alphaGamma)
  betaGamma <- 1 #rnorm(1, mean = 0, sd = hyperParsSig$betaGamma)
  # alphaP <- rnorm(1, mean = 0, sd = hyperParsSig$alphaP)
  # betaP <- rnorm(1, mean = 0, sd = hyperParsSig$betaP)
  # alphaPsi <- rnorm(1, mean = 0, sd = hyperParsSig$alphaPsi)
  # betaPsi <- rnorm(1, mean = 0, sd = hyperParsSig$betaPsi)
  # alphaGamma <- rnorm(1, mean = 0, sd = hyperParsSig$alphaGamma)
  # betaGamma <- rnorm(1, mean = 0, sd = hyperParsSig$betaGamma)

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
  persistenceProb  <- colonisationProb <- matrix(NA,  nrow = nsites, ncol = nyears)
  for(year.tag in 1:nyears){
    for(site.tag in 1:nsites){
      persistenceProb[site.tag,  year.tag] <- plogis(alphaPsi + betaPsi* springPrecipitation[site.tag, year.tag])
       colonisationProb[site.tag,  year.tag] <- plogis(alphaGamma + betaGamma* sizeOfBeak[site.tag, year.tag])
    }
  }

  #colonisationProb <- 0.5

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

# Simulating and saving data
simData <- dynOccupancyModels(nyears = 20,
                              nsites = 100,
                              nvisits = 6)
#
#save(simData, file = "simDataDynamicOccupancy.RData")

#load("simDataDynamicOccupancy.RData")

# NIMBLE code
# Try to do this first before updating, so I can use it for the bootstrap updating
#for(i in seq_along(iNodePrevAll)){
i = 2
iNodePrev <- iNodePrevAll[i]
dynOccupancyCode <- nimbleCode({

  # Prior distributions of hyperparameters
  # alphaPSig ~ T(dnorm(0, 1), 0.001, )
  # betaPSig ~ T(dnorm(0, 1), 0.001, )
  # alphaPsiSig ~ T(dnorm(0, 1), 0.001, )
  # betaPsiSig ~ T(dnorm(0, 1), 0.001, )

  # alphaPSig ~ dunif(0.01, 2)
  # betaPSig ~ dunif(0.01, 2)
  # alphaPsiSig ~ dunif(0.01, 2)
  # betaPsiSig ~ dunif(0.01, 2)
  # alphaGammaSig ~ dunif(0.01, 2)
  # betaGammaSig ~ dunif(0.01, 2)

  alphaPhi ~ dnorm(0, sd = 1)
  betaPhi ~ dnorm(0, sd = 1)

  # Prior distributions
  # alphaP ~ dnorm(mean = 0, sd = alphaPSig)
  # betaP ~ dnorm(mean = 0, sd =  betaPSig)
  # alphaPsi ~ dnorm(mean = 0, sd =  alphaPsiSig)
  # betaPsi ~ dnorm(mean = 0, sd =  betaPsiSig)
  # alphaGamma ~ dnorm(mean = 0, sd =  alphaGammaSig)
  # betaGamma ~ dnorm(mean = 0, sd =  betaGammaSig)
  #
  alphaP ~ dnorm(mean = 0, sd = 1)
  betaP ~ dnorm(mean = 0, sd =  1)
  alphaPsi ~ dnorm(mean = 0, sd =  1)
  betaPsi ~ dnorm(mean = 0, sd =  1)
  alphaGamma ~ dnorm(mean = 0, sd =  1)
  betaGamma ~ dnorm(mean = 0, sd =  1)

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
       logit(colonisationProb[site.tag,  year.tag]) <- alphaGamma + betaGamma* sizeOfBeak[site.tag, year.tag]
    }
  }

  #colonisationProb ~ dunif(0.01, 1)


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

  # create a dummy Variable for sampler Times per iteration
  for(i in 1:4){
    samplerTimes[i] <- 1000
  }

})


# ## define data, constants, and initial values
z <- apply(simData$y, c(1, 3), function(x){
  ifelse(sum(x) > 0, 1, NA)
})
#
#
#   #####################
#   #   Reduced Model
#   ########################
#
data <- list(
  y = simData$y[,,-c((iNodePrev +1):nyears)],
  windSpeed = simData$covariates$windSpeed[,,-c((iNodePrev +1):nyears)],
  elevation = simData$covariates$elevation,
  springPrecipitation = simData$covariates$springPrecipitation[,-c((iNodePrev +1):nyears)],
  sizeOfBeak = simData$covariates$sizeOfBeak,
  z = z[,-c((iNodePrev +1):nyears)]
)
#
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
  #colonisationProb = 0.01,
  alphaPhi=  rnorm(1),
  betaPhi= rnorm(1)
)
#
#
#
#   # NIMBLE model for reduced model
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
                                             MCMCconfiguration = list(target = c(#'alphaPSig', 'betaPSig',
                                                                                 #'alphaPsiSig', 'betaPsiSig',
                                                                                 'alphaGamma', 'betaGamma',
                                                                                 #'alphaGammaSig', 'betaGammaSig',
                                                                                 'alphaPhi', 'betaPhi',
                                                                                 'alphaP', 'betaP',
                                                                                 'alphaPsi', 'betaPsi'
                                             ),
                                             additionalPars = c( "psi.fs",
                                                                 "samplerTimes"
                                             ),
                                             n.iter = nIterations,
                                             n.chains = nChains,
                                             n.burnin = nBurnin,
                                             n.thin = nThin)
)
#save results
# save(example2ReducedModelTrue, file= paste0("example4ReducedBootstrapTrueM",numParticles,"Ind",iNodePrev,".RData"))

#load(paste0("example4ReducedBootstrapTrue",i,".RData"))
replicatedReducedModel <- nimMCMCSMCupdates::replicateMCMCchains(example2ReducedModelTrue,
                                                                 N = nIterations,
                                                                 nChains = nChains)
################
# Updated Model
################
#load(paste0("example4ReducedBootstrapTrueM1000Ind",iNodePrev,".RData"))
# z <- simData$z
# z1 <- apply(simData$y, c(1, 3), function(x){
#   ifelse(sum(x) > 0, 1, NA)
# })
#
# z[,49:50] <- z1[, 49:50]

data <- list(
  y = simData$y,
  windSpeed = simData$covariates$windSpeed,
  elevation = simData$covariates$elevation,
  springPrecipitation = simData$covariates$springPrecipitation,
  sizeOfBeak = simData$covariates$sizeOfBeak,
  z = z
)
#
constants <- list(
  nyears = nyears,
  nsites = nsites,
  nvisits = nvisits#,
 # alphaPhi  = example2ReducedModelTrue$summary$all.chains["alphaPhi",1 ],
  #betaPhi  = example2ReducedModelTrue$summary$all.chains["betaPhi",1 ]
)


inits <- list(
 # alphaPSig = example2ReducedModelTrue$summary$all.chains["alphaPSig",1 ],
 # betaPSig  = example2ReducedModelTrue$summary$all.chains["betaPSig",1 ],
 # alphaPsiSig = example2ReducedModelTrue$summary$all.chains["alphaPsiSig",1 ],
 # betaPsiSig  = example2ReducedModelTrue$summary$all.chains["betaPsiSig",1 ],
  alphaGamma = example2ReducedModelTrue$summary$all.chains["alphaGamma",1 ],
 betaGamma = example2ReducedModelTrue$summary$all.chains["betaGamma",1 ],
  #alphaPhi  = 0,
  #betaPhi  = 0,
  alphaP =  example2ReducedModelTrue$summary$all.chains["alphaP",1 ],
  betaP =  example2ReducedModelTrue$summary$all.chains["betaP",1 ],
  alphaPsi=  example2ReducedModelTrue$summary$all.chains["alphaPsi",1 ],
  alphaPhi  = example2ReducedModelTrue$summary$all.chains["alphaPhi",1 ],
  betaPhi  = example2ReducedModelTrue$summary$all.chains["betaPhi",1 ],
  betaPsi= example2ReducedModelTrue$summary$all.chains["betaPsi",1 ]#,
 # colonisationProb = example2ReducedModelTrue$summary$all.chains["colonisationProb",1 ]
)

# NIMBLE model for updated model
newModelUpdated <- nimbleModel(dynOccupancyCode,
                               data = data,
                               constants = constants,
                               inits = inits,
                               check = FALSE)
#
#
#
#   # ## Checking predictions
#   # notDataNodes <- newModelUpdated$expandNodeNames("z")[!newModelUpdated$isData("z")]
#   # notDataNode1 <- notDataNodes[1:44]
#   # notDataNode2 <- notDataNodes[45: length(notDataNodes)]
#   # newModelUpdated$persistenceProb
#   # newModelUpdated$psi.fs
#   # newModelUpdated$simulate(nodes = notDataNode1)
#   # newModelUpdated$calculate("psi.fs")
#
#   # Fitting the model
example2UpdatedModelTrue <- spartaNimUpdates(model = newModelUpdated, #nimble model
                                             reducedModel = newModelReduced,
                                             latent = "z", #latent variable
                                             nParFiltRun = numParticles,
                                             #nParFiltRun = 1000,
                                             pfType = "bootstrap",
                                             #propCov = c(1,1,1,1,1,1,0.01)*diag(7),
                                             mcmcScale = 1,
                                             adaptiveSampler = TRUE,
                                             adaptScaleOnly = TRUE,
                                             initsList = inits,
                                             extraVars = c('alphaPhi', 'betaPhi'),
                                             MCMCconfiguration = list(target = c(#'alphaPSig', 'betaPSig',
                                                                                 #'alphaPsiSig', 'betaPsiSig',
                                                                                 'alphaGamma', 'betaGamma',
                                                                                # 'alphaPhi', 'betaPhi',
                                                                                # 'alphaGammaSig', 'betaGammaSig',
                                                                                 'alphaP', 'betaP',
                                                                                 'alphaPsi', 'betaPsi'),
                                                                      additionalPars = c("z", "psi.fs", "samplerTimes"),
                                                                      n.iter = nIterations,
                                                                      n.chains = 1,
                                                                      n.burnin = nBurnin,
                                                                      n.thin = nThin
                                             ),  #saved loglikelihoods from reduced model
                                             postReducedMCMC = replicatedReducedModel,# MCMC summary to use as initial values
                                             pfControl = list(saveAll = TRUE,
                                                              timeIndex = 2,
                                                              smoothing = TRUE,
                                                              mcmc = TRUE,
                                                              M = nyears - iNodePrev,
                                                              iNodePrev = iNodePrev),
                                             nCores = 1
)


example2UpdatedModelTrue$summary$all.chains[grepl("psi", rownames(example2UpdatedModelTrue$summary$all.chains)),1 ]
example2UpdatedModelTrue$samples$chain1[, grepl("psi", colnames(example2UpdatedModelTrue$samples$chain1))]

model = newModelUpdated #nimble model
reducedModel = newModelReduced
latent = "z" #latent variable
nParFiltRun = numParticles
#nParFiltRun = 1000,
pfType = "bootstrap"
#propCov = c(1,1,1,1,1,1,0.01)*diag(7),
mcmcScale = 0.5
adaptiveSampler = TRUE
adaptScaleOnly = TRUE
extraVars = NULL#c('alphaPhi', 'betaPhi'),
MCMCconfiguration = list(target = c('alphaPSig', 'betaPSig',
                                    'alphaPsiSig', 'betaPsiSig',
                                    'alphaGamma', 'betaGamma',
                                    'alphaPhi', 'betaPhi',
                                    'alphaGammaSig', 'betaGammaSig',
                                    'alphaP', 'betaP',
                                    'alphaPsi', 'betaPsi'),
                         additionalPars = c("z", "psi.fs", "samplerTimes"),
                         n.iter = (nIterations - nBurnin)/nThin,
                         n.chains = nChains,
                         n.burnin = 10,
                         n.thin = 1
)  #saved loglikelihoods from reduced model
postReducedMCMC = example2ReducedModelTrue# MCMC summary to use as initial values
pfControl = list(saveAll = TRUE,
                 timeIndex = 2,
                 smoothing = TRUE,
                 mcmc = TRUE,
                 M = nyears - iNodePrev,
                 iNodePrev = iNodePrev)
nCores = 1
#
save(example2UpdatedModelTrue, file= paste0("example4UpdatedBootstrapTrueM",numParticles,"Ind",iNodePrev,".RData"))
#
#
#   # Fitting the model
example2UpdatedModelTrue <- spartaNimUpdates(model = newModelUpdated, #nimble model
                                             reducedModel = newModelReduced,
                                             latent = "z", #latent variable
                                             nParFiltRun = numParticles,
                                             #nParFiltRun = 1000,
                                             pfType = "auxiliary",
                                             #propCov = c(1,1,1,1,1,1,0.01)*diag(7),
                                             mcmcScale = 0.1,
                                             adaptiveSampler = TRUE,
                                             adaptScaleOnly = TRUE,
                                             extraVars = c('alphaPhi', 'betaPhi'),
                                             MCMCconfiguration = list(target = c(#'alphaPSig', 'betaPSig',
                                               #'alphaPsiSig', 'betaPsiSig',
                                               'alphaGamma', 'betaGamma',
                                               #'alphaPhi', 'betaPhi',
                                               # 'alphaGammaSig', 'betaGammaSig',
                                               'alphaP', 'betaP',
                                               'alphaPsi', 'betaPsi'),
                                               additionalPars = c("z", "psi.fs", "samplerTimes"),
                                               n.iter = (nIterations - nBurnin)/nThin,
                                               n.chains = nChains,
                                               n.burnin = 1000,
                                               n.thin = 1
                                             ),  #saved loglikelihoods from reduced model
                                             postReducedMCMC = example2ReducedModelTrue,# MCMC summary to use as initial values
                                             pfControl = list(saveAll = TRUE,
                                                              timeIndex = 2,
                                                              smoothing = TRUE,
                                                              mcmc = TRUE,
                                                              M = nyears - iNodePrev,
                                                              iNodePrev = iNodePrev),
                                             nCores = 1
)


save(example2UpdatedModelTrue, file= paste0("example4UpdatedAuxiliaryTrueM",numParticles,"Ind",iNodePrev,".RData"))
# #}
#
#   dynOccupancyCode <- nimbleCode({
#
#     # Prior distributions of hyperparameters
#     alphaPSig ~ T(dnorm(0, 0.1), 0.001, )
#     betaPSig ~ T(dnorm(0, 0.1), 0.001, )
#     alphaPsiSig ~ T(dnorm(0, 0.1), 0.001, )
#     betaPsiSig ~ T(dnorm(0, 0.1), 0.001, )
#
#
#     alphaPhi ~ dnorm(0, sd = 10)
#     betaPhi ~ dnorm(0, sd = 10)
#
#     # Prior distributions
#     alphaP ~ dnorm(mean = 0, sd = alphaPSig)
#     betaP ~ dnorm(mean = 0, sd =  betaPSig)
#     alphaPsi ~ dnorm(mean = 0, sd =  alphaPsiSig)
#     betaPsi ~ dnorm(mean = 0, sd =  betaPsiSig)
#
#     # Detection Probability
#     for(year.tag in 1:nyears){
#       for(site.tag in 1:nsites){
#         for(visit.tag in 1:nvisits){
#           logit(detectProb[site.tag, visit.tag, year.tag]) <- alphaP + betaP* windSpeed[site.tag, visit.tag, year.tag]
#         }
#       }
#     }
#
#     # Initial occupancy probability psi
#     for(site.tag in 1:nsites){
#       logit(initOccuProb[site.tag]) <- alphaPhi + betaPhi*elevation[site.tag]
#     }
#
#     # Persistence and colonisation probability
#     for(year.tag in 1:nyears){
#       for(site.tag in 1:nsites){
#         logit(persistenceProb[site.tag,  year.tag]) <- alphaPsi + betaPsi* springPrecipitation[site.tag, year.tag]
#         # logit(colonisationProb[site.tag,  year.tag]) <- alphaGamma[year.tag] + betaGamma[year.tag]* sizeOfBeak[site.tag, year.tag]
#       }
#     }
#
#     colonisationProb ~ dunif(0.01, 1)
#
#
#     # Initial Presence/Absence
#     for(site.tag in 1:nsites){
#       z[site.tag, 1] ~ dbin(prob = initOccuProb[site.tag], 1)
#     }
#
#     # True presence/absence
#     for(year.tag in 2:nyears){
#       for(site.tag in 1:nsites){
#         z[site.tag, year.tag] ~ dbin(prob = (z[site.tag, (year.tag -1)] * persistenceProb[site.tag,  year.tag] + (1 - z[site.tag, (year.tag -1)])*colonisationProb), 1)
#       }
#     }
#
#     #observations
#     for(year.tag in 1:nyears){
#       for(site.tag in 1:nsites){
#         for(visit.tag in 1:nvisits){
#           y[site.tag, visit.tag, year.tag] ~ dbin(prob = z[site.tag, year.tag] * detectProb[site.tag, visit.tag, year.tag], 1)
#         }
#       }
#     }
#
#     # Derived quantities
#     for(year.tag in 1:nyears){
#       psi.fs[year.tag] <- sum(z[1:nsites, year.tag])/nsites
#     }
#
#   })
#
#
# # ## define data, constants, and initial values
# z <- apply(simData$y, c(1, 3), function(x){
#   ifelse(sum(x) > 0, 1, NA)
# })
#
#
data <- list(
  y = simData$y,
  windSpeed = simData$covariates$windSpeed,
  elevation = simData$covariates$elevation,
  springPrecipitation = simData$covariates$springPrecipitation,
  sizeOfBeak = simData$covariates$sizeOfBeak,
  z = z
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
  alphaGamma = rnorm(1, mean = 0, sd = 1),
  betaGamma = rnorm(1, mean = 0, sd = 1),
  alphaPhi  = 0,
  betaPhi  = 0,
  alphaP =  rnorm(1, mean = 0, sd = 1),
  betaP =  rnorm(1, mean = 0, sd = 1),
  alphaPsi=  rnorm(1, mean = 0, sd = 1),
  betaPsi= rnorm(1, mean = 0, sd = 1)#,
  #colonisationProb = 0.01
)

# NIMBLE model for reduced model
baselineModel <- nimbleModel(dynOccupancyCode,
                             data = data,
                             constants = constants,
                             inits = inits,
                             check = FALSE)

# Fitting the model
baselineModelEst <- spartaNimWeights(model = baselineModel,
                                     latent = "z",
                                     nParFiltRun = numParticles,
                                     mcmc = TRUE,
                                     pfType = pfTypeRun,
                                     #pfType = pfTypeRun,
                                     pfControl = list(saveAll = TRUE,
                                                      smoothing = TRUE,
                                                      timeIndex = 2),
                                     MCMCconfiguration = list(target = c(#'alphaPSig', 'betaPSig',
                                                                         #'alphaPsiSig', 'betaPsiSig',
                                                                         'alphaPhi', 'betaPhi',
                                                                         'alphaP', 'betaP',
                                                                         'alphaPsi', 'betaPsi',
                                                                         'alphaGamma', 'betaGamma'
                                                                         #"colonisationProb"
                                     ),
                                     additionalPars = c( "psi.fs"
                                     ),
                                     n.iter = nIterations,
                                     n.chains = nChains,
                                     n.burnin = nBurnin,
                                     n.thin = nThin)
)
# #save results
# save(baselineModelEst, file= "baselineModel.RData")
