# This script fits

# load packages
library(nimble)
library(nimbleSMC)
library(nimMCMCSMCupdates)

### Parameters for running script
numParticles <- 20
pfTypeRun <- "bootstrap"
nIterations = 32000
nChains = 3
nThin = 2
nBurnin = 30000

# Extract values
# simData = out$model$data()
#
# mcmcOutput <- out$BUGSoutput$sims.list
# inits <- out$model$state()
# save(simData, file = "simData.RData")
# save(inits, file = "inits.RData")
# Load data and initial values from SPARTA function
setwd("/Volumes/kwakupa/particleFilters/Example3")
load("simDataFullModel.RData")
load("initsFullModel.RData")

## Extracting NIMBLE code


# Defining NIMBLE model
occupancyCode <- nimbleCode({
  # JAGS code for SPARTA model plus random walk prior
  # on the year effect of the state model + intercept + halfcauchy hyperpriors

  # State model
  for (i in 1:nsite){
    z[i,1] ~ dbern(muZ[i,1])
    logit(muZ[i,1]) <- a[1] + eta[i]
    for (t in 2:nyear){
      z[i,t] ~ dbern(muZ[i,t])
      logit(muZ[i,t]) <- a[t] + eta[i]
    }}


  # State model priors
  a[1] ~ dnorm(0, 0.001)
  for(t in 2:nyear){
    a[t] ~ dnorm(a[t-1], tau.a)
  }

  tau.a <- 1/(sd.a * sd.a)

  #sd.a ~ dt(0, 1, 1)T(0,)
  sd.a ~  T(dt(0, 1, 1),0,)

  for (i in 1:nsite) {
    eta[i] ~ dnorm(0, tau2)
  }

  tau2 <- 1/(sigma2 * sigma2)
  sigma2 ~  T(dt(0, 1, 1),0,)

  # Observation model priors
  for (t in 1:nyear) {
    alpha.p[t] ~ dnorm(mu.lp, tau.lp)
  }

  mu.lp ~ dnorm(0, 0.01)

  tau.lp <- 1 / (sd.lp * sd.lp)
  sd.lp ~ T(dt(0, 1, 1),0,)

  # Derived parameters"
  for (t in 1:nyear) {
    psi.fs[t] <- sum(z[1:nsite, t])/nsite
  }

  dtype2.p ~ dnorm(0, 0.01)
  dtype3.p ~ dnorm(0, 0.01)

  ### Observation Model
  for(j in 1:nvisit) {
    y[j] ~ dbern(Py[j])
    Py[j]<- z[Site[j],Year[j]]*p[j]
    logit(p[j]) <-  alpha.p[Year[j]] + dtype2.p*DATATYPE2[j] + dtype3.p*DATATYPE3[j]
  }

  # deviance <- 1
})

#####################
### Reduced Model
###################
simData <- simDataFullModel
iNodePrev = 50 # year 2021
Year = simData$Year
index <- Year %in% seq(1:iNodePrev) # index to subset the years
reducedVisits <- sum(index) # number of visit for reduced years
Site = simData$Site
Year = simData$Year
y = simData$y

z <- matrix(NA, nrow = simData$nsite, ncol = simData$nyear)
for(j in 1:simData$nvisit){
  if(y[j] == 1){
    z[Site[j], Year[j]] <- 1
  }
}

data <- list(
  y = simData$y[index],
  DATATYPE2 = simData$DATATYPE2[index],
  DATATYPE3 = simData$DATATYPE3[index],
  z = z[, 1:iNodePrev]
)

constants <- list(
  nyear = iNodePrev,
  nsite = simData$nsite,
  nvisit = reducedVisits,
  Site = simData$Site[index],
  Year = simData$Year[index]
)


#inits <- out$model$state()
inits <- inits[[1]]
inits$a <- inits$a[1:iNodePrev]
inits$alpha.p <- inits$alpha.p[1:iNodePrev]
init <- inits[1:9]


# ## build the model
reducedModel <- nimbleModel(occupancyCode,
                            data = data,
                            constants = constants,
                            inits = init,
                            check = FALSE)


targetsForSMC <- c("a","eta",
                   "sd.a", "sigma2",
                   "mu.lp", "sd.lp",
                   "dtype2.p", "dtype3.p", "alpha.p")


# example2ReducedModelTrue <- spartaNimWeights(model = reducedModel,
#                                              latent = "z",
#                                              nParFiltRun = numParticles,
#                                              mcmc = TRUE,
#                                              #pfType = pfTypeRun,
#                                              #pfType = pfTypeRun,
#                                              pfControl = list(saveAll = TRUE,
#                                                               smoothing = TRUE,
#                                                               timeIndex = 2),
#                                              MCMCconfiguration = list(target = targetsForSMC,
#                                                                       additionalPars = c("z",
#                                                                                          "psi.fs"),
#                                                                       n.iter = 200,
#                                                                       n.chains = 2,
#                                                                       n.burnin = 50,
#                                                                       n.thin = 1)
# )

#save(example2ReducedModelTrue, file = "reducedModelResults.RData")

load("reducedModelResults.RData")


#example2ReducedModelTrue <- replicateMCMCchains(example2ReducedModelTrue, N = 5000, nChains = nChains)
#####################
# Updated model
#####################
#index <-Year %in% 1:7 # index to subset the years

#reducedVisits <- sum(index) # number of visit for reduced years

occupancyCodeUpdated <- nimbleCode({
  # JAGS code for SPARTA model plus random walk prior
  # on the year effect of the state model + intercept + halfcauchy hyperpriors

  # State model
  for (i in 1:nsite){
    z[i,1] ~ dbern(muZ[i,1])
    logit(muZ[i,1]) <- a[1] + eta[i]
    for (t in 2:nyear){
      z[i,t] ~ dbern(muZ[i,t])
      logit(muZ[i,t]) <- a[t] + eta[i]
    }}


  # State model priors
  a[1] ~ dnorm(0, 0.001)
  for(t in 2:nyear){
    a[t] ~ dnorm(a[t-1], tau.a)
  }

  tau.a <- 1/(sd.a * sd.a)

  #sd.a ~ dt(0, 1, 1)T(0,)
  sd.a ~  T(dt(0, 1, 1),0,)

  #for (i in 1:nsite) {
  #  eta[i] ~ dnorm(0, tau2)
  #}

  #tau2 <- 1/(sigma2 * sigma2)
  #sigma2 ~  T(dt(0, 1, 1),0,)

  # Observation model priors
  for (t in 1:nyear) {
    alpha.p[t] ~ dnorm(mu.lp, tau.lp)
  }

  mu.lp ~ dnorm(0, 0.01)

  tau.lp <- 1 / (sd.lp * sd.lp)
  sd.lp ~ T(dt(0, 1, 1),0,)

  # Derived parameters"
  for (t in 1:nyear) {
    psi.fs[t] <- sum(z[1:nsite, t])/nsite
  }

  dtype2.p ~ dnorm(0, 0.01)
  dtype3.p ~ dnorm(0, 0.01)

  ### Observation Model
  for(j in 1:nvisit) {
    y[j] ~ dbern(Py[j])
    Py[j]<- z[Site[j],Year[j]]*p[j]
    logit(p[j]) <-  alpha.p[Year[j]] + dtype2.p*DATATYPE2[j] + dtype3.p*DATATYPE3[j]
  }

  # deviance <- 1
})

x <- example2ReducedModelTrue
extractNames <- rownames(x$summary$all.chains)[grepl("eta", rownames(x$summary$all.chains))]
eta <- x$summary$all.chains[extractNames, 1]

data1 <- list(
  y = simData$y,
  DATATYPE2 = simData$DATATYPE2,
  DATATYPE3 = simData$DATATYPE3,
  z = z,
  eta = eta
)

constants1 <- list(
  nyear = simData$nyear,#(simData$nyear - iNodePrev + 1),
  #nyear = iNodePrev,
  nsite = simData$nsite,
  nvisit = simData$nvisit,
  Site = simData$Site, #[index],
  Year = simData$Year #- iNodePrev + 1)
)


load("initsFullModel.RData")
inits1 <- inits[[1]]
#inits1$a <- inits1$a#[(iNodePrev) : simData$nyear]
#inits1$alpha.p <- inits1$alpha.p[c((iNodePrev) : simData$nyear)]
init <- inits1[-c(5,10)]
#init$a <- init$a
#init$alpha.p <- init$alpha.p

# ## build the model
updatedModel <- nimbleModel(occupancyCodeUpdated,
                            data = data1,
                            constants = constants1,
                            inits = init,
                            check = FALSE)


### Parameters to monitor
latent = "z"
targetsForSMC <- c("a",
                   "sd.a",
                   "mu.lp", "sd.lp",
                   "dtype2.p", "dtype3.p", "alpha.p")


### run state space model with reduced model
#nimbleOptions(MCMCusePredictiveDependenciesInCalculations = TRUE)
#nimbleOptions(MCMCorderPosteriorPredictiveSamplersLast = FALSE)

example2UpdatedModelTrue <- spartaNimUpdates(model = updatedModel , #nimble model
                                             reducedModel = reducedModel,
                                             latent = "z", #latent variable
                                             nParFiltRun = numParticles,
                                             #nParFiltRun = 1000,
                                             pfType = "bootstrap",
                                             extraVars = c("a", "alpha.p"),
                                             mcmcScale = 0.2,
                                             MCMCconfiguration = list(target = targetsForSMC,
                                                                      additionalPars = c("z",
                                                                                         "psi.fs",
                                                                                         "eta"),
                                                                      n.iter = 100,#(nIterations - nBurnin)/nThin,
                                                                      n.chains = 2,#nChains,
                                                                      n.burnin = 1,
                                                                      n.thin = nThin),  #saved loglikelihoods from reduced model
                                             postReducedMCMC = example2ReducedModelTrue,# MCMC summary to use as initial values
                                             pfControl = list(saveAll = TRUE,
                                                              timeIndex = 2,
                                                              smoothing = TRUE,
                                                              #mcmc = TRUE,
                                                              M = 18 - iNodePrev,
                                                              iNodePrev = iNodePrev),
                                             nCores = 1
)


save(example2UpdatedModelTrue, file = "updatedModelResultsAuxiliary1.RData")

### run state space model with reduced model
nimbleOptions(MCMCusePredictiveDependenciesInCalculations = TRUE)
nimbleOptions(MCMCorderPosteriorPredictiveSamplersLast = FALSE)

example2UpdatedModelTrue <- spartaNimUpdates(model = updatedModel , #nimble model
                                             reducedModel = reducedModel,
                                             latent = "z", #latent variable
                                             nParFiltRun = numParticles,
                                             #nParFiltRun = 1000,
                                             pfType = "bootstrap",
                                             extraVars = c("a", "alpha.p"),
                                             mcmcScale = 0.2,
                                             MCMCconfiguration = list(target = targetsForSMC,
                                                                      additionalPars = c("z",
                                                                                         "psi.fs",
                                                                                         "eta"),
                                                                      n.iter = 3000,#(nIterations - nBurnin)/nThin,
                                                                      n.chains = nChains,
                                                                      n.burnin = 1000,
                                                                      n.thin = nThin),  #saved loglikelihoods from reduced model
                                             postReducedMCMC = example2ReducedModelTrue,# MCMC summary to use as initial values
                                             pfControl = list(saveAll = TRUE,
                                                              timeIndex = 2,
                                                              smoothing = TRUE,
                                                              #mcmc = TRUE,
                                                              M = 18 - iNodePrev,
                                                              iNodePrev = iNodePrev)
)


save(example2UpdatedModelTrue, file = "updatedModelResultsBootstrap1.RData")


