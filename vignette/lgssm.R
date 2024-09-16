# This is the general function that is used to fit the model with MCMC
# and simulate the required data


# Load the required packages
library('nimble')
library(nimbleSMC)
library(nimMCMCSMCupdates)
nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
nimbleOptions(MCMCsaveHistory = TRUE)

# MCMC iterations set-up
nIterations = 50000
nChains = 2
nBurnin = 20000
nThin = 10

################################
#  Function to simulate data
###############################
sim2 <- function(a, b, c, mu0, t, seed){
  set.seed(seed)
  x <- y <- numeric(t)
  x[1] <- rnorm(1, mu0, sd = 1 )
  y[1] <- rnorm(1, x[1]*c, sd=  1)


  for(k in 2:t){
    x[k] <- rnorm(1, a*x[k -1] + b, sd = 1)
    y[k] <- rnorm(1, x[k]*c, sd= 1)
  }
  return(list(x=x, y=y))
}

# simulate data
simData <- sim2(0.5, 1,1,1,50,1)

# NIMBLE CODE
stateSpaceCode <- nimbleCode({
  x[1] ~ dnorm(mu, sd = 1)
  #x[1] <- 0
  y[1] ~ dnorm(c*x[1], sd = 1)
  for(i in 2:t){
    x[i] ~ dnorm(a*x[i-1] + b, sd = 1)
    y[i] ~ dnorm(x[i] * c, 1)
  }
  a ~ dunif(0, 0.999)
  b ~ dnorm(0, 1)
  c ~ dnorm(1,1)
  mu ~ dnorm(0, 1)


  #estimate biases
  ahat <- a - aTrue
  bhat <- b - bTrue
  chat <- c - cTrue
  muhat <- mu - muTrue

  # create a dummy Variable for sampler Times per iteration
  for(i in 1:4){
    samplerTimes[i] <- 1000
  }

})

# ## define data
data <- list(
  y = simData$y
)

# define constants
constants <- list(
  t = 50,
  aTrue = 0.5,
  bTrue = 1,
  cTrue = 1,
  muTrue = 1
)

#define initial values
inits <- list(
  a = 0.2,
  b = 0.8,
  mu= 0.8,
  c = 0.5
)


# ## build the model
stateSpaceModel <- nimbleModel(stateSpaceCode,
                               data = data,
                               constants = constants,
                               inits = inits,
                               check = FALSE)

data <- list(
  y = simData$y[-c(21:50)]
)

constants <- list(
  t = 20,
  aTrue = 0.5,
  bTrue = 1,
  cTrue = 1,
  muTrue = 1
)


#Define model for reduced model
newModelReduced <- nimbleModel(stateSpaceCode,
                               data = data,
                               constants = constants,
                               inits = inits,
                               check = FALSE)


# Fit reduced model with MCMC
example1ReducedModelTrue <- nimMCMCSMCupdates::spartaNimWeights(model = newModelReduced,
                                                                             latent = "x",
                                                                             nParFiltRun = numParticles,
                                                                             mcmc = TRUE,
                                                                             block = FALSE,
                                                                             pfType = pfTypeRun,
                                                                             MCMCconfiguration = list(target = c('a', 'b',  'c', "mu"),
                                                                                                      additionalPars = c("x", "ahat", "bhat","chat", "muhat", "samplerTimes"),
                                                                                                      n.iter = nIterations,
                                                                                                      n.chains = nChains,
                                                                                                      n.burnin = nBurnin,
                                                                                                      n.thin = nThin)
)


################
# Updated Model
################

# data for updated model
data <- list(
  y = simData$y
)


# constants for updated model
constants <- list(
  t = 50,
  aTrue = 0.5,
  bTrue = 1,
  cTrue = 1,
  muTrue = 1
)


inits <- list(
  a = example1ReducedModelTrue$summary$all.chains["a",1 ],
  b = example1ReducedModelTrue$summary$all.chains["b", 1],
  c = example1ReducedModelTrue$summary$all.chains["c",1 ],
  mu = example1ReducedModelTrue$summary$all.chains["mu", 1]
)


# Define model for updated model
newModelUpdated <- nimbleModel(stateSpaceCode,
                               data = data,
                               constants = constants,
                               inits = inits,
                               check = FALSE)


# Fit model with results from fitted reduced MCMC model fit
## set options to make history accessible

# Resample the MCMC results
replicatedReducedModel <- nimMCMCSMCupdates::replicateMCMCchains(example1ReducedModelTrue,
                                                   N = nIterations,
                                                   nChains = nChains)

example1UpdatedModelTrueBF <- nimMCMCSMCupdates::spartaNimUpdates(model = newModelUpdated, #nimble model
                                                                               reducedModel = newModelReduced,
                                                                               latent = "x", #latent variable
                                                                               nParFiltRun = 100,
                                                                               mcmcScale = 2,
                                                                              #extraVars = c("mu"),
                                                                               adaptiveSampler = TRUE,
                                                                                adaptScaleOnly = TRUE,
                                                                  initsList = inits,
                                                                               MCMCconfiguration = list(target = c('a', 'c', 'b', "mu"),
                                                                                                        additionalPars = c("x", "ahat", "chat", "bhat", "samplerTimes", "muhat"),
                                                                                                        n.iter = nIterations,
                                                                                                        n.chains = nChains,
                                                                                                        n.burnin = nBurnin,
                                                                                                        n.thin = nThin),  #saved loglikelihoods from reduced model
                                                                               postReducedMCMC = replicatedReducedModel,# MCMC summary to use as initial values
                                                                               pfControl = list(saveAll = TRUE,
                                                                                                smoothing = TRUE,
                                                                                                mcmc = TRUE,
                                                                                                M = 30,
                                                                                                iNodePrev = 20),
                                                                               nCores = 1
)




 example1UpdatedModelTrueAux <- nimMCMCSMCupdates::spartaNimUpdates(model = newModelUpdated, #nimble model
                                                                                reducedModel = newModelReduced,
                                                                                latent = "x", #latent variable
                                                                                nParFiltRun = 100,
                                                                                mcmcScale = 0.5,
                                                                                #propCov = diag(2)* c(0.1,1),
                                                                                pfType = "auxiliary",
                                                                                extraVars = NULL,
                                                                                MCMCconfiguration = list(target = c('a', 'c'),
                                                                                                         additionalPars = c("x", "ahat", "chat"),
                                                                                                         n.iter = (nIterations - nBurnin)/nThin,
                                                                                                         n.chains = nChains,
                                                                                                         n.burnin = 0,
                                                                                                         n.thin = 1),  #saved loglikelihoods from reduced model
                                                                                postReducedMCMC = example1ReducedModelTrue,# MCMC summary to use as initial values
                                                                                pfControl = list(saveAll = TRUE,
                                                                                                 smoothing = TRUE,
                                                                                                 mcmc = TRUE,
                                                                                                 lookahead = "mean",
                                                                                                 M = 30,
                                                                                                 iNodePrev = 20),
                                                                   nCores = 1
)















  ############
  # Baseline Model
  ##############
  # ## define data, constants, and initial values
  data <- list(
    y = simData$y
  )
  constants <- list(
    t = nyears,
    aTrue = 0.5,
    bTrue = 0.5,
    cTrue = 1
  )
  inits <- list(
    a = 0.2,
    b = 1,
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




  #Define a new model
  newModel <- stateSpaceModel$newModel(replicate = TRUE)


  # Function to run the baseline model


  # Fit baseline model with MCMC
  baselineModel[[1]] <- nimMCMCSMCupdates::baselineSpartaEstimation(model = newModel,
                                                                    latent = "x",
                                                                    pfType = "bootstrap",
                                                                    nParFiltRun = numParticles,
                                                                    pfControl = list(saveAll = TRUE,
                                                                                     lookahead = "mean",
                                                                                     smoothing = FALSE),
                                                                    MCMCconfiguration = list(target = c('a', 'c'),
                                                                                             additionalPars = c("x", "ahat", "chat"),
                                                                                             n.iter = nIterations,
                                                                                             n.chains = nChains,
                                                                                             n.burnin = nBurnin,
                                                                                             n.thin = nThin))


  baselineModel[[2]] <- nimMCMCSMCupdates::baselineSpartaEstimation(model = newModel,
                                                                    latent = "x",
                                                                    pfType = "auxiliary",
                                                                    nParFiltRun = numParticles,
                                                                    pfControl = list(saveAll = TRUE,
                                                                                     lookahead = "mean",
                                                                                     smoothing = FALSE),
                                                                    MCMCconfiguration = list(target = c('a', 'c'),
                                                                                             additionalPars = c("x", "ahat", "chat"),
                                                                                             n.iter = nIterations,
                                                                                             n.chains = nChains,
                                                                                             n.burnin = nBurnin,
                                                                                             n.thin = nThin))




  # Fit baseline model with SMC
  baselineModel[[3]]  <- nimMCMCSMCupdates::spartaNimWeights(model = stateSpaceModel,
                                                             latent = "x",
                                                             nParFiltRun = numParticles,
                                                             mcmc = TRUE,
                                                             block = TRUE,
                                                             pfType = pfTypeRun,
                                                             MCMCconfiguration = list(target = c('a', 'c'),
                                                                                      additionalPars = c("x", "ahat", "chat"),
                                                                                      n.iter = nIterations,
                                                                                      n.chains = nChains,
                                                                                      n.burnin = nBurnin,
                                                                                      n.thin = nThin)
  )






#Function to simulate data
sim2 <- function(a, b, c, mu0, t, seed){
  set.seed(seed)
  x <- y <- numeric(t)
  x[1] <- rnorm(1, mu0, sd = 1 )
  y[1] <- rnorm(1, x[1]*c, sd=  1)


  for(k in 2:t){
    x[k] <- rnorm(1, a*x[k -1] + b, sd = 1)
    y[k] <- rnorm(1, x[k]*c, sd= 1)
  }
  return(list(x=x, y=y))
}
