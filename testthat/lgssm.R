nyears = 50
numParticles = 5000
#Function to simulate data
sim2 <- function(a, b, c, t, seed){
  set.seed(seed)
  x <- y <- numeric(t)
  x[1] <- rnorm(1, a, sd = 1 )
  y[1] <- rnorm(1, c*x[1], 1)

  for(k in 2:t){
    x[k] <- rnorm(1, a*x[k -1], sd = 1)
    y[k] <- rnorm(1, x[k]*c, 1)
  }
  return(list(x=x, y=y))
}

simData = sim2(0.5, 1,1, 50, 1)
simData


stateSpaceCode <- nimbleCode({
  x[1] ~ dnorm(a, sd = 1)
  #x[1] <- 0
  y[1] ~ dnorm(c*x[1], sd = 1)
  for(i in 2:t){
    x[i] ~ dnorm(a*x[i-1] , sd = 1)
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
  t = 50,
  aTrue = 0.5,
  bTrue = 0.5,
  cTrue = 1,
  dTrue = 1
)

#define initial values
inits <- list(
  a = 0.2,
  b = 0.8,
  # mu0= 0.8,
  c = 1
)

# ## build the model
stateSpaceModel <- nimbleModel(stateSpaceCode,
                               data = data,
                               constants = constants,
                               inits = inits,
                               check = FALSE)


iNodePrev = 40
iNodetag = 1
pfTypeRun = "auxiliary"


  data <- list(
    y = simData$y[-c((iNodePrev[iNodetag]+1):50)]
  )

  constants <- list(
    t = iNodePrev[iNodetag],
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
  nIterations = 100
  nBurnin = 10
  nChains = 2
  nThin = 1
  # Fit reduced model with MCMC
  example1ReducedModelTrue  <- nimMCMCSMCupdates::spartaNimWeights(model = newModelReduced,
                                                                               latent = "x",
                                                                               nParFiltRun = numParticles,
                                                                               mcmc = TRUE,
                                                                               block = FALSE,
                                                                               pfType = pfTypeRun,
                                                                               MCMCconfiguration = list(target = c('a',  'c'),
                                                                                                        additionalPars = c("x", "ahat",  "chat"),
                                                                                                        n.iter = nIterations,
                                                                                                        n.chains = nChains,
                                                                                                        n.burnin = nBurnin,
                                                                                                        n.thin = nThin)
  )





  ################
  # Updated Model
  ################
  #message(paste("Running updated model for iNodePrev = ", iNodePrev[iNodetag], "and a = ", aVars[aVarstag]))

  # data for updated model
  data <- list(
    y = simData$y#[c((iNodePrev[iNodetag]):50)]
  )

  # constants for updated model
  constants <- list(
    t = 50 ,#- iNodePrev[iNodetag] + 1,
    aTrue = 0.5,
    bTrue = 0.5,
    cTrue = 1
  )

  # Define model for updated model
  newModelUpdated <- nimbleModel(stateSpaceCode,
                                 data = data,
                                 constants = constants,
                                 inits = inits,
                                 check = FALSE)

  # Fit model with results from fitted reduced MCMC model fit
  example1UpdatedModelTrue <- nimMCMCSMCupdates::spartaNimUpdates(model = newModelUpdated, #nimble model
                                                                               reducedModel = newModelReduced,
                                                                               latent = "x", #latent variable
                                                                               nParFiltRun = numParticles,
                                                                               propCov = diag(2)* c(0.1, 1),
                                                                               pfType = pfTypeRun,
                                                                               MCMCconfiguration = list(target = c('a',  'c'),
                                                                                                        additionalPars = c("x", "ahat", "chat"),
                                                                                                        n.iter = (nIterations - nBurnin)/nThin,
                                                                                                        n.chains = nChains,
                                                                                                        n.burnin = 0,
                                                                                                        n.thin = 1),  #saved loglikelihoods from reduced model
                                                                               postReducedMCMC = example1ReducedModelTrue,# MCMC summary to use as initial values
                                                                               pfControl = list(saveAll = TRUE,
                                                                                                smoothing = TRUE,
                                                                                                mcmc = TRUE,
                                                                                                M = nyears - iNodePrev[iNodetag],
                                                                                                iNodePrev = iNodePrev[iNodetag])
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
                                                                  pfType = pfTypeRun,
                                                                  nParFiltRun = numParticles,
                                                                  MCMCconfiguration = list(target = c('a', 'b', 'c'),
                                                                                           additionalPars = c("x", "ahat", "bhat", "chat"),
                                                                                           n.iter = nIterations,
                                                                                           n.chains = nChains,
                                                                                           n.burnin = nBurnin,
                                                                                           n.thin = nThin))


# Fit baseline model with SMC
baselineModel[[2]]  <- nimMCMCSMCupdates::spartaNimWeights(model = newModel,
                                                           latent = "x",
                                                           nParFiltRun = numParticles,
                                                           mcmc = TRUE,
                                                           pfType = pfTypeRun,
                                                           MCMCconfiguration = list(target = c('a', 'b', 'c'),
                                                                                    additionalPars = c("x", "ahat", "bhat", "chat"),
                                                                                    n.iter = nIterations,
                                                                                    n.chains = nChains,
                                                                                    n.burnin = nBurnin,
                                                                                    n.thin = nThin)
)
