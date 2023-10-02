## load the nimble library and set seed
library('nimble')
library(nimbleSMC)
devtools::load_all(".")
set.seed(1)
mcmcRun <- FALSE #use mcmc or nimbleSMC for reduced Model
pfTypeRun = "auxiliary"
#source("spartaUpdatingFunctions.R")
#source("~/Documents/GitHub/myphdthesis/R/MCMCSamplersUpdate.R")
#load("case2SimData.RData")
#load("reducedCase2.RData")

# Setting up MCMC configuration values and variation of parameters
nIterations = 5000
nBurnin = 2000
nChains = 2
nThin = 2
nyears = 50
aVars <- c(0.1, 0.5) # changing the intercept
#High and small values of a
iNodePrev <- 49  # The number of years for reduced model

aVarstag = 2
#for(aVarstag in 1:length(aVars)){

#Function to simulate data
sim2 <- function(a, b, c, t, mu0){
  x <- y <- numeric(t)
  x[1] <- rnorm(1, a, 1 )
  y[1] <- rnorm(1, x[1], 1)

  for(k in 2:t){
    x[k] <- rnorm(1, a*x[k -1] , 1)
    y[k] <- rnorm(1, x[k]*c, 1)# + (sigOE * (sqrt(df -2)/df) * rt(1, df))
  }
  return(list(x=x, y=y))
}

message("simulating data for a = ", aVars[aVarstag])
simData <- sim2(a = aVars[aVarstag],
                b = 1,
                c = 1,
                t = nyears,
                mu0 = 1)

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
#
#
# ## build the model
stateSpaceModel <- nimbleModel(stateSpaceCode,
                               data = data,
                               constants = constants,
                               inits = inits,
                               check = FALSE)
#


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
  message(paste("Running updated model for iNodePrev = ", iNodePrev, "and a = ", aVars[aVarstag]))
  message(paste("Running reduced model for iNodePrev = ", iNodePrev, "and a = ", aVars[aVarstag]))
  message(paste("Reduced model configuration for iNodePrev = ", iNodePrev, "and a = ", aVars[aVarstag]))
  data <- list(
    y = simData$y#[c((iNodePrev[iNodetag]):50)]
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
                                                              nParFiltRun = 100,
                                                              propCov = diag(2)* c(0.1,0.5),
                                                              pfType = "auxiliary",
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
   ggmcmc::ggs_traceplot(., family = c("a"))

 ggmcmc::ggs(example1ReducedModel$samples)%>%
   ggmcmc::ggs_traceplot(., family = c("a"))

  #save(example1UpdatedModel, file= paste0("example1UpdatedIn",iNodePrev[iNodetag],"A",aVarstag,".RData"))

#}

#}
