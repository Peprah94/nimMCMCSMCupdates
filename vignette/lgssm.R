# This is the general function that is used to fit the model with MCMC
# and simulate the required data


# Load the required packages
library('nimble')
library(nimbleSMC)
library(nimMCMCSMCupdates)

nIterations = 50000
nChains = 2
nBurnin = 20000
nThin = 1

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

# simulate data
simData <- sim2(0.5, 1,1,1,50,1)

# NIMBLE CODE
stateSpaceCode <- nimbleCode({
  x[1] ~ dnorm(1, sd = 1)
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
  mu0= 0.8,
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
                                                                             MCMCconfiguration = list(target = c('a', 'b',  'c'),
                                                                                                      additionalPars = c("x", "ahat", "bhat","chat"),
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
  t = 50,
  aTrue = 0.5,
  bTrue = 1,
  cTrue = 1,
  muTrue = 1#,
  # mu = example1ReducedModelTrue[[iNodetag]]$summary$all.chains["mu", 1]
)


inits <- list(
  a = example1ReducedModelTrue$summary$all.chains["a",1 ],
  b = example1ReducedModelTrue$summary$all.chains["b", 1],
  # x= c(example1ReducedModelTrue[[iNodetag]]$summary$all.chains[grepl("x", rownames(example1ReducedModelTrue[[iNodetag]]$summary$all.chains)), 1], rep(0, (constants$t - (iNodePrev[iNodetag])))),
  c = example1ReducedModelTrue$summary$all.chains["c",1 ]
)


# Define model for updated model
newModelUpdated <- nimbleModel(stateSpaceCode,
                               data = data,
                               constants = constants,
                               inits = inits,
                               check = FALSE)


# Fit model with results from fitted reduced MCMC model fit
example1UpdatedModelTrueBF <- nimMCMCSMCupdates::spartaNimUpdates(model = newModelUpdated, #nimble model
                                                                               reducedModel = newModelReduced,
                                                                               latent = "x", #latent variable
                                                                               nParFiltRun = 100,
                                                                               mcmcScale = 1,
                                                                               extraVars = NULL,
                                                                               adaptiveSampler = TRUE,
                                                                               MCMCconfiguration = list(target = c('a', 'c', 'b'),
                                                                                                        additionalPars = c("x", "ahat", "chat", "bhat"),
                                                                                                        n.iter = (nIterations - nBurnin)/2,
                                                                                                        n.chains = nChains,
                                                                                                        n.burnin = 500,
                                                                                                        n.thin = 1),  #saved loglikelihoods from reduced model
                                                                               postReducedMCMC = example1ReducedModelTrue,# MCMC summary to use as initial values
                                                                               pfControl = list(saveAll = TRUE,
                                                                                                smoothing = TRUE,
                                                                                                mcmc = TRUE,
                                                                                                M = 20,
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














runFunction <- function(iter,
                        simDataAll,
                        numParticles,
                        mcmcRun = TRUE,
                        pfTypeRun,
                        nIterations,
                        nBurnin,
                        nChains,
                        nThin,
                        nyears,
                        iNodePrev ){


  library('nimble')
  library(nimbleSMC)
  library(nimMCMCSMCupdates)


  # create empty lists to store the results
  example1ReducedModelTrue <- list()
  #example1ReducedModelFalse <- list()
  example1UpdatedModelTrueBF <- list()
  example1UpdatedModelTrueAux <- list()
  #example1UpdatedModelFalse <- list()
  baselineModel <- list()




  #load data
  simData <- simDataAll[[iter]]


  # NIMBLE CODE
  stateSpaceCode <- nimbleCode({
    x[1] ~ dnorm(1, sd = 1)
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
  })


  # ## define data
  data <- list(
    y = simData$y
  )


  # define constants
  constants <- list(
    t = nyears,
    aTrue = 0.8,
    bTrue = 1,
    cTrue = 1,
    muTrue = 1
  )


  #define initial values
  inits <- list(
    a = 0.2,
    b = 0.8,
    mu0= 0.8,
    c = 0.5
  )


  # ## build the model
  stateSpaceModel <- nimbleModel(stateSpaceCode,
                                 data = data,
                                 constants = constants,
                                 inits = inits,
                                 check = FALSE)


  for(iNodetag in seq_along(iNodePrev)){


    data <- list(
      y = simData$y[-c((iNodePrev[iNodetag]+1):50)]
    )


    constants <- list(
      t = iNodePrev[iNodetag],
      aTrue = 0.8,
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
    example1ReducedModelTrue[[iNodetag]]  <- nimMCMCSMCupdates::spartaNimWeights(model = newModelReduced,
                                                                                 latent = "x",
                                                                                 nParFiltRun = numParticles,
                                                                                 mcmc = TRUE,
                                                                                 block = FALSE,
                                                                                 pfType = pfTypeRun,
                                                                                 MCMCconfiguration = list(target = c('a', 'b',  'c'),
                                                                                                          additionalPars = c("x", "ahat", "bhat","chat"),
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
      t = nyears,
      aTrue = 0.8,
      bTrue = 1,
      cTrue = 1,
      muTrue = 1#,
      # mu = example1ReducedModelTrue[[iNodetag]]$summary$all.chains["mu", 1]
    )


    inits <- list(
      a = example1ReducedModelTrue[[iNodetag]]$summary$all.chains["a",1 ],
      b = example1ReducedModelTrue[[iNodetag]]$summary$all.chains["b", 1],
      # x= c(example1ReducedModelTrue[[iNodetag]]$summary$all.chains[grepl("x", rownames(example1ReducedModelTrue[[iNodetag]]$summary$all.chains)), 1], rep(0, (constants$t - (iNodePrev[iNodetag])))),
      c = example1ReducedModelTrue[[iNodetag]]$summary$all.chains["c",1 ]
    )


    # Define model for updated model
    newModelUpdated <- nimbleModel(stateSpaceCode,
                                   data = data,
                                   constants = constants,
                                   inits = inits,
                                   check = FALSE)


    # Fit model with results from fitted reduced MCMC model fit
    example1UpdatedModelTrueBF[[iNodetag]]  <- nimMCMCSMCupdates::spartaNimUpdates(model = newModelUpdated, #nimble model
                                                                                   reducedModel = newModelReduced,
                                                                                   latent = "x", #latent variable
                                                                                   nParFiltRun = numParticles,
                                                                                   mcmcScale = 0.1,
                                                                                   extraVars = NULL,
                                                                                   adaptiveSampler = FALSE,
                                                                                   MCMCconfiguration = list(target = c('a', 'c', 'b'),
                                                                                                            additionalPars = c("x", "ahat", "chat", "bhat"),
                                                                                                            n.iter = (nIterations - nBurnin)/nThin,
                                                                                                            n.chains = nChains,
                                                                                                            n.burnin = 1000,
                                                                                                            n.thin = 1),  #saved loglikelihoods from reduced model
                                                                                   postReducedMCMC = example1ReducedModelTrue[[iNodetag]],# MCMC summary to use as initial values
                                                                                   pfControl = list(saveAll = TRUE,
                                                                                                    smoothing = TRUE,
                                                                                                    mcmc = TRUE,
                                                                                                    M = nyears - iNodePrev[iNodetag],
                                                                                                    iNodePrev = iNodePrev[iNodetag]),
                                                                                   nCores = 1
    )




    example1UpdatedModelTrueAux[[iNodetag]]  <- nimMCMCSMCupdates::spartaNimUpdates(model = newModelUpdated, #nimble model
                                                                                    reducedModel = newModelReduced,
                                                                                    latent = "x", #latent variable
                                                                                    nParFiltRun = numParticles,
                                                                                    propCov = diag(2)* c(0.1,1),
                                                                                    pfType = "auxiliary",
                                                                                    extraVars = NULL,
                                                                                    MCMCconfiguration = list(target = c('a', 'c'),
                                                                                                             additionalPars = c("x", "ahat", "chat"),
                                                                                                             n.iter = (nIterations - nBurnin)/nThin,
                                                                                                             n.chains = nChains,
                                                                                                             n.burnin = 0,
                                                                                                             n.thin = 1),  #saved loglikelihoods from reduced model
                                                                                    postReducedMCMC = example1ReducedModelTrue[[iNodetag]],# MCMC summary to use as initial values
                                                                                    pfControl = list(saveAll = TRUE,
                                                                                                     smoothing = TRUE,
                                                                                                     mcmc = TRUE,
                                                                                                     lookahead = "mean",
                                                                                                     M = nyears - iNodePrev[iNodetag],
                                                                                                     iNodePrev = iNodePrev[iNodetag])
    )


  }




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


  #return results
  ret <- c(baselineModel,
           example1ReducedModelTrue,
           #example1ReducedModelFalse,
           example1UpdatedModelTrueBF,
           example1UpdatedModelTrueAux
           # example1UpdatedModelFalse
  )


  return(ret)
}




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
