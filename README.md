# nimMCMCSMCupdates
[![DOI](https://zenodo.org/badge/616969908.svg)](https://zenodo.org/doi/10.5281/zenodo.10495691)

## Introduction

This package comes with a sequential Monte Carlo algorithm for data
assimilation problems in ecology. This package is edited from the
*nimbleSMC* package Michaud et al. (2021), so the users of this package
need to be familiar with the NIMBLE software through its R-package
*nimble* (de Valpine et al. 2022) and then its R-package for sequential
Monte Carlo problems *nimbleSMC* (the reader is referred to Chapter 8 of
de Valpine et al. (2022) and Michaud et al. (2021) for details on how to
fit SSMs using SMC approach in NIMBLE).

## Installation

You can install the development version of *nimMCMCSMCupdates* from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("Peprah94/nimMCMCSMCupdates")
```

## Functions in the R-package

*nimMCMCSMCupdates* includes a set of functions that are used to update
posterior distribution of latent states and model parameters using the
bootstap and auxiliary particle filters. Refer to the [main
paper](https://peprah94.github.io/#Publications) and vignette on
[method](https://github.com/Peprah94/nimMCMCSMCupdates/blob/main/vignette/methodUnderstanding.pdf)
for further details on the methodology.

| Function name                  | Function description                                                                                 |
|--------------------------------|------------------------------------------------------------------------------------------------------|
| `buildAuxiliaryFilterUpdate()` | Create an updated auxiliary particle filter algorithm to estimate log-likelihood.                    |
| `buildBootstrapFilterUpdate()` | Create an updated bootstrap particle filter algorithm to estimate log-likelihood.                    |
| `sampler_RW_PF_blockUpdate()`  | The particle filter block sampler to perform particle MCMC.                                          |
| `spartaNimWeights()`           | Fit reduced model using MCMC.                                                                        |
| `spartaNimUpdates()`           | Fit updated model using the smc algorithm in [main paper](https://peprah94.github.io/#Publications). |
| `updateUtils()`                | Create a modelValues object for the MCMC samples from the reduced model.                             |

## Updating models with the R-package

We load the package needed to run the script. This package loads
*nimble* and *nimbleSMC* as dependencies.

``` r
library(nimMCMCSMCupdates)
```

We simulate a dynamic occupancy model, described in simulation study two
in the [main article](https://peprah94.github.io/#Publications), with
further detains in Kéry and Royle (2020).

``` r
# set the configurations for fitting the model
nyears = 30
nsites = 300
nvisits = 3
iNodePrevAll = c(29, 25) # the time steps used to fit the reduced models

# Set the configurations for MCMC
nIterations = 5000
nBurnin = 2000
nChains = 2
nThin = 15
numParticles = 10

dynOccupancyModels <- function(nyears, #number of years
                               nsites, #number of sites
                               nvisits, #number of visits
                               #intercept and covariate effect for initial occupancy probability
                               fixedPars = list(alphaPhi = - 2,
                                                betaPhi = 1.5),
                               # standard deviation hyperparameters
                               hyperParsSig = list(alphaP = 2,
                                                   betaP = 3,
                                                   alphaPsi = 2,
                                                   betaPsi = 3)){

  set.seed(1994)

  # Simulating covariates
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
  alphaP <- rnorm(1, mean = 4, sd = hyperParsSig$alphaP)
  betaP <- rnorm(1, mean = -2, sd = hyperParsSig$betaP)
  alphaPsi <- rnorm(1, mean = 3, sd = hyperParsSig$alphaPsi)
  betaPsi <- rnorm(1, mean = -2.5, sd = hyperParsSig$betaPsi)

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
                                  betaPsi = betaPsi)

  retList$occSites = psi.fs

  return(retList)
}

# Simulating and saving data
simData <- dynOccupancyModels(nyears = 20,
                               nsites = 30,
                               nvisits = 3)

 save(simData, file = "data/simDataDynamicOccupancy.RData")
```

The simulated occupancy data can be loaded from the R-package for this
example.

``` r
#data("simDataDynamicOccupancy")
str(simData)
```

    List of 7
     $ y               : int [1:30, 1:3, 1:20] 0 0 0 0 0 0 0 0 1 0 ...
     $ z               : int [1:30, 1:20] 0 0 0 0 0 0 0 0 1 0 ...
     $ covariates      :List of 4
      ..$ windSpeed          : num [1:30, 1:3, 1:20] -0.802 0.566 0.227 0.914 0.911 ...
      ..$ elevation          : num [1:30] 0.477 0.596 -0.669 -0.792 -0.591 ...
      ..$ springPrecipitation: num [1:30, 1:20] -0.4735 0.1146 -0.7951 -0.3556 0.0118 ...
      ..$ sizeOfBeak         : num [1:30, 1:20] 0.386 -0.623 0.424 -0.539 0.307 ...
     $ trueSigma       :List of 4
      ..$ alphaP  : num 2
      ..$ betaP   : num 3
      ..$ alphaPsi: num 2
      ..$ betaPsi : num 3
     $ truePars        :List of 2
      ..$ alphaPhi: num -2
      ..$ betaPhi : num 1.5
     $ covariateEffects:List of 4
      ..$ alphaP  : num 5.68
      ..$ betaP   : num -3.76
      ..$ alphaPsi: num 3.59
      ..$ betaPsi : num -1.27
     $ occSites        : num [1:20] 0.167 0.167 0.2 0.233 0.267 ...

Next, we define the nimble model for the dynamic occupancy model.

``` r
  dynOccupancyCode <- nimbleCode({

    # Prior distributions of hyperparameters
    alphaPSig ~ T(dnorm(0, 0.1), 0.001, )
    betaPSig ~ T(dnorm(0, 0.1), 0.001, )
    alphaPsiSig ~ T(dnorm(0, 0.1), 0.001, )
    betaPsiSig ~ T(dnorm(0, 0.1), 0.001, )

    # Prior distribution for intercept and covariate effect
    alphaPhi ~ dnorm(0, sd = 10)
    betaPhi ~ dnorm(0, sd = 10)

    # Prior distributions
    alphaP ~ dnorm(mean = 4, sd = alphaPSig)
    betaP ~ dnorm(mean = -2, sd =  betaPSig)
    alphaPsi ~ dnorm(mean = 3, sd =  alphaPsiSig)
    betaPsi ~ dnorm(mean = -2.5, sd =  betaPsiSig)

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

    colonisationProb ~ dunif(0.01, 1)

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
```

We assume reduced model has already been fittef with MCMC, and we aim to
use that result. We show how to use results from JAGS (Depaoli, Clifton,
and Cobb 2016) as the reduced model in the this algorithm
[here](https://github.com/Peprah94/nimMCMCSMCupdates/blob/main/vignette/fittingMethodUsingJAGSsample.pdf).
Although we have not tried it here, we think the same can be applied to
results from STAN (Carpenter et al. 2017).

With this example, we fit the reduced model with the `spartaNimWeights`
function.

``` r
# set the configurations for fitting the model
nyears = 20
nsites = 30
nvisits = 3
iNodePrev <- 18 # the time steps used to fit the reduced models

# Set the configurations for MCMC
nIterations = 5000
nBurnin = 2000
nChains = 2
nThin = 5
numParticles = 10
## define data, constants, and initial values
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

#   # NIMBLE model for reduced model
  newModelReduced <- nimbleModel(dynOccupancyCode,
                                 data = data,
                                 constants = constants,
                                 inits = inits,
                                 check = FALSE)

  # Fitting the model
  mcmcResults <- spartaNimWeights(model = newModelReduced,
                                               latent = "z",
                                               nParFiltRun = numParticles,
                                               mcmc = TRUE,
                                               pfType = pfTypeRun,
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
```

    ===== Monitors =====
    thin = 1: alphaP, alphaPhi, alphaPsi, alphaPSig, alphaPsiSig, betaP, betaPhi, betaPsi, betaPSig, betaPsiSig, colonisationProb, psi.fs, z
    ===== Samplers =====
    RW sampler (11)
      - alphaPSig
      - betaPSig
      - alphaPsiSig
      - betaPsiSig
      - alphaPhi
      - betaPhi
      - colonisationProb
      - alphaP
      - betaP
      - alphaPsi
      - betaPsi
    binary sampler (367)
      - z[]  (367 elements)
    |-------------|-------------|-------------|-------------|
    |-------------------------------------------------------|
    |-------------|-------------|-------------|-------------|
    |-------------------------------------------------------|

Note that the `spartaNimWeights` function can also be used to fit the
full model by settung $t = nyears$. or $iNodePrev = nyears$.

Now that we have results from the reduced model, we can then proceed to
fit the updated model using importance sampling.

``` r
numParticles = 10
## define data, constants, and initial values
  z <- apply(simData$y, c(1, 3), function(x){
    ifelse(sum(x) > 0, 1, NA)
  })

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

  # NIMBLE model for updated model
  newModelUpdated <- nimbleModel(dynOccupancyCode,
                                 data = data,
                                 constants = constants,
                                 inits = inits,
                                 check = FALSE)

#   # Fitting the model
  updatedModelResults <- nimMCMCSMCupdates::spartaNimUpdates(model = newModelUpdated, #nimble model
                                               reducedModel = newModelReduced,
                                               latent = "z",
                                               nParFiltRun = numParticles,
                                               pfType = "bootstrap", #change to auxiliary for auxiliary PF
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
                                               postReducedMCMC = mcmcResults,# MCMC summary to use as initial values
                                               pfControl = list(saveAll = TRUE,
                                                                timeIndex = 2,
                                                                smoothing = TRUE,
                                                                mcmc = TRUE,
                                                                iNodePrev = iNodePrev),
                                               nCores = 1
  )
```

    NULL
    ===== Monitors =====
    thin = 1: alphaP, alphaPhi, alphaPsi, alphaPSig, alphaPsiSig, betaP, betaPhi, betaPsi, betaPSig, betaPsiSig, colonisationProb, psi.fs, z
    ===== Samplers =====
    (no samplers assigned)
    ===== Comments =====
    thin = 1: alphaP, alphaPhi, alphaPsi, alphaPSig, alphaPsiSig, betaP, betaPhi, betaPsi, betaPSig, betaPsiSig, colonisationProb, psi.fs, z
    [1] "alphaPSig"   "betaPSig"    "alphaPsiSig" "betaPsiSig" 
    [1] "alphaP"           "betaP"            "alphaPsi"         "betaPsi"         
    [5] "alphaPhi"         "betaPhi"          "colonisationProb"
    [1] TRUE
    [1] "z"
    |-------------|-------------|-------------|-------------|
    |-------------------------------------------------------|
    NULL
    ===== Monitors =====
    thin = 1: alphaP, alphaPhi, alphaPsi, alphaPSig, alphaPsiSig, betaP, betaPhi, betaPsi, betaPSig, betaPsiSig, colonisationProb, psi.fs, z
    ===== Samplers =====
    (no samplers assigned)
    ===== Comments =====
    thin = 1: alphaP, alphaPhi, alphaPsi, alphaPSig, alphaPsiSig, betaP, betaPhi, betaPsi, betaPSig, betaPsiSig, colonisationProb, psi.fs, z
    [1] "alphaPSig"   "betaPSig"    "alphaPsiSig" "betaPsiSig" 
    [1] "alphaP"           "betaP"            "alphaPsi"         "betaPsi"         
    [5] "alphaPhi"         "betaPhi"          "colonisationProb"
    [1] TRUE
    [1] "z"
    |-------------|-------------|-------------|-------------|
    |-------------------------------------------------------|

The samples from the algorithm are of class `mcmc.list`. Therefore,
packages like *ggmcmc* (Fernández-i-Marín 2016) and *coda* (Plummer et
al. 2006) can be used to do the posterior analysis.

``` r
#data("reducedModelResults")
#data("updatedModelResults")
class(mcmcResults$samples)
```

    [1] "mcmc.list"

``` r
class(updatedModelResults$samples)
```

    [1] "mcmc.list"

The results also contain the summary of the parameters monitored. Here,
we present the summary of the updated model.

``` r
head(updatedModelResults$summary$all.chains)
```

                     Mean    Median   St.Dev.   95%CI_low  95%CI_upp
    alphaP       4.384379  4.342556 0.4350456  3.68843274  5.3769865
    alphaPSig    1.587326  1.158108 1.3808694  0.08402534  5.2705133
    alphaPhi    -2.311743 -2.201483 0.8924351 -4.31790956 -0.8886312
    alphaPsi     5.022900  4.931248 0.9029854  3.42794720  6.9861129
    alphaPsiSig  2.791914  2.428226 1.5747295  0.68731967  6.5483746
    betaP       -2.067356 -2.039567 0.5640460 -3.35449672 -1.0032202

The function `spartaNimUpdates` has been been set up to run in parallel
between chains. The output contains the time taken for each chain to
run.

``` r
#Time taken
updatedModelResults$timeRun
```

    $chain1
    Time difference of 1.518683 secs

    $chain2
    Time difference of 1.543938 secs

    $all.chains
    Time difference of 3.062621 secs

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-carpenter2017stan" class="csl-entry">

Carpenter, Bob, Andrew Gelman, Matthew D Hoffman, Daniel Lee, Ben
Goodrich, Michael Betancourt, Marcus Brubaker, Jiqiang Guo, Peter Li,
and Allen Riddell. 2017. “Stan: A Probabilistic Programming Language.”
*Journal of Statistical Software* 76 (1).

</div>

<div id="ref-nimblepackage" class="csl-entry">

de Valpine, Perry, Christopher Paciorek, Daniel Turek, Nick Michaud,
Cliff Anderson-Bergman, Fritz Obermeyer, Claudia Wehrhahn Cortes, Abel
Rodrìguez, Duncan Temple Lang, and Sally Paganin. 2022. *NIMBLE User
Manual* (version 0.13.1). <https://doi.org/10.5281/zenodo.1211190>.

</div>

<div id="ref-depaoli2016just" class="csl-entry">

Depaoli, Sarah, James P Clifton, and Patrice R Cobb. 2016. “Just Another
Gibbs Sampler (JAGS) Flexible Software for MCMC Implementation.”
*Journal of Educational and Behavioral Statistics* 41 (6): 628–49.

</div>

<div id="ref-ggmcmc" class="csl-entry">

Fernández-i-Marín, Xavier. 2016. “<span class="nocase">ggmcmc</span>:
Analysis of MCMC Samples and Bayesian Inference.” *Journal of
Statistical Software* 70 (9): 1–20.
<https://doi.org/10.18637/jss.v070.i09>.

</div>

<div id="ref-kery2020applied" class="csl-entry">

Kéry, Marc, and J Andrew Royle. 2020. *Applied Hierarchical Modeling in
Ecology: Analysis of Distribution, Abundance and Species Richness in r
and BUGS: Volume 2: Dynamic and Advanced Models*. Academic Press.

</div>

<div id="ref-michaud2021sequential" class="csl-entry">

Michaud, Nicholas, Perry de Valpine, Daniel Turek, Christopher J
Paciorek, and Dao Nguyen. 2021. “Sequential Monte Carlo Methods in the
Nimble and nimbleSMC r Packages.” *Journal of Statistical Software* 100:
1–39.

</div>

<div id="ref-coda" class="csl-entry">

Plummer, Martyn, Nicky Best, Kate Cowles, and Karen Vines. 2006. “CODA:
Convergence Diagnosis and Output Analysis for MCMC.” *R News* 6 (1):
7–11. <https://journal.r-project.org/archive/>.

</div>

</div>
