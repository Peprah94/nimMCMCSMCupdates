# nimMCMCSMCupdates

# Introduction

This package comes with a sequential Monte Carlo algorithm for data
assimilation problems in ecology. This package is editted from the
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

``` r
library(dplyr)
library(nimMCMCSMCupdates)
```

we load the simulated occupancy data in the package for this example.

``` r
data("simDataDynamicOccupancy")
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

The nimble code

``` r
iNodePrev <- 18
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
use that result. In the case of this simulation, we can fit the reduced
model with the `spartaNimWeights` function.

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

``` r
   #save results
  save(mcmcResults, file= paste0("data/reducedModelResults.RData")) 
```

## Other Materials

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-nimblepackage" class="csl-entry">

de Valpine, Perry, Christopher Paciorek, Daniel Turek, Nick Michaud,
Cliff Anderson-Bergman, Fritz Obermeyer, Claudia Wehrhahn Cortes, Abel
Rodrìguez, Duncan Temple Lang, and Sally Paganin. 2022. *NIMBLE User
Manual* (version 0.13.1). <https://doi.org/10.5281/zenodo.1211190>.

</div>

<div id="ref-michaud2021sequential" class="csl-entry">

Michaud, Nicholas, Perry de Valpine, Daniel Turek, Christopher J
Paciorek, and Dao Nguyen. 2021. “Sequential Monte Carlo Methods in the
Nimble and nimbleSMC r Packages.” *Journal of Statistical Software* 100:
1–39.

</div>

</div>
