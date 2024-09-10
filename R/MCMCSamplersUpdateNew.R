#' This sampler is modified from the algorithms in nimbleSMC. We only provide algorithms for
#' Block sampling. Other samplers will be explored in the future.

#' Particle Filtering MCMC Sampling Algorithms
#'
#' Details of the particle filtering MCMC sampling algorithms provided in nimbleSMC.
#'
#' @param model (uncompiled) model on which the MCMC is to be run
#' @param mvSaved \code{modelValues} object to be used to store MCMC samples
#' @param target node(s) on which the sampler will be used
#' @param control named list that controls the precise behavior of the sampler, with elements specific to \code{samplertype}.  The default values for control list are specified in the setup code of each sampling algorithm.  Descriptions of each sampling algorithm, and the possible customizations for each sampler (using the \code{control} argument) appear below.
#' @section RW_PF_block_Update sampler:
#'
#' The particle filter block sampler allows the user to perform particle MCMC (PMCMC) (Andrieu et al., 2010) for multiple parameters jointly, primarily for state-space or hidden Markov models of time-series data.  This method uses Metropolis-Hastings block samplers for top-level parameters but uses the likelihood approximation of a particle filter (sequential Monte Carlo) to integrate over latent nodes in the time-series.  The \code{RW_PF} sampler uses an adaptive Metropolis-Hastings algorithm with a multivariate normal proposal distribution.  Note that samples of the latent states can be retained as well, but the top-level parameter being sampled must be a scalar.   A bootstrap, auxiliary, or user defined particle filter can be used to integrate over latent states.
#'
#' For more information about user-defined samplers within a PMCMC sampler, see the NIMBLE User Manual.
#'
#' The \code{RW_PF_block} sampler accepts the following control list elements:
#' \itemize{
#' \item adaptive. A logical argument, specifying whether the sampler should adapt the proposal covariance throughout the course of MCMC execution. (default = TRUE)
#' \item adaptScaleOnly. A logical argument, specifying whether adaptation should be done only for \code{scale} (TRUE) or also for \code{provCov} (FALSE).  This argument is only relevant when \code{adaptive = TRUE}.  When \code{adaptScaleOnly = FALSE}, both \code{scale} and \code{propCov} undergo adaptation; the sampler tunes the scaling to achieve a theoretically good acceptance rate, and the proposal covariance to mimic that of the empirical samples.  When \code{adaptScaleOnly = TRUE}, only the proposal scale is adapted. (default = FALSE)
#' \item adaptInterval. The interval on which to perform adaptation. (default = 200)
#' \item scale. The initial value of the scalar multiplier for \code{propCov}.  If \code{adaptive = FALSE}, \code{scale} will never change. (default = 1)
#' \item adaptFactorExponent. Exponent controling the rate of decay of the scale adaptation factor.  See Shaby and Wells, 2011, for details. (default = 0.8)
#' \item propCov. The initial covariance matrix for the multivariate normal proposal distribution.  This element may be equal to the \code{'identity'}, in which case the identity matrix of the appropriate dimension will be used for the initial proposal covariance matrix. (default is \code{'identity'})
#' \item pfNparticles.  The number of particles to use in the approximation to the log likelihood of the data (default = 1000).
#' \item latents.  Character vector specifying the nodes that are latent states over which the particle filter will operate to approximate the log-likelihood function.
#' \item pfType.  Character argument specifying the type of particle filter that should be used for likelihood approximation.  Choose from \code{"bootstrap"} and \code{"auxiliary"}.  Defaults to \code{"bootstrap"}.
#' \item pfControl.  A control list that is passed to the particle filter function.  For details on control lists for bootstrap or auxiliary particle filters, see \code{\link{buildBootstrapFilter}} or \code{\link{buildAuxiliaryFilter}} respectively.  Additionally, this can be used to pass custom arguments into a user defined particle filter.
#' \item pfOptimizeNparticles.  A logical argument, specifying whether to automatically determine the optimal number of particles to use, based on Pitt and Shephard (2011).  This will override any value of \code{pfNparticles} specified above.
#' \item pf.  A user-defined particle filter object, if a bootstrap or auxiliary particle filter is not adequate.
#' }
#'
#' @name SMCsamplers
#' @aliases samplers sampler RW_PF RW_PF_block sampler_RW_PF sampler_RW_PF_block
#'
#'

#######################################################################################
### RW_PF_block, does a block RW, but using a particle filter likelihood function #####
#######################################################################################
#' @rdname samplers
#' @export
sampler_RW_PF_blockUpdate <- nimbleFunction(
  name = 'sampler_RW_PF_blockUpdate',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## control list extraction
    adaptive            <- extractControlElement(control, 'adaptive',             FALSE)
    adaptScaleOnly      <- extractControlElement(control, 'adaptScaleOnly',       FALSE)
    adaptInterval       <- extractControlElement(control, 'adaptInterval',        200)
    adaptFactorExponent <- extractControlElement(control, 'adaptFactorExponent',  0.8)
    scale               <- extractControlElement(control, 'scale',                1)
    propCov             <- extractControlElement(control, 'propCov',              'identity')
    propCov1             <- extractControlElement(control, 'propCov1',              'identity')
    existingPF          <- extractControlElement(control, 'pf',                   NULL)
    m                   <- extractControlElement(control, 'pfNparticles',         1000)
    filterType          <- extractControlElement(control, 'pfType',               'bootstrapUpdate')
    filterControl       <- extractControlElement(control, 'pfControl',            list())
    optimizeM           <- extractControlElement(control, 'pfOptimizeNparticles', FALSE)
    maxDimCovHistory    <- extractControlElement(control, 'maxDimCovHistory',    10)
    latents             <- extractControlElement(control, 'latents',              error = 'RW_PF sampler missing required control argument: latents')
    # postSamples <- extractControlElement(control, 'postSamples', double())
    extraVars <- extractControlElement(control, 'extraVars', double())
    mvSamplesEst <- extractControlElement(control, 'mvSamplesEst', double())
    reducedModel <- extractControlElement(control, 'reducedModel', double())
    iNodePrev <- extractControlElement(control, 'iNodePrev', 1)
    additionalPars <- extractControlElement(control, 'additionalPars', double())
    #target <- extractControlElement(control, 'target', double())

    if('pfLookahead' %in% names(control)) {
      warning("The `pfLookahead` control list argument is deprecated and will not be supported in future versions of `nimbleSMC`. Please specify the lookahead function via the pfControl argument instead.")
      filterControl$lookahead <- control[['pfLookahead']]
    }
    else if(!('lookahead' %in% names(filterControl))) {
      filterControl$lookahead <- 'simulate'
    }

    ##############################################################
    # Obtaining latent states, hyperparameters (topParams)
    # and intermediate parameters and the
    # their dependencies
    #############################################################

    ## node list generation
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)

    extraVars <- model$expandNodeNames(extraVars)
   # print(extraVars)
    if(length(extraVars) > 0){
      targetAsScalar <- targetAsScalar[!targetAsScalar %in% extraVars]
    }

    # Get dependencies of the target variables
    calcNodes <- model$getDependencies(target)

    # We always sample the latent states with this application
    # and get their dependencies
    latentSamp <- TRUE
    latentDep <- model$getDependencies(latents)
    latentAsScalar <- model$expandNodeNames(latents,
                                            returnScalarComponents = TRUE)

    ## Now we try to seperate the hyperparameters from the other model parameters
    ## We do this by looking at the model graph and going one step back at a time.
    ## Currently, we want to update model parameter that are greater than 5
    ## This is not an efficient approach, but it works for the application in the manuscript.

    # Get top parameters
    #topParams <- model$getNodeNames(stochOnly=TRUE, includeData=FALSE, topOnly=TRUE)

    targetSamplesRun <- lapply(model$getVarNames(nodes = targetAsScalar), function(x){
      ret <- model$getDependencies(x, stochOnly = TRUE,  downstream = TRUE, self = FALSE)
      nodesToReturn <- model$getVarNames(nodes = ret)
      return(nodesToReturn)
    })

    ##########################################################
    # Last run: Those that depend on the observation data
    #########################################################
    lastRun <- sapply(targetSamplesRun, function(x){
      ret <- all(model$isData(x))
      return(ret)
    })

    lastRunPars <- model$expandNodeNames(model$getVarNames(nodes = target)[lastRun])
    if(length(extraVars) > 0){
      lastRunPars <- lastRunPars[!lastRunPars %in% extraVars]
    }
    print(paste("Parameters that affect observed data:", c(lastRunPars)))

    #####################################
    # Third run: Hyperparameters of the last run parameters
    # I will call it third run
    #####################################
    thirdRun <- sapply(targetSamplesRun, function(x){
      any(x %in% model$getVarNames(nodes = lastRunPars))
    })
    thirdRunPars <-model$expandNodeNames(model$getVarNames(nodes = target)[thirdRun])
    if(length(extraVars) > 0){
      thirdRunPars <- thirdRunPars[!thirdRunPars %in% extraVars]
    }
    print(paste("Hyperparameters of parameters that affect observed data:", c(thirdRunPars)))

    ############################################
    # Second run: The parameters that directly affect the latent state
    ####################################
    suppressWarnings(secondRun <- sapply(targetSamplesRun, function(x){
      ret <- x == c(model$getVarNames(nodes = latents), model$getVarNames(nodes = model$getDependencies(lastRunPars, self = FALSE, stochOnly = TRUE)))
    ret <- all(ret == TRUE)
    return(ret)
      }))
    secondRunPars <- model$expandNodeNames(model$getVarNames(nodes = target)[secondRun])
    if(length(extraVars) > 0){
      secondRunPars <- secondRunPars[!secondRunPars %in% extraVars]
    }
    print(paste("Parameters that affect latent state:", c(secondRunPars)))

    ############################################
    # First run: The parameters that affect second run
    ####################################
    firstRun <- sapply(targetSamplesRun, function(x){
      any(x %in% model$getVarNames(nodes = secondRunPars))
    })
    firstRunPars <- model$expandNodeNames(model$getVarNames(nodes = target)[firstRun])
    if(length(extraVars) > 0){
      firstRunPars <- firstRunPars[!firstRunPars %in% extraVars]
    }
    print(paste("Hyperparameters of parameters that affect latent state:", c(firstRunPars)))

    #####################################
    # Set the conditions for the MCMC run
    #####################################
    simLastPars <- TRUE
    simThirdPars <- TRUE

    if(length(lastRunPars) == 0) simLastPars <- FALSE
    if(length(thirdRunPars) == 0) simThirdPars <- FALSE


    if(length(firstRunPars) == 0){
      simFirstPars <- FALSE
    } else{
      simFirstPars <- TRUE
      topParams <- firstRunPars
    }

    if(length(secondRunPars) == 0){
      simSecondPars <- FALSE
      }else{
      simSecondPars <- TRUE
      topParamsInter <- secondRunPars
    }

    parsToSim <- c("Hyperparameters of latent state parameters",
               "Latent state parameeters",
               "Hyperparameters of observed data",
               "Observed data parameters")
   print(paste("Parameters sampled:",
               parsToSim[c(simFirstPars, simSecondPars, simThirdPars, simLastPars)]))

    # Extrace derived quantities we need to store
    extraScalarPars <- targetAsScalar[!targetAsScalar %in% model$expandNodeNames(c(firstRunPars, secondRunPars, thirdRunPars, lastRunPars))]


    #####################################################
    # Set samplers for first, third and last parameters to un
    ######################################################
    if(simFirstPars){
      isDiscrete <- any(model$isDiscrete(topParams) == TRUE)
      if(isDiscrete){
        if(length(topParams) > 1){
          my_topParamsSampler <- sampler_AF_slice(model,
                                                  mvSaved = mvSaved,
                                                  topParams,
                                                  control = list(maxSteps = 100))
        } else {
          my_topParamsSampler <- sampler_slice(model,
                                               mvSaved = mvSaved,
                                               topParams, control = list(maxSteps = 100,
                                                                         adaptive = adaptive))
        }
      } else {
        if(length(topParams) > 1){
          my_topParamsSampler <- sampler_RW_block(model,
                                                  mvSaved = mvSaved,
                                                  topParams,
                                                  control = list(scale = scale,
                                                                 adaptive = adaptive,
                                                                 propCov = "identity",
                                                                 adaptInterval = adaptInterval))
        } else {
          my_topParamsSampler <- sampler_RW(model,
                                            mvSaved = mvSaved,
                                            topParams,
                                            control = list(scale = scale,
                                                           adaptive = adaptive,
                                                           adaptInterval = adaptInterval,
                                                           adaptScaleOnly = adaptScaleOnly))
        }
      }
    } else{
      my_topParamsSampler <-  sampler_RW_block(model,
                                               mvSaved = mvSaved,
                                               topParamsInter,
                                               control = list(scale = scale,
                                                              adaptive = adaptive,
                                                              adaptInterval = adaptInterval,
                                                              adaptScaleOnly = adaptScaleOnly))
    }

    if(simThirdPars){
      isDiscrete <- any(model$isDiscrete(thirdRunPars) == TRUE)
      if(isDiscrete){
        if(length(thirdRunPars) > 1){
          my_thirdSimParamsSampler <- sampler_AF_slice(model,
                                                       mvSaved = mvSaved,
                                                       thirdRunPars,
                                                       control = list(maxSteps = 100))
        } else {
          my_thirdSimParamsSampler <- sampler_slice(model,
                                                    mvSaved = mvSaved,
                                                    thirdRunPars,
                                                    control = list(maxSteps = 100,
                                                                   adaptive = adaptive))
        }
      } else {
        if(length(thirdRunPars) > 1){
          my_thirdSimParamsSampler <- sampler_RW_block(model,
                                                       mvSaved = mvSaved,
                                                       thirdRunPars,
                                                       control = list(scale = scale,
                                                                      adaptive = adaptive,
                                                                      adaptInterval = adaptInterval,
                                                                      adaptScaleOnly = adaptScaleOnly))
        } else {
          my_thirdSimParamsSampler <- sampler_RW(model,
                                                 mvSaved = mvSaved,
                                                 thirdRunPars,
                                                 control = list(scale = scale,
                                                                adaptive = adaptive,
                                                                adaptInterval = adaptInterval,
                                                                adaptScaleOnly = adaptScaleOnly))
        }
      }
    } else {
      my_thirdSimParamsSampler <-  sampler_RW_block(model,
                                                    mvSaved = mvSaved,
                                                    topParamsInter,
                                                    control = list(scale = scale,
                                                                   adaptive = adaptive,
                                                                   adaptInterval = adaptInterval,
                                                                   adaptScaleOnly = adaptScaleOnly))
    }

    if(simLastPars){
      isDiscrete <- any(model$isDiscrete(lastRunPars) == TRUE)
      if(isDiscrete){
        if(length(lastRunPars) > 1){
          my_lastSimParamsSampler <- sampler_AF_slice(model,
                                                      mvSaved = mvSaved,
                                                      lastRunPars,
                                                      control = list(maxSteps = 100))
        } else {
          my_lastSimParamsSampler <- sampler_slice(model,
                                                   mvSaved = mvSaved,
                                                   lastRunPars,
                                                   control = list(maxSteps = 100,
                                                                  adaptive = adaptive))
        }
      } else {
        if(length(lastRunPars) > 1){
          my_lastSimParamsSampler <- sampler_RW_block(model,
                                                      mvSaved = mvSaved,
                                                      lastRunPars,
                                                      control = list(scale = scale,
                                                                     adaptive = adaptive,
                                                                     adaptInterval = adaptInterval,
                                                                     adaptScaleOnly = adaptScaleOnly))
        } else {
          my_lastSimParamsSampler <- sampler_RW(model,
                                                mvSaved = mvSaved,
                                                lastRunPars,
                                                control = list(scale = scale,
                                                               adaptive = adaptive))
        }
      }
    } else{
      my_lastSimParamsSampler <-  sampler_RW_block(model,
                                                   mvSaved = mvSaved,
                                                   topParamsInter,
                                                   control = list(scale = scale,
                                                                  adaptive = adaptive,
                                                                  adaptInterval = adaptInterval,
                                                                  adaptScaleOnly = adaptScaleOnly))
    }

    # If there are extraVars, we need to add them to the target as scalar
    # to correctly estimate the log-likelihood
    extraVars <- model$expandNodeNames(extraVars)
   # print(extraVars)
    if(length(secondRunPars) == 0){
      extras <- FALSE
    } else {
      extras <- TRUE
      targetAsScalar <- c(targetAsScalar,
                          extraVars)
    }
   # print(extras)
    #print(targetAsScalar)

    ################
    # What about if the topParamsInter are time dependent?
    ########################
    # I need to seperate the values that are in the reduced model
    # from those in the full model
    # I save the values from the reduced model at each iteration
    # but simulate values for those that are yet to be sorted out
    timeDependent <- FALSE

    #####################################################
    # Set parameterizations for the MCMC run
    ####################################################
    optimizeM     <- as.integer(optimizeM)
    scaleOriginal <- scale
    timesRan      <- 0
    timesAccepted <- 0
    timesAdapted  <- 0
    prevLL        <- 0
    nVarEsts      <- 0
    itCount       <- 0
    d <- length(topParamsInter)
    scaleHistory <- c(0, 0)                                                                  ## scaleHistory
    acceptanceHistory <- c(0, 0)                                                             ## scaleHistory
    propCovHistory <- if(d <= maxDimCovHistory) array(0, c(2,d,d)) else array(0, c(2,2,2))   ## scaleHistory
    saveMCMChistory <- if(getNimbleOption('MCMCsaveHistory')) TRUE else FALSE
    if(is.character(propCov) && propCov == 'identity')     propCov <- diag(d)
    propCovOriginal <- propCov
    chol_propCov <- chol(propCov)
    chol_propCov_scale <- scale * chol_propCov
    empirSamp <- matrix(0, nrow=adaptInterval, ncol=d)

    storeParticleLP <- -Inf
    iterRan <- 1
    storeLLVar  <- 0
    nVarReps <- 7    ## number of LL estimates to compute to get each LL variance estimate for m optimization
    mBurnIn  <- 15   ## number of LL variance estimates to compute before deciding optimal m
    if(optimizeM)   m <- 3000


    my_decideAndJump <- myDecideAndJump(model, mvSaved, topParamsInter)
    my_calcAdaptationFactor <- calcAdaptationFactor(d, adaptFactorExponent)


    ###########################################################
    #  Particle filter Configurations
    ###########################################################
    if(!is.null(existingPF)) {
      my_particleFilter <- existingPF
    } else {
      if(latentSamp == TRUE) {
        filterControl$saveAll <- TRUE
        filterControl$smoothing <- TRUE
      } else {
        filterControl$saveAll <- FALSE
        filterControl$smoothing <- FALSE
      }
      filterControl$initModel <- FALSE
      if(is.character(filterType) && filterType == 'auxiliary') {
        my_particleFilter <- buildAuxiliaryFilterNew(model,
                                                     latents,
                                                     control = filterControl)
      }
      else if(is.character(filterType) && filterType == 'auxiliaryUpdate' ) {
        my_particleFilter <- buildAuxiliaryFilterUpdate(model,
                                                        latents,
                                                        mvSamplesEst = mvSamplesEst,
                                                        target = target,
                                                        control = filterControl)

      }
      else if(is.character(filterType) && filterType == 'bootstrapUpdate' ) {
        my_particleFilter <- buildBootstrapFilterUpdate(model,
                                                        latents,
                                                        mvSamplesEst = mvSamplesEst,
                                                        target = target,
                                                        control = filterControl)


      }
      else if(is.character(filterType) && filterType == 'bootstrap') {
        my_particleFilter <- buildBootstrapFilterNew(model,
                                                     latents,
                                                     control = filterControl)

      }
      else if(is.nfGenerator(filterType)){
        my_particleFilter <- filterType(model, latents, control = filterControl)
      }
      else stop('filter type must be either "bootstrap", "auxiliary", or a
                  user defined filtering algorithm created by a call to
                  nimbleFunction(...).')
    }

    my_particleFilter$setLastTimeRan(0)
    particleMV <- my_particleFilter$mvEWSamples

    # set and calculate
    my_setAndCalculate <- setAndCalculate(model, topParamsInter)
    my_setAndCalculateUpdate <-  mySetAndCalculateUpdate(model = model,
                                                         target = targetAsScalar,
                                                         latents = latents,
                                                         mvSamplesEst = mvSamplesEst,
                                                         my_particleFilter = my_particleFilter,
                                                         m = m,
                                                         topParamsInter = topParamsInter,
                                                         mvSaved = mvSaved
                                                         )

    # Create a vector to store the parameter values
    pfModelValues <- rep(0, length(latentAsScalar))
    storeModelVals <- rep(0, length(targetAsScalar))
    samplerTime <- numeric(length = 4)
    tt <- 0

    # stochastic and deterministic nodes to store
    #ccList1 <- myMcmc_determineCalcAndCopyNodes(model, topParamsInter)
    #copyNodesDeterm1 <- ccList1$copyNodesDeterm; copyNodesStoch1 <- ccList1$copyNodesStoch
    ccList <- myMcmc_determineCalcAndCopyNodes(model, topParamsInter)
    calcNodes <- ccList$calcNodes; copyNodesDeterm1 <- ccList$copyNodesDeterm; copyNodesStoch1 <- ccList$copyNodesStoch   # not used: calcNodesNoSelf
    finalTargetIndex <- max(match(model$expandNodeNames(topParamsInter), calcNodes))
    if(!is.integer(finalTargetIndex) | length(finalTargetIndex) != 1 | is.na(finalTargetIndex[1]))   stop('problem with target node in RW_block sampler')
    calcNodesProposalStage <- calcNodes[1:finalTargetIndex]
    calcNodesDepStage <- calcNodes[-(1:finalTargetIndex)]

    # Make sure calcNodesProposalStage & calcNodesDepStage is not in latent node
    calcNodesProposalStage <- calcNodesProposalStage[!calcNodesProposalStage %in% latentAsScalar]
    calcNodesDepStage <- calcNodesDepStage[!calcNodesDepStage %in% latentAsScalar]

   # print(extraTargetStore)

    ## checks
    if(!inherits(propCov, 'matrix'))                    stop('propCov must be a matrix\n')
    if(!inherits(propCov[1,1], 'numeric'))              stop('propCov matrix must be numeric\n')
    if(!all(dim(propCov) == d))                         stop('propCov matrix must have dimension ', d, 'x', d, '\n')
    if(!isSymmetric(propCov))                           stop('propCov matrix must be symmetric')
    if(length(topParamsInter) < 2)                      stop('less than two top-level targets; cannot use RW_PF_block sampler, try RW_PF sampler')
    if(any(target%in%model$expandNodeNames(latents)))   stop('PMCMC \'target\' argument cannot include latent states')
  },
  run = function() {
    iterRan <<- my_particleFilter$getLastTimeRan()

    if(extras) {
    nimCopy(from = mvSamplesEst, to = model, nodes = extraVars, nodesTo = extraVars, row = iterRan)
    }
    ########################################
    # Top Parameters
    ######################################
    # Copy values of topParams from the previous model
    if(simFirstPars){
      #nimCopy(from = mvSamplesEst, to = model, nodes = topParams, nodesTo = topParams, row = iterRan)
     samplerTime[1] <<- run.time(my_topParamsSampler$run())

    }


    ####################################################
    # Third parameters
    ###################################################
    if(simThirdPars){
      samplerTime[3] <<-  run.time(my_thirdSimParamsSampler$run())
    }

    if(simLastPars){
      samplerTime[4] <<- run.time(my_lastSimParamsSampler$run())
    }


    #############################
    # Intermediate top parameters
    #############################

    # Estimate the loglikehood given the saved latent state samples and
   # tt <<- tt + run.time(storeParticleLP <<- my_setAndCalculateUpdate$run(iterRan))

    tt <<- tt + run.time(storeParticleLP <<- my_particleFilter$getLastLogLik())
    #store latent values and target model parameters from old values
    pfModelValues <<- values(model, latentAsScalar)
    #storeModelVals <<- values(model, targetAsScalar)

    # calculate log likelihood
    tt <<- tt + run.time( modelLP0 <- storeParticleLP + getLogProb(model, calcNodesProposalStage) + getLogProb(model, calcNodesDepStage))

    # generate proposal value
    tt <<- tt + run.time( propValueVector <- generateProposalVector())
    #values(model, topParamsInter) <<- propValueVector
    tt <<- tt + run.time(my_setAndCalculate$run(propValueVector))

    tt <<- tt + run.time(newLP <- getLogProb(model, topParamsInter))

    if(!is.nan(newLP) & (newLP != -Inf)) {
      tt <<- tt + run.time( particleLP <- my_particleFilter$run(m = m, iterRun = iterRan, storeModelValues = values(model, targetAsScalar)))
      tt <<- tt + run.time( modelLP1 <- particleLP + getLogProb(model, calcNodesProposalStage) + getLogProb(model, calcNodesDepStage))
    } else {
      tt <<- tt + run.time(  modelLP1 <- -Inf)
    }

    tt <<- tt + run.time(jump <- my_decideAndJump$run(modelLP1, modelLP0, 0, 0))#, iterRan)#, pfModelValues, predVals, topParamsVals)
#print(jump)
#print(storeParticleLP)
#print(particleLP)
#print(modelLP1)
#print(modelLP0)
    if(jump){
      tt <<- tt + run.time(nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesProposalStage, logProb = TRUE))
      tt <<- tt + run.time(nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm1, logProb = FALSE))
      tt <<- tt + run.time(nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch1, logProbOnly = TRUE))
      tt <<- tt + run.time( index <- ceiling(runif(1, 0, m)))
      tt <<- tt + run.time(copy(particleMV, model, latents, latents, index))
      tt <<- tt + run.time( calculate(model, latentDep))
      tt <<- tt + run.time( copy(from = model, to = mvSaved, nodes = latentDep, row = 1, logProb = TRUE))
    }else{
      tt <<- tt + run.time( my_particleFilter$setLastLogLik(storeParticleLP))
      ## if we don't jump, replace model latent nodes with saved latent nodes
      tt <<- tt + run.time(values(model, latentAsScalar) <<- pfModelValues)
      #values(model, targetAsScalar) <<- storeModelVals
      tt <<- tt + run.time(model$calculate())
      tt <<- tt + run.time(nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesProposalStage, logProb = TRUE))
      tt <<- tt + run.time(nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesDeterm1, logProb = FALSE))
      tt <<- tt + run.time(nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesStoch1, logProbOnly = TRUE))
      #nimCopy(from = mvSaved, to = model, row = 1, nodes = latents, logProb = TRUE)
      tt <<- tt + run.time(nimCopy(from = mvSaved, to = model, nodes = latentDep, row = 1, logProb = TRUE))
      tt <<- tt + run.time(nimCopy(from = model, to = mvSaved, row = 1, nodes = latents, logProb = TRUE))
      # nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm1, logProb = FALSE)
      # nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch1, logProbOnly = TRUE)
      # nimCopy(from = model, to = mvSaved, nodes = latentDep, row = 1, logProb = TRUE)
    }


    if(jump & optimizeM) tt <<- tt + run.time(optimM())
    if(adaptive)    tt <<- tt + run.time( adaptiveProcedure(jump))

    my_particleFilter$setLastTimeRan(iterRan)



    samplerTime[2] <<- tt
    values(model, "samplerTimes") <<- samplerTime[1:4]
    nimCopy(from = model, to = mvSaved, nodes = "samplerTimes", row = 1, logProb = FALSE)
    tt <<- 0

  },
  methods = list(
    optimM = function() {
      tempM <- 15000
      declare(LLEst, double(1, nVarReps))
      if(nVarEsts < mBurnIn) {  # checks whether we have enough var estimates to get good approximation
        for(i in 1:nVarReps)
          LLEst[i] <- my_particleFilter$run(m = tempM, iterRun = iterRan, storeModelValues = values(model, targetAsScalar))
        ## next, store average of var estimates
        if(nVarEsts == 1)
          storeLLVar <<- var(LLEst)/mBurnIn
        else {
          LLVar <- storeLLVar
          LLVar <- LLVar + var(LLEst)/mBurnIn
          storeLLVar <<- LLVar
        }
        nVarEsts <<- nVarEsts + 1
      }
      else {  # once enough var estimates have been taken, use their average to compute m
        m <<- m*storeLLVar/(0.92^2)
        m <<- ceiling(m)
        storeParticleLP <<- my_particleFilter$run(m, iterRun = iterRan, storeModelValues = values(model, targetAsScalar))
        optimizeM <<- 0
      }
    },
    generateProposalVector = function() {
      propValueVector <- rmnorm_chol(1, values(model,topParamsInter), chol_propCov_scale, 0)  ## last argument specifies prec_param = FALSE
      returnType(double(1))
      return(propValueVector)
    },
    generateReducedVals = function(){

      lp <- my_setAndCalculateUpdate$run(iterRan)
      returnType(double())
      return(lp)
    },
    adaptiveProcedure = function(jump = logical()) {
      timesRan <<- timesRan + 1
      if(jump)     timesAccepted <<- timesAccepted + 1
      if(!adaptScaleOnly)     empirSamp[timesRan, 1:d] <<- values(model, topParamsInter)
      if(timesRan %% adaptInterval == 0) {
        acceptanceRate <- timesAccepted / timesRan
        timesAdapted <<- timesAdapted + 1
        if(saveMCMChistory) {
          setSize(scaleHistory, timesAdapted)                 ## scaleHistory
          scaleHistory[timesAdapted] <<- scale                ## scaleHistory
          setSize(acceptanceHistory, timesAdapted)            ## scaleHistory
          acceptanceHistory[timesAdapted] <<- acceptanceRate  ## scaleHistory
          if(d <= maxDimCovHistory) {
            propCovTemp <- propCovHistory                                           ## scaleHistory
            setSize(propCovHistory, timesAdapted, d, d)                             ## scaleHistory
            if(timesAdapted > 1)                                                    ## scaleHistory
              for(iTA in 1:(timesAdapted-1))                                      ## scaleHistory
                propCovHistory[iTA, 1:d, 1:d] <<- propCovTemp[iTA, 1:d, 1:d]    ## scaleHistory
            propCovHistory[timesAdapted, 1:d, 1:d] <<- propCov[1:d, 1:d]            ## scaleHistory
          }
        }
        adaptFactor <- my_calcAdaptationFactor$run(acceptanceRate)
        scale <<- scale * adaptFactor
        ## calculate empirical covariance, and adapt proposal covariance
        if(!adaptScaleOnly) {
          gamma1 <- my_calcAdaptationFactor$getGamma1()
          for(i in 1:d)     empirSamp[, i] <<- empirSamp[, i] - mean(empirSamp[, i])
          empirCov <- (t(empirSamp) %*% empirSamp) / (timesRan-1)
          propCov <<- propCov + gamma1 * (empirCov - propCov)
          chol_propCov <<- chol(propCov)
        }
        chol_propCov_scale <<- chol_propCov * scale
        timesRan <<- 0
        timesAccepted <<- 0
      }
    },
    setScale = function(newScale = double()) {
      scale         <<- newScale
      scaleOriginal <<- newScale
      chol_propCov_scale <<- chol_propCov * scale
    },
    setPropCov = function(newPropCov = double(2)) {
      propCov         <<- newPropCov
      propCovOriginal <<- newPropCov
      chol_propCov <<- chol(propCov)
      chol_propCov_scale <<- chol_propCov * scale
    },
    getScaleHistory = function() {  ## scaleHistory
      if(!saveMCMChistory)   print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC.")
      returnType(double(1))
      return(scaleHistory)
    },
    getAcceptanceHistory = function() {  ## scaleHistory
      returnType(double(1))
      if(!saveMCMChistory)   print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC.")
      return(acceptanceHistory)
    },
    getPropCovHistory = function() { ## scaleHistory
      if(!saveMCMChistory | d > maxDimCovHistory)   print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC.  Note that to reduce memory use, proposal covariance histories are only saved for parameter vectors of length <= 10; this value can be modified using the 'maxDimCovHistory' control list element.")
      returnType(double(3))
      return(propCovHistory)
    },
    ##getScaleHistoryExpanded = function() {                                                              ## scaleHistory
    ##    scaleHistoryExpanded <- numeric(timesAdapted*adaptInterval, init=FALSE)                         ## scaleHistory
    ##    for(iTA in 1:timesAdapted)                                                                      ## scaleHistory
    ##        for(j in 1:adaptInterval)                                                                   ## scaleHistory
    ##            scaleHistoryExpanded[(iTA-1)*adaptInterval+j] <- scaleHistory[iTA]                      ## scaleHistory
    ##    returnType(double(1)); return(scaleHistoryExpanded) },                                          ## scaleHistory
    ##getPropCovHistoryExpanded = function() {                                                            ## scaleHistory
    ##    propCovHistoryExpanded <- array(dim=c(timesAdapted*adaptInterval,d,d), init=FALSE)              ## scaleHistory
    ##    for(iTA in 1:timesAdapted)                                                                      ## scaleHistory
    ##        for(j in 1:adaptInterval)                                                                   ## scaleHistory
    ##            propCovHistoryExpanded[(iTA-1)*adaptInterval+j,1:d,1:d] <- propCovHistory[iTA,1:d,1:d]  ## scaleHistory
    ##    returnType(double(3)); return(propCovHistoryExpanded) },                                        ## scaleHistory
    reset = function() {
      scale   <<- scaleOriginal
      propCov <<- propCovOriginal
      chol_propCov <<- chol(propCov)
      chol_propCov_scale <<- chol_propCov * scale
      timesRan      <<- 0
      timesAccepted <<- 0
      timesAdapted  <<- 0
      if(saveMCMChistory) {
        scaleHistory  <<- c(0, 0)    ## scaleHistory
        acceptanceHistory  <<- c(0, 0)
        if(d <= maxDimCovHistory)
          propCovHistory <<- nimArray(0, dim = c(2,d,d))
      }
      my_calcAdaptationFactor$reset()
    }
  )
)
