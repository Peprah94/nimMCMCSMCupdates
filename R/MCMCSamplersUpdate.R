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
    adaptive            <- extractControlElement(control, 'adaptive',             TRUE)
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

    #############################
    # Obtaining latent states, hyperparameters (topParams) and intermediate parameters and the
    #their dependencies
    ##########################
    ## node list generation
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)

    # Get dependencies of the target variables
    calcNodes <- model$getDependencies(target)

    # We always sample the latent states with this application and get their dependencies
    latentSamp <- TRUE
    latentDep <- model$getDependencies(latents)
    # Get the paramters we are sampling
    # MCMCmonitors <- tryCatch(parent.frame(2)$conf$monitors, error = function(e) e)
    # if(identical(MCMCmonitors, TRUE))
    #   latentSamp <- TRUE
    # else if(any(model$expandNodeNames(latents) %in% model$expandNodeNames(MCMCmonitors)))
    #   latentSamp <- TRUE


    ## Now we try to seperate the hyperparameters from the other model parameters
    ## We do this by looking at the model graph and going one step back at a time.
    ## Currently, we want to update model parameter that are greater than 5
    ## This is not an efficient approach, but it works for the application in the manuscript.

    # Get top parameters
    topParams <- model$getNodeNames(stochOnly=TRUE, includeData=FALSE, topOnly=TRUE)


    if(length(topParams) <5){
  topParams <- topParams
  trueParamsFALSE <- topParams
}else{
    #Extract the true top parameters (hyperparameters). Some of the top parameters depend
    # are intermediate
    trueTopParams <- sapply(topParams, function(x){
      # Is the target variable time dependent? Its a logical value with FALSE for time dependence
      nodeExpandedTimeDependent <- length(model$expandNodeNames(model$getVarNames(nodes = x))) == length(x)
      #get dependencies of each top Parameter
      if(nodeExpandedTimeDependent){
      specTargetDeps <- model$getDependencies(x, includeData = FALSE, self = FALSE, stochOnly = TRUE)
      #get the variable name
      varNames <- model$getVarNames(nodes = specTargetDeps)
      #Find out of it is part of topParameters or latent state
      isLatent <- varNames[varNames%in%c(latents, topParams[!topParams %in% x])]
      #Number of intermediate parameters (other model parameters)
      numRep <- length(specTargetDeps)
      ret <- ifelse(numRep == 0 | length(isLatent) == 0, TRUE, FALSE)
      }else{
        ret <- FALSE
      }
      return(ret)
    })

    ## The above gives an index of the parameters that are hyperparameters and those that are not
    topParamsTRUE <- topParams[unlist(trueTopParams )]
    trueParamsFALSE <- topParams[!unlist(trueTopParams )]

    # Is the intermediate variable in extra vars? Extra vars are other stochastic variables we do not intend to update
    isInterInExtraVars <- trueParamsFALSE %in% model$expandNodeNames(extraVars)

    # For consistency with the nimbleSMC work, we assign the hyperparameters as the topParameters
    topParams <- topParamsTRUE

    print(topParams)


    # Get dependencies of top Parameters that does not include itself and data components
    topParamsDeps <- model$getDependencies(topParams, self = FALSE, includeData = FALSE, stochOnly = TRUE)

     # Get dependencies of top Parameters that does not include itself and data components
    calcNodesTopParams <- model$getDependencies(topParams)
    #Get the data Node that includes the vars and the data itself
    dataNodes <- model$getVarNames(nodes = model$getDependencies(latents, downstream = TRUE, stochOnly = TRUE, self = FALSE))

    # Now we focus on the extra variables
    if(!is.null(extraVars)){
      #extraVars <- model$expandNodeNames(node = extraVars)
      #seqExtraVars <- paste0("[[", 1:iNodePrev, "]]")
      extraTargetVarsWithAdditionalPars <- lapply(seq_along(extraVars), function(x){
        extraVarsScalar <- model$expandNodeNames(node = extraVars[[x]])
        extraVarsScalar[1 : iNodePrev]
        #extraVars[!grepl(seqExtraVars[[x]], extraVars)]
      })%>%
        do.call("c",.)
    }
    #Intermediary parameters
    if(is.null(extraVars)){
    topParamsInter <- unique(c(topParamsDeps[!topParamsDeps %in% c(model$expandNodeNames(latents), model$expandNodeNames(dataNodes), model$expandNodeNames(additionalPars))], trueParamsFALSE))
    }else{
      if(!isInterInExtraVars){
      topParamsInter <- unique(c(topParamsDeps[!topParamsDeps %in% c(model$expandNodeNames(latents),
                                                                     model$expandNodeNames(dataNodes),
                                                                     model$expandNodeNames(additionalPars),
                                                                     extraTargetVarsWithAdditionalPars
                                                                     )],
                                 trueParamsFALSE))}else{
                                   topParamsInter <- unique(c(topParamsDeps[!topParamsDeps %in% c(model$expandNodeNames(latents),
                                                                                                 model$expandNodeNames(dataNodes),
                                                                                                 model$expandNodeNames(additionalPars),
                                                                                                  extraTargetVarsWithAdditionalPars
                                   )]))
                                 }
    }
    print(topParamsInter)
}

    #Extra vars to simulate
    #extraTargetVars <- topParamsInter[!grepl("[[1]]", topParamsInter)]
    # extraTargetVars <- extraTargetVars[!extraTargetVars %in% topParamsInter]

    #target <- model$expandNodeNames(target)

    #Whether we should use multiple
    #if(all(topParams %in% target) == TRUE){
    # if(){
    #  multiple <- TRUE
    # }else{
    #   multiple <- FALSE
    # }

    # check whether topParams is the same as length
    NotMultiple <- length(topParams) == length(target)

    if(NotMultiple == FALSE){
      multiple = TRUE
    }else{
      multiple = FALSE
    }

    print(multiple)

    #get the extra time dependent vars
    #uncomment to revert to old
    #reducedTarget <- reducedModel$expandNodeNames(target)


    #if we have don't have multiple,
    #let topParamsInter = targetAsScalar
    if(NotMultiple == TRUE) topParamsInter <- targetAsScalar
    topParamsInterDep <- model$getDependencies(topParamsInter, self = FALSE,  includeData = FALSE, stochOnly = TRUE)
    calcNodesTopParamsInter <- model$getDependencies(topParamsInterDep)
    storeOtherVals <- FALSE
    if(multiple){
      if(!is.null(extraVars)){
      #extraVars <- model$expandNodeNames(node = extraVars)
      #seqExtraVars <- paste0("[[", 1:iNodePrev, "]]")
      extraTargetVars <- lapply(seq_along(extraVars), function(x){
        extraVarsScalar <- model$expandNodeNames(node = extraVars[[x]])
        ret <- extraVarsScalar[(iNodePrev +1) : length(extraVarsScalar)]
        return(ret[!is.na(ret)])
        #extraVars[!grepl(seqExtraVars[[x]], extraVars)]
      })%>%
        do.call("c",.)
      # We do not intend to update these parameters, although they are time dependent.
      extraAdditionalPars <- model$expandNodeNames(additionalPars[!additionalPars %in% latents])
  # Get extra variables that are stochastic
      stochExtraPars <- extraAdditionalPars[model$isStoch(extraAdditionalPars)]
     # Get extra variables that are deterministic but time dependent
       detExtraPars <- extraAdditionalPars[!model$isStoch(extraAdditionalPars)]#[1:iNodePrev]
       detExtraParsVarNames <- model$getVarNames(nodes = detExtraPars)
       #filterControl$timeIndex > 1 &

         nDimLatents <- dim(mvSamplesEst[[latents]][[1]])[2]
       detExtras <- lapply(detExtraParsVarNames, function(x){
         namesNodes <- model$expandNodeNames(x)
         if(filterControl$timeIndex>1 & length(namesNodes)==nDimLatents){
           ret <- namesNodes[1:iNodePrev]
         }else{
           ret <- namesNodes
         }
         return(ret)
       })%>%
         do.call("c", .)

       #stochExtraParsNodeName <- model$getVarNames(nodes = stochExtraPars)
      #extraTargetVars <- extraVars[!grepl(, extraVars)]
     extraTargetStore <-  c(extraTargetVarsWithAdditionalPars,
                            stochExtraPars,
                            detExtras)
     storeOtherVals <- TRUE
     extras <- TRUE
      }else{
        extraTargetVars <- latents
        extraTargetStore <- latents
        extras <- FALSE
      }
    }else{
      extraTargetVars <- targetAsScalar
      extraTargetStore <- targetAsScalar
      extras <- FALSE
    }
#print(extraTargetVars)
#print(storeOtherVals)
#print(extraTargetStore)

    #if(multiple) target <- c(topParams, topParamsInter)
    ## numeric value generation
    optimizeM     <- as.integer(optimizeM)
    scaleOriginal <- scale
    timesRan      <- 0
    timesAccepted <- 0
    timesAdapted  <- 0
    prevLL        <- 0
    nVarEsts      <- 0
    itCount       <- 0
    # iterRun <- 1
    d <- length(topParamsInter)
    # d1 <- length(topParams)
    #dTopPars <- length(topParams)
    if(is.character(propCov) && propCov == 'identity')     propCov <- diag(d)
    propCovOriginal <- propCov
    chol_propCov <- chol(propCov)
    chol_propCov_scale <- scale * chol_propCov
    #if(multiple){

    #}

    empirSamp <- matrix(0, nrow=adaptInterval, ncol=d)
    storeParticleLP <- -Inf
    iterRan <- 1
    storeLLVar  <- 0
    nVarReps <- 7    ## number of LL estimates to compute to get each LL variance estimate for m optimization
    mBurnIn  <- 15   ## number of LL variance estimates to compute before deciding optimal m
    if(optimizeM)   m <- 3000
    ## nested function and function list definitions
    # newModel <- model$newModel(replicate = TRUE)
    #my_setAndCalculate <- setAndCalculate(model, topParamsInter)
    #my_setAndCalculateUpdate <-  mySetAndCalculateUpdate(model, target, latents, mvSamplesEst, my_particleFilter, m, topParamsInter, extraTargetVars)
    my_decideAndJump <-  myDecideAndJump(model, mvSaved, topParamsInter,latentAsScalar, mvSamplesEst)
    my_calcAdaptationFactor <- calcAdaptationFactor(d, adaptFactorExponent)
    #if(multiple)
    #mvSavedNew <- modelValues(model)


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
        my_particleFilter <- buildAuxiliaryFilterNew(model, latents, control = filterControl)
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
        my_particleFilter <- buildBootstrapFilterNew(model, latents, control = filterControl)

      }
      else if(is.nfGenerator(filterType)){
        my_particleFilter <- filterType(model, latents, control = filterControl)
      }
      else stop('filter type must be either "bootstrap", "auxiliary", or a
                  user defined filtering algorithm created by a call to
                  nimbleFunction(...).')
    }

    #Save the latents as scalars for the model fitted with the reduced model parameters

    # Latent states
    latentAsScalar <- model$expandNodeNames(latents, returnScalarComponents = TRUE)
    #predictivePars <- target[!target %in% topParams]

    # if(!multiple){
    #  topParamsDeps <- target
    #}
    #extraParsToSave <- calcNodes[!calcNodes %in% model$expandNodeNames(c(target, latents))]
    #my_setAndCalculate <- setAndCalculate(model, target)
    my_setAndCalculate <- setAndCalculate(model, topParamsInter)
    my_setAndCalculateUpdate <-  mySetAndCalculateUpdate(model = model,
                                                         target = target,
                                                         latents = latents,
                                                         mvSamplesEst = mvSamplesEst,
                                                         my_particleFilter = my_particleFilter,
                                                         m = m,
                                                         topParamsInter = topParamsInter,
                                                         mvSaved = mvSaved,
                                                         extraTargetVars = extraTargetVars,
                                                         extras = extras)
    #my_decideAndJump <-  myDecideAndJump(model, mvSaved, topParamsInter,latentAsScalar, mvSamplesEst)
    #my_calcAdaptationFactor <- calcAdaptationFactor(d, adaptFactorExponent)
    #if(multiple)
    #my_sampleTopPars <- sampleTopPars(model, mvSaved, topParams, mvSamplesEst, scale, latents)

    #set time run = 0
    my_particleFilter$setLastTimeRan(0)
    #my_particleFilter1$setLastTimeRan(0)
    particleMV <- my_particleFilter$mvEWSamples

    # particleMVold <- my_particleFilter$mvEWSamples
    #pfModelValues <- rep(0, length(model$expandNodeNames(latentAsScalar)))
    pfModelValues <- rep(0, length(latentDep))
    targetModelValues <- rep(0, length(topParamsInter))
    topParamsValues <- rep(0, length(topParams))
    storeModelVals <- rep(0, length(targetAsScalar))
    # topParsModelValues <- rep(0, length(topParams))
    #predVals <- rep(0, length(predictivePars))
    #topParamsVals <- rep(0, length(topParams))

    #predictive nodes
    #simNodes  <- model$getDependencies(target, downstream = TRUE, includePredictive = TRUE)
    # calcNodes <- model$getDependencies(target, downstream = TRUE, includePredictive = TRUE, stochOnly = TRUE)

    ccList1 <- myMcmc_determineCalcAndCopyNodes(model, topParamsInter)
    copyNodesDeterm1 <- ccList1$copyNodesDeterm; copyNodesStoch1 <- ccList1$copyNodesStoch
print(extraTargetStore)
    #initialize the decision process
    #jump <- 0
    #pfModelValues <- rep(0, length(latents))
    # if(extraSave){
    #   pfModelValuesForExtra <- rep(0, length(calNodesStoch))
    # }
# print(target)
# print(targetAsScalar)
# print(topParamsInter)
# print(extraTargetVars)
    #ccList <- myMcmc_determineCalcAndCopyNodes(model, target)
    # copyNodesDeterm <- ccList$copyNodesDeterm; copyNodesStoch <- ccList$copyNodesStoch
    #pfNewModelValues <- rep(0, length(model$expandNodeNames(latents)))

    # saveOldVars <- modelValues(modelValuesConf(vars = my_particleFilter$mvEWSamples$getVarNames(),
    #                                           types = "double",
    #                                           sizes = length(model$expandNodeNames(latents))))
    # targetNodesAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    #storeModelValues <- values(model, targetNodesAsScalar)
    ## checks
    if(!inherits(propCov, 'matrix'))                    stop('propCov must be a matrix\n')
    if(!inherits(propCov[1,1], 'numeric'))              stop('propCov matrix must be numeric\n')
    if(!all(dim(propCov) == d))                         stop('propCov matrix must have dimension ', d, 'x', d, '\n')
    if(!isSymmetric(propCov))                           stop('propCov matrix must be symmetric')
    if(length(targetAsScalar) < 2)                      stop('less than two top-level targets; cannot use RW_PF_block sampler, try RW_PF sampler')
    if(any(target%in%model$expandNodeNames(latents)))   stop('PMCMC \'target\' argument cannot include latent states')

    d1 <- length(topParams)
    #propCov1 <- diag(d1)
    if(is.character(propCov1) && propCov1 == 'identity')     propCov1 <- diag(d)
    my_sampleTopPars <- sampleTopPars(model, mvSaved, topParams, mvSamplesEst, scale, latents, target,propCov1)
    },
  run = function() {
    iterRan <<- my_particleFilter$getLastTimeRan()
    #print(iterRan)
    nimCopy(from = mvSamplesEst, to = model, nodes = target, nodesTo = target, row = iterRan)
    #print(values(model, targetAsScalar))
    # Update top Pars
    if(storeOtherVals)  nimCopy(from = mvSamplesEst, to = model, nodes = extraTargetStore, nodesTo = extraTargetStore, row = iterRan, rowTo = 1)
    # Update top Pars
    #propValueVectorTopPars <- generateProposalVector(topParams)
    #MHAR for top Level pars
    if(multiple) my_sampleTopPars$run(iterRan)
    #print(2)
    #MHAR for additional pars
  storeParticleLP <<- my_setAndCalculateUpdate$run(iterRan)

    #print(storeParticleLP)
    #store latent values and target model parameters from old values
    pfModelValues <<- values(model, latentDep)
    storeModelVals <<- values(model, targetAsScalar)
   # print(all(storeModelVals == )
    targetModelValues <<- values(model, topParamsInter)
    # print(4)
    if(multiple) topParamsValues <<- values(model, topParams)
    #print(topParamsValues)
    #predVals <<- values(model, predictivePars)
    #topParamsVals <<- values(model, topParams)
    #print(5)
    modelLP0 <- storeParticleLP + getLogProb(model, topParamsInter)
    propValueVector <- generateProposalVector()
   # print(propValueVector)
    #print(6)
    my_setAndCalculate$run(propValueVector)
    #add the topParams values to the proposal vector since it will enter the particle filter
    #if(multiple) storeModelVals <<- c(topParamsValues, propValueVector)
    #if(!multiple) storeModelVals <<- propValueVector
    #print(7)
    particleLP <- my_particleFilter$run(m = m, iterRun = iterRan, storeModelValues = values(model, targetAsScalar))
   #print(particleLP)
    modelLP1 <- particleLP + getLogProb(model, topParamsInter)
    #print(9)
    jump <- my_decideAndJump$run(modelLP1, modelLP0, 0, 0, iterRan)#, pfModelValues, predVals, topParamsVals)
#print(jump)
    # if(!jump) {
    #   my_particleFilter$setLastLogLik(storeParticleLP)
    # }
    #print(10)
    #print(jump)
    if(jump){#& latentSamp) {
      ## if we jump, randomly sample latent nodes from pf output and put
      ## into model so that they can be monitored
      nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm1, logProb = FALSE)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch1, logProbOnly = TRUE)
      index <- ceiling(runif(1, 0, m))
      copy(particleMV, model, latents, latents, index)
      calculate(model, latentDep)
      copy(from = model, to = mvSaved, nodes = latentDep, row = 1, logProb = TRUE)
      #my_particleFilter$setLastLogLik(particleLP)
      }else{#& latentSamp) {
      ## if we don't jump, replace model latent nodes with saved latent nodes
      values(model, latentDep) <<- pfModelValues
      values(model, targetAsScalar) <<- storeModelVals
      if(multiple) nimCopy(from = mvSamplesEst, to = model, nodes = topParams,nodesTo = topParams, row = iterRan, rowTo = 1)
      #calculate(model)
      model$calculate()
      #copy(from = model, to = mvSaved, nodes = latent)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm1, logProb = FALSE)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch1, logProbOnly = TRUE)
   #   nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
      nimCopy(from = model, to = mvSaved, nodes = latentDep, row = 1, logProb = TRUE)
      #copy(from = model, to = mvSaved, nodes = latents, row = 1, logProb = TRUE)
      #nimCopy(from = mvSaved, to = model, nodes = latentDep, row = 1, logProb = TRUE)
      #storeParticleLP <- my_particleFilter$run(m = m, iterRun = (iterRan+1), storeModelValues = values(model, targetAsScalar))
     # my_particleFilter$setLastLogLik(storeParticleLP)
      }
    if(storeOtherVals) nimCopy(from = mvSamplesEst, to = mvSaved, row = iterRan, rowTo = 1, nodes = extraTargetStore, nodesTo = extraTargetStore)
    ##if(jump & !resample)  storeParticleLP <<- particleLP
    if(jump & optimizeM) optimM()
    if(adaptive)     adaptiveProcedure(jump)

    my_particleFilter$setLastTimeRan(iterRan)
    #generateIterRun()

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
      # if(!jump & !multiple){
      #   index <- ceiling(runif(1, 0, m))
      #   lp <- my_particleFilter$run(m = m,
      #                               iterRun = iterRan,
      #                               storeModelValues = values(model, targetAsScalar))
      #   copy(particleMV, model, latents, latents, index)
      # }else if(jump & !multiple){
      #   lp <- my_particleFilter$getLastLogLik()
      #   }
      lp <- my_setAndCalculateUpdate$run(iterRan)
     # if(jump & !multiple)
       # lp <- my_particleFilter$getLastLogLik()


      #}
      returnType(double())
      return(lp)
    },
    adaptiveProcedure = function(jump = logical()) {
      timesRan <<- timesRan + 1
      if(jump)     timesAccepted <<- timesAccepted + 1
      if(!adaptScaleOnly)     empirSamp[timesRan, 1:d] <<- values(model, target)
      if((timesRan-1) %% adaptInterval == 0) {
        acceptanceRate <- timesAccepted / (timesRan-1)
        timesAdapted <<- timesAdapted + 1
        adaptFactor <- my_calcAdaptationFactor$run(acceptanceRate)
        scale <<- scale * adaptFactor
        ## calculate empirical covariance, and adapt proposal covariance
        if(!adaptScaleOnly) {
          gamma1 <- my_calcAdaptationFactor$getGamma1()
          for(i in 1:d)     empirSamp[, i] <<- empirSamp[, i] - mean(empirSamp[, i])
          empirCov <- (t(empirSamp) %*% empirSamp) / (timesRan-2)
          propCov <<- propCov + gamma1 * (empirCov - propCov)
          chol_propCov <<- chol(propCov)
        }
        chol_propCov_scale <<- chol_propCov * scale
        timesRan <<- 0
        #iterRun <<- 1
        timesAccepted <<- 0
      }
    },
    reset = function() {
      scale   <<- scaleOriginal
      propCov <<- propCovOriginal
      chol_propCov <<- chol(propCov)
      chol_propCov_scale <<- chol_propCov * scale
      storeParticleLP <<- -Inf
      timesRan      <<- 1
      timesAccepted <<- 0
      timesAdapted  <<- 0
      iterRan <<- 0
      my_calcAdaptationFactor$reset()
    }
  )
)
#


#' Particle Filtering MCMC Sampling Algorithms
#'
#' Details of the particle filtering MCMC sampling algorithms provided in nimbleSMC.
#'
#' @param model (uncompiled) model on which the MCMC is to be run
#' @param mvSaved \code{modelValues} object to be used to store MCMC samples
#' @param target node(s) on which the sampler will be used
#' @param control named list that controls the precise behavior of the sampler, with elements specific to \code{samplertype}.  The default values for control list are specified in the setup code of each sampling algorithm.  Descriptions of each sampling algorithm, and the possible customizations for each sampler (using the \code{control} argument) appear below.
#'
#' @section RW_PF sampler:
#'
#' The particle filter sampler allows the user to perform particle MCMC (PMCMC) (Andrieu et al., 2010), primarily for state-space or hidden Markov models of time-series data. This method uses Metropolis-Hastings samplers for top-level parameters but uses the likelihood approximation of a particle filter (sequential Monte Carlo) to integrate over latent nodes in the time-series.  The \code{RW_PF} sampler uses an adaptive Metropolis-Hastings algorithm with a univariate normal proposal distribution for a scalar parameter.  Note that samples of the latent states can be retained as well, but the top-level parameter being sampled must be a scalar.   A bootstrap, auxiliary, or user defined particle filter can be used to integrate over latent states.
#'
#' For more information about user-defined samplers within a PMCMC sampler, see the NIMBLE User Manual.
#'
#' The \code{RW_PF} sampler accepts the following control list elements:
#' \itemize{
#' \item adaptive. A logical argument, specifying whether the sampler should adapt the scale (proposal standard deviation) throughout the course of MCMC execution to achieve a theoretically desirable acceptance rate. (default = TRUE)
#' \item adaptInterval. The interval on which to perform adaptation.  Every adaptInterval MCMC iterations (prior to thinning), the RW sampler will perform its adaptation procedure.  This updates the scale variable, based upon the sampler's achieved acceptance rate over the past adaptInterval iterations. (default = 200)
#' \item scale. The initial value of the normal proposal standard deviation.  If \code{adaptive = FALSE}, scale will never change. (default = 1)
#' \item pfNparticles.  The number of particles to use in the approximation to the log likelihood of the data (default = 1000).
#' \item latents.  Character vector specifying the nodes that are latent states over which the particle filter will operate to approximate the log-likelihood function.
#' \item pfType.  Character argument specifying the type of particle filter that should be used for likelihood approximation.  Choose from \code{"bootstrap"} and \code{"auxiliary"}.  Defaults to \code{"bootstrap"}.
#' \item pfControl.  A control list that is passed to the particle filter function.  For details on control lists for bootstrap or auxiliary particle filters, see \code{\link{buildBootstrapFilter}} or \code{\link{buildAuxiliaryFilter}} respectively.  Additionally, this can be used to pass custom arguments into a user-defined particle filter.
#' \item pfOptimizeNparticles.  A logical argument, specifying whether to use an experimental procedure to automatically determine the optimal number of particles to use, based on Pitt and Shephard (2011).  This will override any value of \code{pfNparticles} specified above.
#' \item pf.  A user-defined particle filter object, if a bootstrap or auxiliary particle filter is not adequate.
#' }
#'
#' @section RW_PF_block sampler:
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
### RW_PF, does a univariate RW, but using a particle filter likelihood function ######
#######################################################################################

#' @rdname samplers
#' @export
# sampler_RW_PFUpdate <- nimbleFunction(
#   name = 'sampler_RW_PFUpdate',
#   contains = sampler_BASE,
#   setup = function(model, mvSaved, target, mvWSamplesWTSaved,
#                    mvWSamplesXSaved, mvEWSamplesXSaved, logLikeVals,
#                    control) {
#     ## control list extraction
#     adaptive       <- extractControlElement(control, 'adaptive',             TRUE)
#     adaptInterval  <- extractControlElement(control, 'adaptInterval',        200)
#     scale          <- extractControlElement(control, 'scale',                1)
#     m              <- extractControlElement(control, 'pfNparticles',         1000)
#     existingPF     <- extractControlElement(control, 'pf',                   NULL)
#     filterType     <- extractControlElement(control, 'pfType',               'bootstrap')
#     filterControl  <- extractControlElement(control, 'pfControl',            list())
#     optimizeM      <- extractControlElement(control, 'pfOptimizeNparticles', FALSE)
#     latents        <- extractControlElement(control, 'latents',              error = 'RW_PF sampler missing required control argument: latents')
#     if('pfLookahead' %in% names(control)) {
#       warning("The `pfLookahead` control list argument is deprecated and will not be supported in future versions of `nimbleSMC`. Please specify the lookahead function via the pfControl argument instead.")
#       filterControl$lookahead <- control[['pfLookahead']]
#     }
#     else if(!('lookahead' %in% names(filterControl))) {
#       filterControl$lookahead <- 'simulate'
#     }
#     ## node list generation
#     targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
#     calcNodes <- model$getDependencies(target)
#     latentSamp <- FALSE
#     MCMCmonitors <- tryCatch(parent.frame(2)$conf$monitors, error = function(e) e)
#     if(identical(MCMCmonitors, TRUE))
#       latentSamp <- TRUE
#     else if(any(model$expandNodeNames(latents) %in% model$expandNodeNames(MCMCmonitors)))
#       latentSamp <- TRUE
#     latentDep <- model$getDependencies(latents)
#     topParams <- model$getNodeNames(stochOnly=TRUE, includeData=FALSE, topOnly=TRUE)
#     ## numeric value generation
#     optimizeM       <- as.integer(optimizeM)
#     scaleOriginal   <- scale
#     times <- 1
#     timesRan        <- 1
#     timesAccepted   <- 0
#     timesAdapted    <- 0
#     prevLL          <- 0
#     nVarEsts        <- 0
#     itCount         <- 0
#     optimalAR       <- 0.44
#     gamma1          <- 0
#     storeParticleLP <- -Inf
#     storeLLVar      <- 0
#     nVarReps        <- 7    ## number of LL estimates to compute to get each LL variance estimate for m optimization
#     mBurnIn         <- 15   ## number of LL variance estimates to compute before deciding optimal m
#     d               <- length(targetAsScalar)
#     if(optimizeM) m <- 3000
#     ## Nested function and function list definitions.
#     my_setAndCalculate <- setAndCalculateOne(model, target)
#     my_decideAndJump <- decideAndJump(model, mvSaved, target, calcNodes)
#     if(!is.null(existingPF)) {
#       my_particleFilter <- existingPF
#     } else {
#       if(latentSamp == TRUE) {
#         filterControl$saveAll <- TRUE
#         filterControl$smoothing <- TRUE
#       } else {
#         filterControl$saveAll <- FALSE
#         filterControl$smoothing <- FALSE
#       }
#       filterControl$initModel <- TRUE
#       if(is.character(filterType) && filterType == 'auxiliary') {
#         my_particleFilter <- buildAuxiliaryFilter(model, latents,
#                                                   control = filterControl)
#       }
#       else if(is.character(filterType) && filterType == 'bootstrap') {
#         my_particleFilter <- buildBootstrapFilter(model, latents,
#                                                   control = filterControl)
#       }
#       else if(is.nfGenerator(filterType)){
#         my_particleFilter <- filterType(model, latents,
#                                         control = filterControl)
#       }
#       else stop('filter type must be either "bootstrap", "auxiliary", or a
#                   user defined filtering algorithm created by a call to
#                   nimbleFunction(...).')
#     }
#     particleMV <- my_particleFilter$mvEWSamples
#     ## checks
#     if(any(target%in%model$expandNodeNames(latents)))   stop('PMCMC \'target\' argument cannot include latent states')
#     if(length(targetAsScalar) > 1)                      stop('more than one top-level target; cannot use RW_PF sampler, try RW_PF_block sampler')
#   },
#   run = function() {
#     #
#     storeParticleLP <<- my_particleFilter$getLastLogLik()
#     modelLP0 <- storeParticleLP + getLogProb(model, target)
#     propValue <- rnorm(1, mean = model[[target]], sd = scale)
#     my_setAndCalculate$run(propValue)
#     particleLP <- my_particleFilter$run(m, iterRun = timesRan, storeModelValues = propValue)
#     modelLP1 <- particleLP + getLogProb(model, target)
#     jump <- my_decideAndJump$run(modelLP1, modelLP0, 0, 0)
#     #times <<- times + 1
#     if(!jump) {
#       my_particleFilter$setLastLogLik(storeParticleLP)
#     }
#     if(jump & latentSamp){
#       ## if we jump, randomly sample latent nodes from pf output and put into model so that they can be monitored
#       index <- ceiling(runif(1, 0, m))
#       copy(particleMV, model, latents, latents, index)
#       calculate(model, latentDep)
#       copy(from = model, to = mvSaved, nodes = latentDep, row = 1, logProb = TRUE)
#     }
#     else if(!jump & latentSamp){
#       ## if we don't jump, replace model latent nodes with saved latent nodes
#       copy(from = mvSaved, to = model, nodes = latentDep, row = 1, logProb = TRUE)
#     }
#     ##if(jump & !resample)  storeParticleLP <<- particleLP
#     if(jump & optimizeM) optimM()
#     if(adaptive)     adaptiveProcedure(jump)
#   },
#   methods = list(
#     optimM = function() {
#       tempM <- 15000
#       declare(LLEst, double(1, nVarReps))
#       if(nVarEsts < mBurnIn) {  # checks whether we have enough var estimates to get good approximation
#         for(i in 1:nVarReps)
#           LLEst[i] <- my_particleFilter$run(tempM, iterRun = 1)
#         ## next, store average of var estimates
#         if(nVarEsts == 1)
#           storeLLVar <<- var(LLEst)/mBurnIn
#         else {
#           LLVar <- storeLLVar
#           LLVar <- LLVar + var(LLEst)/mBurnIn
#           storeLLVar<<- LLVar
#         }
#         nVarEsts <<- nVarEsts + 1
#       }
#       else {  # once enough var estimates have been taken, use their average to compute m
#         m <<- m*storeLLVar/(0.92^2)
#         m <<- ceiling(m)
#         storeParticleLP <<- my_particleFilter$run(m, iterRun = 1)
#         optimizeM <<- 0
#       }
#     },
#     adaptiveProcedure = function(jump = logical()) {
#       timesRan <<- timesRan + 1
#       if(jump)     timesAccepted <<- timesAccepted + 1
#       if(timesRan %% adaptInterval == 0) {
#         acceptanceRate <- timesAccepted / timesRan
#         timesAdapted <<- timesAdapted + 1
#         gamma1 <<- 1/((timesAdapted + 3)^0.8)
#         gamma2 <- 10 * gamma1
#         adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
#         scale <<- scale * adaptFactor
#         timesRan <<- 0
#         timesAccepted <<- 0
#       }
#     },
#     reset = function() {
#       scale <<- scaleOriginal
#       timesRan      <<- 0
#       timesAccepted <<- 0
#       timesAdapted  <<- 0
#       #iterRun <<- 0
#       storeParticleLP <<- -Inf
#       gamma1 <<- 0
#     }
#   )
# )



