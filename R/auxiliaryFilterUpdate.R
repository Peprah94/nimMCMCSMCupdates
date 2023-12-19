
##  Contains code to run auxiliary particle filters.
##  We have a build function (buildAuxiliaryFilter),
##  and step function (auxFStep)
##
##  This version of the APF is based on
##  Pitt et al., 2012

auxStepVirtual2 <- nimble::nimbleFunctionVirtual(
  run = function(m = integer(),iterRun = integer(0),storeModelValues = double(1)) {
    returnType(double())
  },
  methods = list(
    returnESS = function() {
      returnType(double())
    }
  )
)

auxFuncVirtualUpdate <- nimble::nimbleFunctionVirtual(
  methods = list(
    lookahead = function(){}
  )
)

# Loook ahead functions
auxLookFuncUpdate  = nimble::nimbleFunction(
  name = 'auxLookFuncUpdate',
  contains = auxFuncVirtualUpdate ,
  setup = function(model, node){
  },
  methods = list(
    lookahead = function(){
      model[[node]] <<- model$getParam(node, 'mean')
    })
)


auxSimFuncUpdate  = nimble::nimbleFunction(
  name = 'auxSimFuncUpdate',
  contains = auxFuncVirtualUpdate ,
  setup = function(model, node){},
  methods = list(
    lookahead = function(){
      model$simulate(node)
    })
)

# Function to run for each time step
auxFStepUpdate <- nimbleFunction(
  name = 'auxFStepUpdate',
  contains = auxStepVirtual2,
  setup = function(model, mvEWSamples, mvWSamples, nodes, iNode, names,
                   saveAll, smoothing, lookahead, resamplingMethod,
                   silent = TRUE,
                   iNodePrev, target,
                   latent,
                   mvSamplesEst) {
    notFirst <- iNode != 1
    last <- iNode == length(nodes)
    prevNode <- nodes[if(notFirst) iNode-1 else iNode]

    modelSteps <- particleFilter_splitModelSteps(model, nodes, iNode, notFirst)
    prevDeterm <- modelSteps$prevDeterm
    calc_thisNode_self <- modelSteps$calc_thisNode_self
    calc_thisNode_deps <- modelSteps$calc_thisNode_deps
    targetNodesAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    ## prevDeterm <- model$getDependencies(prevNode, determOnly = TRUE)
    thisNode <- nodes[iNode]

    isAllData <- all(model$isBinary(latent) == TRUE)
    #if(isAllData){
    calc_thisNode_self1 <- calc_thisNode_self[!model$isData(calc_thisNode_self)]
    calc_thisNode_self2 <- calc_thisNode_self[model$isData(calc_thisNode_self)]
    calc_thisNode_self2Vals <- rep(1, length(calc_thisNode_self2))
    ## thisDeterm <- model$getDependencies(thisNode, determOnly = TRUE)
    ## thisData   <- model$getDependencies(thisNode, dataOnly = TRUE)
    ## t is the current time point.
    t <- iNode
    ## Get names of x and xs node for current and previous time point,
    ## will be different depending on whether we are saving all time points
    ## or only the most recent.
    if(saveAll == 1){
      allPrevNodes <- model$expandNodeNames(nodes[1:(iNode-1)])
      prevXName <- prevNode
      thisXName <- thisNode
      currInd <- t
      prevInd <- t-1
      if(smoothing == TRUE){
        currInd <- 1
        prevInd <- 1
      }
    }
    else{
      allPrevNodes <- names # dummy value -- not used.
      prevXName <- names
      thisXName <- names
      currInd <- 1
      prevInd <- 1
    }
#iterRun <- 1
    auxFuncList <- nimbleFunctionList(auxFuncVirtualUpdate)
    allLatentNodes <- model$expandNodeNames(calc_thisNode_self, sort = TRUE) ## They should already be sorted, but sorting here is a failsafe.
    numLatentNodes <- length(allLatentNodes)
    if(lookahead == "mean"){
      for(i in 1:numLatentNodes)
        auxFuncList[[i]] <- auxLookFuncUpdate(model, allLatentNodes[i])
    }else{
      for(i in 1:numLatentNodes)
        auxFuncList[[i]] <- auxSimFuncUpdate(model,  allLatentNodes)
    }
    ess <- 0
    resamplerFunctionList <- nimbleFunctionList(resamplerVirtual)
    defaultResamplerFlag <- FALSE
    if(resamplingMethod == 'default'){
      resamplerFunctionList[[1]] <- residualResampleFunction()
      defaultResamplerFlag <- TRUE
    }
    if(resamplingMethod == 'residual')
      resamplerFunctionList[[1]] <- residualResampleFunction()
    if(resamplingMethod == 'multinomial')
      resamplerFunctionList[[1]] <- multinomialResampleFunction()
    if(resamplingMethod == 'stratified')
      resamplerFunctionList[[1]] <- stratifiedResampleFunction()
    if(resamplingMethod == 'systematic')
      resamplerFunctionList[[1]] <- systematicResampleFunction()
  },
  run = function(m = integer(),iterRun = integer(0),storeModelValues = double(1)) {
    returnType(double())
    auxll <- numeric(m, init=FALSE)
    auxWts <- numeric(m, init=FALSE)
    wts <- numeric(m, init=FALSE)
    ids <- integer(m, 0)
    ll <- numeric(m, init=FALSE)
    #wts <- ids <- auxWts <- auxll <- ll <- numeric(m)
    #out <- numeric(2)
    ## This is the look-ahead step, not conducted for first time-point.
    if(t > iNodePrev){
   values(model, targetNodesAsScalar) <<- storeModelValues
    if(notFirst){
      for(i in 1:m) {
        if(smoothing == 1){
          ## smoothing is only allowed if saveAll is TRUE, so this should be ok.
          ## i.e., mvEWSamples have been resampled.
          nimCopy(mvEWSamples, mvWSamples, nodes = allPrevNodes,
               nodesTo = allPrevNodes, row = i, rowTo=i)
        }
        nimCopy(mvWSamples, model, prevXName, prevNode, row=i)
        model$calculate(prevDeterm)
        ## The lookahead steps may include determ and stoch steps.
        if(lookahead == "mean"){
          for(j in 1:numLatentNodes)
            auxFuncList[[j]]$lookahead()
        } else auxFuncList[[1]]$lookahead()

        ## Get p(y_t+1 | x_t+1).
        auxll[i] <- model$calculate(calc_thisNode_deps)
        if(is.nan(auxll[i])){
          return(-Inf)
        }
        ## Multiply (on log scale) by weight from time t.
        auxWts[i] <- auxll[i] + mvWSamples['wts',i][prevInd]
      }
      ## Normalize weights and resample, using log-sum-exp trick to avoid underflow.
      maxWt <- max(auxWts)
      normAuxWts <- exp(auxWts - maxWt)/sum(exp(auxWts - maxWt))
      if(defaultResamplerFlag == TRUE){
        rankSample(normAuxWts, m, ids, silent)
      }
      else{
        ids <- resamplerFunctionList[[1]]$run(normAuxWts)
      }
    }

    for(i in 1:m) {
      if(notFirst) {
        nimCopy(mvWSamples, model, nodes = prevXName, nodesTo = prevNode,
                row = ids[i])
        model$calculate(prevDeterm)
      }

      # Simulate from x_t+1 | x_t.
      # treat z as data for y>1 and simuate for y = 0
      if(isAllData){
        #calc_thisNode_self1 <- calc_thisNode_self[!model$isData(calc_thisNode_self)]
        model$simulate(calc_thisNode_self1)
        values(model, calc_thisNode_self2) <<- calc_thisNode_self2Vals
        #print(values(model, calc_thisNode_self2))
      }else{
        model$simulate(calc_thisNode_self)
      }
      #model$simulate(calc_thisNode_self)
      nimCopy(model, mvEWSamples, nodes = thisNode, nodesTo = thisXName, row=i)
      ## Get p(y_t+1 | x_t+1).
      ll[i] <- model$calculate(calc_thisNode_deps)
      if(is.nan(ll[i])){
        return(-Inf)
      }
      if(notFirst){
        ## Construct weight following step 4 of paper.
        wts[i] <- ll[i]-auxll[ids[i]]
      }
      else{
        ## First step has no auxiliary weights.
        wts[i] <- ll[i]
      }
    }
    ## Use log-sum-exp trick to avoid underflow.
    maxWt <- max(wts)
    normWts <- exp(wts - maxWt)/sum(exp(wts - maxWt))
    ess <<- 1/sum(normWts^2)
    for(i in 1:m){
      ## Save weights for use in next timepoint's look-ahead step.
      mvWSamples['wts', i][currInd] <<- log(normWts[i])
    }
    if(defaultResamplerFlag == TRUE){
      rankSample(normWts, m, ids, silent)
    }
    else{
      ids <- resamplerFunctionList[[1]]$run(normWts)
    }
    for(i in 1:m){
      copy(mvEWSamples, mvWSamples, thisXName, thisXName, row = i,  rowTo = i)
    }

    if(saveAll | last) {
      for(i in 1:m) {
        if(smoothing == 1){
          copy(mvWSamples, mvEWSamples, nodes = allPrevNodes,
               nodesTo = allPrevNodes, row = ids[i], rowTo=i)
        }

        copy(mvWSamples, mvEWSamples, thisXName, thisXName, ids[i], i)
      }
    }

    ## Calculate likelihood p(y_t+1 | y_1:t) as in equation (3) of paper.
    ## Use log-sum-exp trick to avoid underflow.
    if(notFirst){
      maxWt <- max(wts)
      maxAuxWt <- max(auxWts)
      outLL <- log(sum(exp(wts - maxWt))) + maxWt - log(m) + log(sum(exp(auxWts - maxAuxWt))) + maxAuxWt
    } else {
      maxWt <- max(wts)
      outLL <- log(sum(exp(wts - maxWt))) + maxWt - log(m)
    }
    for(i in 1:m){
      mvWSamples['auxlog',i][currInd] <<- outLL
    }

    return(outLL)
     }else{
# for t < iNodePrev


       nimCopy(from = mvSamplesEst, to = model, nodes = target, row = iterRun, rowTo = 1)
        # if(notFirst){
        #   for(i in 1:m) {
        #    if(smoothing == 1){
        #      ## smoothing is only allowed if saveAll is TRUE, so this should be ok.
        #      ## i.e., mvEWSamples have been resampled.
        #      nimCopy(mvEWSamples, mvWSamples, nodes = allPrevNodes,
        #              nodesTo = allPrevNodes, row = i, rowTo=i)
        #    }
        #    nimCopy(mvWSamples, model, prevXName, prevNode, row=i)
        #    model$calculate(prevDeterm)
        #    ## The lookahead steps may include determ and stoch steps.
        #    if(lookahead == "mean"){
        #      for(j in 1:numLatentNodes)
        #        auxFuncList[[j]]$lookahead()
        #    } else auxFuncList[[1]]$lookahead()
        #
        #    ## Get p(y_t+1 | x_t+1).
        #    auxll[i] <- model$calculate(calc_thisNode_deps)
        #    if(is.nan(auxll[i])){
        #      return(-Inf)
        #    }
        #    ## Multiply (on log scale) by weight from time t.
        #    if(auxll[i] == -Inf | auxll[i] == Inf){
        #      auxWts[i] <- mvWSamples['wts',i][prevInd]
        #      auxll[i] <- 0
        #    }else{
        #    auxWts[i] <- auxll[i] + mvWSamples['wts',i][prevInd]
        #    }
        #  }
        # }

       if(isAllData){
         #calc_thisNode_self1 <- calc_thisNode_self[!model$isData(calc_thisNode_self)]
         #model$simulate(calc_thisNode_self1)
         nimCopy(mvSamplesEst, model, nodes = calc_thisNode_self1, nodesTo = calc_thisNode_self1, row = iterRun, rowTo = 1)
         values(model, calc_thisNode_self2) <<- calc_thisNode_self2Vals
         #print(values(model, calc_thisNode_self2))
       }else{
         nimCopy(mvSamplesEst, model, nodes = calc_thisNode_self, nodesTo = calc_thisNode_self, row = iterRun, rowTo = 1)
       }
       nimCopy(mvSamplesEst, model, nodes = thisNode, nodesTo = thisXName, row = iterRun, rowTo = 1)
       #
       # for(i in 1:m) {
       #   if(notFirst) {
       #     nimCopy(mvWSamples, model, nodes = prevXName, nodesTo = prevNode,
       #             row = i)
       #     model$calculate(prevDeterm)
       #   }

         # Simulate from x_t+1 | x_t.
         # treat z as data for y>1 and simuate for y = 0
       for(i in 1:m){
         nimCopy(mvSamplesEst, mvWSamples, nodes = thisNode, nodesTo = thisNode, row=iterRun, rowTo = i)
         nimCopy(mvWSamples, mvEWSamples, thisXName, thisXName, row = i,  rowTo = i)
         ## Get p(y_t+1 | x_t+1).
       }

         llEst <- model$calculate(calc_thisNode_deps)
         for(i in 1:m){
           ll[i] <- llEst
         if(is.nan(ll[i])){
           return(-Inf)
         }
         if(notFirst){
           ## Construct weight following step 4 of paper.
           wts[i] <- ll[i] #-auxll[i]
         }
         else{
           ## First step has no auxiliary weights.
           wts[i] <- ll[i]
         }

         #wts[i] <- 1/m
         }

       ## Use log-sum-exp trick to avoid underflow.
       maxWt <- max(wts)
       normWts <- exp(wts - maxWt)/sum(exp(wts - maxWt))
       ##normalised weights would be different based on the auxiliary step
       ess <<- 1/sum(normWts^2)
       for(i in 1:m){
         ## Save weights for use in next timepoint's look-ahead step.
         mvWSamples['wts', i][currInd] <<- log(normWts[i])
         #copy(mvWSamples, mvEWSamples, thisXName, thisXName, row = i, rowTo = i)
       }


       if(saveAll | last) {
         for(i in 1:m) {
           if(smoothing == 1){
             copy(mvWSamples, mvEWSamples, nodes = allPrevNodes,
                  nodesTo = allPrevNodes, row = i, rowTo=i)
           }

           copy(mvWSamples, mvEWSamples, thisXName, thisXName, i, i)
         }
       }
       ## Calculate likelihood p(y_t+1 | y_1:t) as in equation (3) of paper.
       ## Use log-sum-exp trick to avoid underflow.
         # maxWt <- max(wts)
         # outLL <- log(sum(exp(wts - maxWt))) + maxWt - log(m)
       if(notFirst){
         maxWt <- max(wts)
         maxAuxWt <- max(auxWts)
         outLL <- log(sum(exp(wts - maxWt))) + maxWt - log(m) + log(sum(exp(auxWts - maxAuxWt))) + maxAuxWt
       } else {
         maxWt <- max(wts)
         outLL <- log(sum(exp(wts - maxWt))) + maxWt - log(m)
       }

       for(i in 1:m){
         mvWSamples['auxlog',i][currInd] <<- outLL
       }
       return(outLL)
     }
  },
  methods = list(
    returnESS = function() {
      returnType(double(0))
      return(ess)
    }
  )
)

#' Create an updated auxiliary particle filter algorithm to estimate log-likelihood.
#'
#' @description Create an updated auxiliary particle filter algorithm for a given NIMBLE state space model.
#'
#' @param model A NIMBLE model object, typically representing a state space model or a hidden Markov model.
#' @param nodes  A character vector specifying the latent model nodes
#'  over which the particle filter will stochastically integrate to
#'  estimate the log-likelihood function.  All provided nodes must be stochastic.
#'  Can be one of three forms: a variable name, in which case all elements in the variable
#'  are taken to be latent (e.g., 'x'); an indexed variable, in which case all indexed elements are taken
#'  to be latent (e.g., 'x[1:100]' or 'x[1:100, 1:2]'); or a vector of multiple nodes, one per time point,
#'  in increasing time order (e.g., c("x[1:2, 1]", "x[1:2, 2]", "x[1:2, 3]", "x[1:2, 4]")).
#' @param mvSamplesEst  A modelValue object contained posterior samples from the reduced model using MCMC.
#' @param control  A list specifying different control options for the particle filter.  Options are described in the details section below.
#' @author  Kwaku Peprah Adjei
#' @details
#' Each of the \code{control()} list options are described in detail here:
#' \describe{
#' \item{lookahead}{The lookahead function used to calculate auxiliary weights.  Can choose between \code{'mean'} and \code{'simulate'}.
#'  Defaults to \code{'simulate'}.}
#'  \item{resamplingMethod}{The type of resampling algorithm to be used within the particle filter.  Can choose between \code{'default'} (which uses NIMBLE's \code{rankSample()} function), \code{'systematic'}, \code{'stratified'}, \code{'residual'}, and \code{'multinomial'}.  Defaults to \code{'default'}. Resampling methods other than \code{'default'} are currently experimental.}
#'  \item{saveAll}{Indicates whether to save state samples for all time points (\code{TRUE}), or only for the most recent time point (\code{FALSE})}
#'  \item{smoothing}{Decides whether to save smoothed estimates of latent states, i.e., samples from f(x[1:t]|y[1:t]) if \code{smoothing = TRUE}, or instead to save filtered samples from f(x[t]|y[1:t]) if \code{smoothing = FALSE}. \code{smoothing = TRUE} only works if \code{saveAll = TRUE}.}
#'  \item{timeIndex}{An integer used to manually specify which dimension of the latent state variable indexes time. This need only be set if the number of time points is less than or equal to the size of the latent state at each time point.}
#'  \item{initModel}{A logical value indicating whether to initialize the model before running the filtering algorithm.  Defaults to TRUE.}
#'  }
#'  \item{iNodePrev}{An integer specifying the number of years used to fit the reduced model.}
#' }
#'
#' The auxiliary particle filter modifies the bootstrap filter (\code{\link{buildBootstrapFilter}})
#' by adding a lookahead step to the algorithm: before propagating particles from one time
#' point to the next via the transition equation, the auxiliary filter calculates a weight
#' for each pre-propogated particle by predicting how well the particle will agree with the
#' next data point.  These pre-weights are used to conduct an initial resampling step before
#' propagation.
#'
#'  The resulting specialized particle filter algorthm will accept a
#'  single integer argument (\code{m}, default 10,000), which specifies the number
#'  of random \'particles\' to use for estimating the log-likelihood.  The algorithm
#'  returns the estimated log-likelihood value, and saves
#'  unequally weighted samples from the posterior distribution of the latent
#'  states in the \code{mvWSamples} modelValues object, with corresponding logged weights in \code{mvWSamples['wts',]}.
#'  An equally weighted sample from the posterior can be found in the \code{mvEWsamp} modelValues object.
#'
#'   The auxiliary particle filter uses a lookahead function to select promising particles before propagation.  This function can eithre be the expected
#'   value of the latent state at the next time point (\code{lookahead = 'mean'}) or a simulation from the distribution of the latent state at the next time point (\code{lookahead = 'simulate'}), conditioned on the current particle.
#'
#'  @section \code{returnESS()} Method:
#'  Calling the \code{returnESS()} method of an auxiliary particle filter after that filter has been \code{run()} for a given model will return a vector of ESS (effective
#'  sample size) values, one value for each time point.

#' @export
#'
#' @family smc update methods
#' @references Michaud, N., de Valpine, P., Turek, D., Paciorek, C. J., & Nguyen, D. (2021). Sequential Monte Carlo methods in the nimble and nimbleSMC R packages. \emph{Journal of Statistical Software}. 100, 1-39.
#'
#' @examples
#' ## For illustration only.
#' stateSpaceCode <- nimbleCode({
#'   x0 ~ dnorm(0, var = 1)
#'   x[1] ~ dnorm(a* x0, var = 1)
#'   y[1] ~ dnorm(x[1]*c, var = .5)
#'   for(t in 2:10){
#'     x[t] ~ dnorm(a * x[t-1], var = 1)
#'     y[t] ~ dnorm(x[t]*c, var = .5)
#'   }
#' })
#'
#' model <- nimbleModel(code = exampleCode, data = list(y = rnorm(10)),
#'                      inits = list(x0 = 0, x = rnorm(10)))
#' my_BootF <- buildAuxiliaryFilterUpdate(estimationModel,
#' latent,
#' mvSamplesEst = mvSamplesEst,
#' target = target,
#' control = pfControl)
#' ## Now compile and run, e.g.,
#' ## targetAsScalar <- estimationModel$expandNodeNames(target, returnScalarComponents = TRUE)
#' ## compiledParticleFilter <- compileNimble(estimationModel,  particleFilterEst)
#' ## logLik <- compiledParticleFilter$particleFilterEst$run(m = 2000, iterRun = 1, storeModelValues = values(estimationModel, targetAsScalar))
#' ## ESS <- compiledParticleFilter$particleFilterEst$returnESS()
#' ## boot_X <- as.matrix(compiledParticleFilter$particleFilterEst$mvEWSamples, 'x')

buildAuxiliaryFilterUpdate <- nimbleFunction(
  name = 'buildAuxiliaryFilterUpdate',
  #contains = auxStepVirtual1Update,
  setup = function(model, nodes, mvSamplesEst, target, control = list()) {

    ## Control list extraction.
    saveAll <- control[['saveAll']]
    smoothing <- control[['smoothing']]
    silent <- control[['silent']]
    timeIndex <- control[['timeIndex']]
    lookahead <- control[['lookahead']]
    iNodePrev <- control[['iNodePrev']]
    initModel <- control[['initModel']]
    M <- control[['M']]
    latent <- nodes
    resamplingMethod <- control[['resamplingMethod']]
    if(is.null(silent)) silent <- TRUE
    if(is.null(saveAll)) saveAll <- FALSE
    if(is.null(smoothing)) smoothing <- FALSE
    if(is.null(lookahead)) lookahead <- 'simulate'
    if(is.null(initModel)) initModel <- FALSE
    if(!saveAll & smoothing) stop("must have saveAll = TRUE for smoothing to
                                  work")
    if(lookahead == "mean"){
      errors <- sapply(model$expandNodeNames(nodes), function(node){
        tryCatch(getParam(model, node, 'mean'), error=function(a){
          return("error")})})
      if(any(errors == "error", na.rm=TRUE))
        stop("cannot use 'mean' lookahead for this model, try 'simulate'")
    }
    else if(lookahead != "simulate"){
      stop("lookahead argument must be either 'simulate' or 'mean'")
    }

    if(is.null(resamplingMethod)) resamplingMethod <- 'default'
    if(!(resamplingMethod %in% c('default', 'multinomial', 'systematic', 'stratified',
                                 'residual')))
      stop('resamplingMethod must be one of: "default", "multinomial", "systematic",
           "stratified", or "residual". ')

    # target information
    #targetNodesAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)

    ## Latent state info.
    nodes <- findLatentNodes(model, nodes, timeIndex)

    dims <- lapply(nodes, function(n) nimDim(model[[n]]))
    if(length(unique(dims)) > 1)
      stop('sizes or dimensions of latent states varies')
    vars <- model$getVarNames(nodes =  nodes)

    my_initializeModel <- initializeModel(model, silent = silent)


    # Create mv variables for x state and sampled x states.  If saveAll=TRUE,
    # the sampled x states will be recorded at each time point.
    modelSymbolObjects <- model$getSymbolTable()$getSymbolObjects()[vars]


    if(saveAll){
      #
      names <- sapply(modelSymbolObjects, function(x)return(x$name))
      type <- sapply(modelSymbolObjects, function(x)return(x$type))
      size <- lapply(modelSymbolObjects, function(x)return(x$size))
      mvEWSamples <- modelValues(modelValuesConf(vars = names,
                                                 types = type,
                                                 sizes = size))

      #latent <- names
      names <- c(names, "wts", "auxlog")
      type <- c(type, "double", "double")
      size$wts <- length(dims)
      size$auxlog <- length(dims)
      ## Only need one weight per particle (at time T) if smoothing == TRUE.
      if(smoothing){
        size$wts <- 1
        size$auxlog <- 1
      }
      mvWSamples  <- modelValues(modelValuesConf(vars = names,
                                                 types = type,
                                                 sizes = size))

    }
    else{
      names <- sapply(modelSymbolObjects, function(x)return(x$name))
      type <- sapply(modelSymbolObjects, function(x)return(x$type))
      size <- lapply(modelSymbolObjects, function(x)return(x$size))
      size[[1]] <- as.numeric(dims[[1]])

      mvEWSamples <- modelValues(modelValuesConf(vars = names,
                                                 types = type,
                                                 sizes = size))
      #latent <- names
      names <- c(names, "wts", "auxlog")
      type <- c(type, "double", "double")
      size$wts <- 1
      size$auxlog <- 1
      mvWSamples  <- modelValues(modelValuesConf(vars = names,
                                                 types = type,
                                                 sizes = size))
    }

    names <- names[1]
    auxStepFunctionsUpdate <- nimbleFunctionList(auxStepVirtual2)
    for(iNode in seq_along(nodes))
      auxStepFunctionsUpdate[[iNode]] <- auxFStepUpdate(model, mvEWSamples, mvWSamples, nodes, iNode, names,
                                                   saveAll, smoothing, lookahead, resamplingMethod,
                                                   silent,
                                                   iNodePrev,
                                                   target,
                                                   latent,
                                                   mvSamplesEst)

    essVals <- rep(0, length(nodes))
    lastLogLik <- -Inf
    runTime <- 1

    #iNodePrev <- iNodePrev - 1
    #index <-
  },
  run = function(m = integer(default = 10000),
                 iterRun = integer(default = 1),
                 storeModelValues = double(1)) {
    returnType(double())
    declare(logL, double())
    if(initModel == TRUE) my_initializeModel$run()
    resize(mvEWSamples, m)
    resize(mvWSamples, m)
    logL <- 0
    for(iNode in seq_along(auxStepFunctionsUpdate)){
      logL <- logL + auxStepFunctionsUpdate[[iNode]]$run(m, iterRun, storeModelValues)
      essVals[iNode] <<- auxStepFunctionsUpdate[[iNode]]$returnESS()

      ## When all particles have 0 weight, likelihood becomes NAN
      ## this happens if top-level params have bad values - possible
      ## during pmcmc for example.
      if(is.nan(logL)) {lastLogLik <<- -Inf; return(-Inf)}
      if(logL == -Inf) {lastLogLik <<- logL; return(logL)}
      if(logL == Inf) {lastLogLik <<- -Inf; return(-Inf)}
    }
    lastLogLik <<- logL
    runTime <<- iterRun
    return(logL)
  },
  methods = list(
    getLastLogLik = function() {
      return(lastLogLik)
      returnType(double())
    },
    setLastLogLik = function(lll = double()) {
      lastLogLik <<- lll
},
getLastTimeRan = function() {
  return(runTime)
  returnType(integer())
},
setLastTimeRan = function(lll = integer()) {
  runTime <<- lll + 1
},
  returnESS = function(){
       returnType(double(1))
       return(essVals)
     }
  )
)



