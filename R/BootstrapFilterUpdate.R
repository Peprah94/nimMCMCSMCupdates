##  Contains code to run bootstrap particle filters.
##  We have a build function (buildBootstrapFilter),
##  and step function.
bootStepVirtualUpdate <- nimbleFunctionVirtual(
  run = function(m = integer(),
                 iterRun = integer(0),
                 storeModelValues = double(1),
                 threshNum=double(),
                 prevSamp = logical()) {
    returnType(double(1))
  },
  methods = list(
    returnESS = function() {
      returnType(double())
    }
  )
)

# Bootstrap filter as specified in Doucet & Johnasen '08,
# uses weights from previous time point to calculate likelihood estimate.
bootFStepUpdate <- nimbleFunction(
  name = 'bootFStepUpdate',
  contains = bootStepVirtualUpdate,
  setup = function(model,
                   mvEWSamples,
                   mvWSamples,
                   nodes,
                   iNode,
                   names,
                   saveAll,
                   smoothing,
                   resamplingMethod,
                   silent = FALSE,
                   iNodePrev,
                   latent,
                   target,
                   mvSamplesEst) {

    # setting up parameters
    notFirst <- iNode != 1
    last <- iNode == length(nodes)
    prevNode <- nodes[if(notFirst) iNode-1 else iNode]


      modelSteps <- particleFilter_splitModelSteps(model, nodes, iNode, notFirst)
      prevDeterm <- modelSteps$prevDeterm
     # print(prevDeterm)
      calc_thisNode_self <- modelSteps$calc_thisNode_self
      calc_thisNode_deps <- modelSteps$calc_thisNode_deps
      targetNodesAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)

      thisNode <- nodes[iNode]


      # set up occupancy models
      # ensure that values y>0 have z's = 1 and those with y=0 have z's = NA
      isAllData <- all(model$isBinary(latent) == TRUE)

        calc_thisNode_self1 <- calc_thisNode_self[!model$isData(calc_thisNode_self)]
        calc_thisNode_self2 <- calc_thisNode_self[model$isData(calc_thisNode_self)]
        calc_thisNode_self2Vals <- rep(1, length(calc_thisNode_self2))

    ## t is the current time point.
    t <- iNode
    ## Get names of xs node for current and previous time point (used in copy)
    if(saveAll == 1){
      allPrevNodes <- model$expandNodeNames(nodes[1:(iNode-1)])
      prevXName <- prevNode
      thisXName <- thisNode
      currInd <- t
      prevInd <- t-1
       if(isTRUE(smoothing)){
        currInd <- 1
        prevInd <- 1
      }
    } else {
      allPrevNodes <- names
      prevXName <- names
      thisXName <- names
      currInd <- 1
      prevInd <- 1
    }

    isLast <- (iNode == length(nodes))
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
  run = function(m = integer(),
                 iterRun = integer(0),
                 storeModelValues = double(1),
                 threshNum = double(),
                 prevSamp = logical()) {
    returnType(double(1))

    wts <- numeric(m, init=FALSE)
    ids <- integer(m, 0)
    llEst <- numeric(m, init=FALSE)
    out <- numeric(2, init=FALSE)

    ################
    # Updating the new time steps: t+1, t+2, ..., T
    ############
    #values(model, targetNodesAsScalar) <<- storeModelValues
    if(t > iNodePrev){

    for(i in 1:m) {
      if(notFirst) {
        if(smoothing == 1){
          nimCopy(mvEWSamples, mvWSamples, nodes = allPrevNodes,
               nodesTo = allPrevNodes, row = i, rowTo=i)
        }
        nimCopy(mvEWSamples, model, nodes = prevXName, nodesTo = prevNode, row = i)
        model$calculate(prevDeterm)
        }
      #}
      #isAllData <- all(model$isData(calc_thisNode_self) == TRUE)
      if(isAllData){
        #calc_thisNode_self1 <- calc_thisNode_self[!model$isData(calc_thisNode_self)]
        model$simulate(calc_thisNode_self1)
        values(model, calc_thisNode_self2) <<- calc_thisNode_self2Vals
       # print(values(model, calc_thisNode_self1))
      }else{
      model$simulate(calc_thisNode_self)
      }

      #model$simulate(calc_thisNode_deps)
      ## The logProbs of calc_thisNode_self are, correctly, not calculated.
      nimCopy(model, mvWSamples, nodes = thisNode, nodesTo = thisXName, row = i)

       #model$calculate()
      wts[i]  <- model$calculate(calc_thisNode_deps)
      #print(wts[i])
      if(is.nan(wts[i])){
        out[1] <- -Inf
        out[2] <- 0
        return(out)
      }
      if(prevSamp == 0){ ## defaults to 1 for first step and then comes from the previous step
        wts[i] <- wts[i] + mvWSamples['wts',i][prevInd] #mvWSamples['wts',i][currInd]
        llEst[i] <- wts[i]
      } else {
        llEst[i] <- wts[i] - log(m)
      }
    }

    maxllEst <- max(llEst)
    #print(maxllEst)
    stepllEst <- maxllEst + log(sum(exp(llEst - maxllEst)))

    if(is.nan(stepllEst)){
      out[1] <- -Inf
      out[2] <- 0
      #return(out)
    }else if(stepllEst == Inf | stepllEst == -Inf){
      out[1] <- -Inf
      out[2] <- 0
      #return(out)
    }else{
    out[1] <- stepllEst
}


    ## Normalize weights and calculate effective sample size, using log-sum-exp trick to avoid underflow.
    maxWt <- max(wts)
    wts <- exp(wts - maxWt)/sum(exp(wts - maxWt))
    ess <<- 1/sum(wts^2)

    ## Determine whether to resample by weights or not.
    if(ess < threshNum){
      if(defaultResamplerFlag){
        rankSample(wts, m, ids, silent)
      } else {
        ids <- resamplerFunctionList[[1]]$run(wts)
      }
      ## out[2] is an indicator of whether resampling takes place.
      ## Resampling affects how ll estimate is calculated at next time point.
      out[2] <- 1
      for(i in 1:m){
        nimCopy(mvWSamples, mvEWSamples, nodes = thisXName, nodesTo = thisXName,
             row = ids[i], rowTo = i)
        if(smoothing == 1){
          nimCopy(mvWSamples, mvEWSamples, nodes = allPrevNodes,
               nodesTo = allPrevNodes, row = ids[i], rowTo=i)
        }
        mvWSamples['wts',i][currInd] <<- log(wts[i])
      }
    } else {
      out[2] <- 0
      for(i in 1:m){
        nimCopy(mvWSamples, mvEWSamples, nodes = thisXName, nodesTo = thisXName,
             row = i, rowTo = i)
        if(smoothing == 1){
          copy(mvWSamples, mvEWSamples, nodes = allPrevNodes,
               nodesTo = allPrevNodes, row = i, rowTo=i)
        }
        mvWSamples['wts',i][currInd] <<- log(wts[i])
      }
    }
    for(i in 1:m){
      mvWSamples['bootLL',i][currInd] <<- out[1]
    }
    #print(stepllEst)
   #print(out)

    return(out)
    }else{
      ########################
      #    Copying information from reduced model
      #######################

      # Copy the target values into the model too
      #nimCopy(from = mvSamplesEst, to = model, nodes = target, row = iterRun, rowTo = 1)

      #update the deterministic vars
      # if(notFirst) {
      #   model$calculate()
      # }
      #copy from saved modelValues to model
      if(isAllData){#for occupancy example
        nimCopy(mvSamplesEst, model, nodes = calc_thisNode_self1, nodesTo = calc_thisNode_self1, row = iterRun, rowTo = 1)
        values(model, calc_thisNode_self2) <<- calc_thisNode_self2Vals
      }else{
        nimCopy(mvSamplesEst, model, nodes = thisNode, nodesTo = thisXName, row = iterRun, rowTo = 1)
      }

      for(i in 1:m) {
        # if(notFirst) {
        #   if(smoothing == 1){
        #     nimCopy(mvEWSamples, mvWSamples, nodes = allPrevNodes,
        #             nodesTo = allPrevNodes, row = i, rowTo=i)
        #   }
        #   nimCopy(mvSamplesEst, model, nodes = prevXName, nodesTo = prevNode, row = iterRun, rowTo = i)
        #   #model$calculate(prevDeterm)
        # }

        # copy results to mvSaved values
        nimCopy(mvSamplesEst, mvWSamples, nodes = thisNode, nodesTo = thisXName, row=iterRun, rowTo = i)
        nimCopy(mvWSamples, mvEWSamples, thisXName, thisXName, row = i,  rowTo = i)
      }

      wtsEst <- model$calculate(calc_thisNode_deps)
     # print(wtsEst)
      #print(wtsEst)
      for(i in 1:m){
        wts[i] <- wtsEst #model$calculate(calc_thisNode_deps)
        if(is.nan(wts[i])){
          out[1] <- -Inf
          out[2] <- 1
          return(out)
        }
        if(prevSamp == 0){ ## defaults to 1 for first step and then comes from the previous step
          wts[i] <- wts[i] + mvWSamples['wts',i][prevInd]
          llEst[i] <- wts[i]
        } else {
          llEst[i] <- wts[i] - log(m)
        }
        #llEst[i] <- wts[i] - log(m)
      }

      # The samples are assumed to be equally weighted. So out[2] = 1
      # calculate log-likelihood to be used by sampler
      maxllEst <- max(llEst)
      #print(llEst[1])
      #print(maxllEst)
      stepllEst <- maxllEst + log(sum(exp(llEst - maxllEst)))
      #print(stepllEst)
      if(is.nan(stepllEst)){
        out[1] <- -Inf
        out[2] <- 1
        #return(out)
      }else if(stepllEst == Inf | stepllEst == -Inf){
        out[1] <- -Inf
        out[2] <- 1
        #return(out)
      }else{
        out[1] <- stepllEst
      }

      # We are not re-sampling
      out[2] <- 1
      ## Normalize weights and calculate effective sample size, using log-sum-exp trick to avoid underflow.
      maxWt <- max(wts)
      wts <- exp(wts - maxWt)/sum(exp(wts - maxWt))
      ess <<- 1/sum(wts^2)
      for(i in 1:m){
        #copy(mvWSamples, mvEWSamples, thisXName, thisXName, row = i,  rowTo = i)
        mvWSamples['wts',i][currInd] <<- log(wts[i])
        if(smoothing == 1){
          copy(mvWSamples, mvEWSamples, nodes = allPrevNodes,
               nodesTo = allPrevNodes, row = i, rowTo=i)
        }
      }


      for(i in 1:m){
        mvWSamples['bootLL',i][currInd] <<- stepllEst
      }
#print(out)
      return(out)
  }
    },
  methods = list(
    returnESS = function(){
      returnType(double(0))
      return(ess)
    }
  )
)

#' Create an updated bootstrap particle filter algorithm to estimate log-likelihood.
#'
#'@description Create an updated bootstrap particle filter algorithm for a given NIMBLE state space model.
#'
#' @param model A nimble model object, typically representing a state
#'  space model or a hidden Markov model.
#' @param nodes  A character vector specifying the latent model nodes
#'  over which the particle filter will stochastically integrate to
#'  estimate the log-likelihood function.  All provided nodes must be stochastic.
#'  Can be one of three forms: a variable name, in which case all elements in the variable
#'  are taken to be latent (e.g., 'x'); an indexed variable, in which case all indexed elements are taken
#'  to be latent (e.g., 'x[1:100]' or 'x[1:100, 1:2]'); or a vector of multiple nodes, one per time point,
#'  in increasing time order (e.g., c("x[1:2, 1]", "x[1:2, 2]", "x[1:2, 3]", "x[1:2, 4]")).
#' @param mvSamplesEst  A modelValue object contained posterior samples from the reduced model using MCMC.
#' @param control  A list specifying different control options for the particle filter.  Options are described in the details section below.
#' @author Kwaku Peprah Adjei
#' @details
#' Each of the \code{control()} list options are described in detail here:
#' \describe{
#'  \item{thresh}{ A number between 0 and 1 specifying when to resample: the resampling step will occur when the
#'   effective sample size is less than \code{thresh} times the number of particles. Defaults to 0.8. Note that at the last time step, resampling will always occur so that the \code{mvEWsamples} \code{modelValues} contains equally-weighted samples.}
#'  \item{resamplingMethod}{The type of resampling algorithm to be used within the particle filter. Can choose between \code{'default'} (which uses NIMBLE's \code{rankSample()} function),  \code{'systematic'}, \code{'stratified'}, \code{'residual'}, and \code{'multinomial'}.  Defaults to \code{'default'}.  Resampling methods other than \code{'default'} are currently experimental.}
#'  \item{saveAll}{Indicates whether to save state samples for all time points (TRUE), or only for the most recent time point (FALSE)}
#'  \item{smoothing}{Decides whether to save smoothed estimates of latent states, i.e., samples from f(x[1:t]|y[1:t]) if \code{smoothing = TRUE}, or instead to save filtered samples from f(x[t]|y[1:t]) if \code{smoothing = FALSE}.  \code{smoothing = TRUE} only works if \code{saveAll = TRUE}.}
#'  \item{timeIndex}{An integer used to manually specify which dimension of the latent state variable indexes time.
#'  Only needs to be set if the number of time points is less than or equal to the size of the latent state at each time point.}
#'  \item{initModel}{A logical value indicating whether to initialize the model before running the filtering algorithm.  Defaults to TRUE.}
#' }
#'  \item{iNodePrev}{An integer specifying the number of years used to fit the reduced model.}
#' }
#'
#'  The updated bootstrap filter starts by copying saved MCMC results into a
#'  modelValue object, and calculate weights for times t<= iNodePrev. For time iNode>iNodePrev,
#'  the updated bootstrap filter generates a sample of estimates from the
#'  prior distribution of the latent states of a state space model.  At each of these time point, these particles are propagated forward
#'  by the model's transition equation.  Each particle is then given a weight
#'  proportional to the value of the observation equation given that particle.
#'  The weights are used to draw an equally-weighted sample of the particles at this time point.
#'  The algorithm then proceeds
#'  to the next time point.  Neither the transition nor the observation equations are required to
#'  be normal for the bootstrap filter to work.
#'
#'  The resulting specialized particle filter algorithm will accept a
#'  single integer argument (\code{m}, default 10,000), which specifies the number
#'  of random \'particles\' to use for estimating the log-likelihood.  The algorithm
#'  returns the estimated log-likelihood value, and saves
#'  unequally weighted samples from the posterior distribution of the latent
#'  states in the \code{mvWSamples} modelValues object, with corresponding logged weights in \code{mvWSamples['wts',]}.
#'  An equally weighted sample from the posterior can be found in the \code{mvEWSamples} \code{modelValues} object.
#'
#'  Note that if the \code{thresh} argument is set to a value less than 1, resampling may not take place at every time point.
#'  At time points for which resampling did not take place, \code{mvEWSamples} will not contain equally weighted samples.
#'  To ensure equally weighted samples in the case that \code{thresh < 1}, we recommend resampling from \code{mvWSamples} at each time point
#'  after the filter has been run, rather than using \code{mvEWSamples}.
#'
#' @section \code{returnESS()} Method:
#'  Calling the \code{returnESS()} method of a bootstrap filter after that filter has been \code{run()} for a given model will return a vector of ESS (effective
#'  sample size) values, one value for each time point.
#'
#' @export
#'
#' @family smc update methods
#' @references Michaud, N., de Valpine, P., Turek, D., Paciorek, C. J., & Nguyen, D. (2021). Sequential Monte Carlo methods in the nimble and nimbleSMC R packages. \emph{Journal of Statistical Software}. 100, 1-39.
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
#' my_BootF <- buildBootstrapFilterUpdate(estimationModel,
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

buildBootstrapFilterUpdate <- nimbleFunction(
  name = 'buildBootstrapFilterUpdate',
  setup = function(model,
                   nodes,
                   mvSamplesEst = list(),
                   target,
                   control = list()) {

    #control list extraction
    thresh <- control[['thresh']]
    saveAll <- control[['saveAll']]
    smoothing <- control[['smoothing']]
    silent <- control[['silent']]
    timeIndex <- control[['timeIndex']]
    initModel <- control[['initModel']]
    iNodePrev <- control[['iNodePrev']] #how many time points were used to fit the reduced model
    #initModel <- control[['initModel']]
    #M <- control[['M']]
    resamplingMethod <- control[['resamplingMethod']]
    if(is.null(thresh)) thresh <- 0.8
    if(is.null(silent)) silent <- TRUE
    if(is.null(saveAll)) saveAll <- FALSE
    if(is.null(smoothing)) smoothing <- FALSE
    if(is.null(initModel)) initModel <- TRUE
    if(is.null(resamplingMethod)) resamplingMethod <- 'default'
    if(!(resamplingMethod %in% c('default', 'multinomial', 'systematic', 'stratified',
                                 'residual')))
      stop('buildBootstrapFilter: "resamplingMethod" must be one of: "default", "multinomial", "systematic", "stratified", or "residual". ')
    ## latent state info
    targetNodesAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
  #storeModelValues <- values(model, targetNodesAsScalar)
    nodes <- findLatentNodes(model, nodes, timeIndex)
    dims <- lapply(nodes, function(n) nimDim(model[[n]]))
    if(length(unique(dims)) > 1)
      stop('buildBootstrapFilter: sizes or dimensions of latent states varies.')
    vars <- model$getVarNames(nodes =  nodes)

    my_initializeModel <- initializeModel(model, silent = silent)

    if(0 > thresh || 1 < thresh || !is.numeric(thresh))
      stop('buildBootstrapFilter: "thresh" must be between 0 and 1.')
    if(!saveAll & smoothing)
      stop("buildBootstrapFilter: must have 'saveAll = TRUE' for smoothing to work.")

    if(!is.numeric(iNodePrev) || iNodePrev < 0)
      stop("buildBootstrapFilter: must have 'iNodePrev' numeric and greater than 0")
    ## Create mv variables for x state and sampled x states.  If saveAll=TRUE,
    ## the sampled x states will be recorded at each time point.
    modelSymbolObjects <- model$getSymbolTable()$getSymbolObjects()[vars]

    if(saveAll){
      names <- sapply(modelSymbolObjects, function(x)return(x$name))
      type <- sapply(modelSymbolObjects, function(x)return(x$type))
      size <- lapply(modelSymbolObjects, function(x)return(x$size))

      mvEWSamples <- modelValues(modelValuesConf(vars = names,
                                                 types = type,
                                                 sizes = size))
      latent <- names
      names <- c(names, "wts", "bootLL")
      type <- c(type, "double", "double")
      size$wts <- length(dims)
      size$bootLL <- length(dims)
      ##  Only need one weight per particle (at time T) if smoothing == TRUE.
      if(smoothing){
      size$wts <- 1
      size$bootLL <- 1
      }
      mvWSamples  <- modelValues(modelValuesConf(vars = names,
                                                 types = type,
                                                 sizes = size))

    }else{
      names <- sapply(modelSymbolObjects, function(x)return(x$name))
      type <- sapply(modelSymbolObjects, function(x)return(x$type))
      size <- lapply(modelSymbolObjects, function(x)return(x$size))
      size[[1]] <- as.numeric(dims[[1]])

      mvEWSamples <- modelValues(modelValuesConf(vars = names,
                                                 types = type,
                                                 sizes = size))
      latent <- names
      names <- c(names, "wts", "bootLL")
      type <- c(type, "double", "double")
      size$wts <- 1
      size$bootLL <- 1
      mvWSamples  <- modelValues(modelValuesConf(vars = names,
                                                 types = type,
                                                 sizes = size))
      names <- names[1]
    }


    bootStepFunctions <- nimbleFunctionList(bootStepVirtualUpdate)
    for(iNode in seq_along(nodes)){
      bootStepFunctions[[iNode]] <- bootFStepUpdate(model, mvEWSamples, mvWSamples,
                                                    nodes, iNode, names, saveAll,
                                                    smoothing, resamplingMethod,
                                                    silent,
                                                    iNodePrev, latent,
                                                    target, mvSamplesEst)
    }

    essVals <- rep(0, length(nodes))
    lastLogLik <- -Inf
    runTime <- 1
  },
  run = function(m = integer(default = 10000),
                 iterRun = integer(default = 1),
                 storeModelValues = double(1)
                 ) {
    returnType(double())

    if(initModel) my_initializeModel$run()
    resize(mvWSamples, m)
    resize(mvEWSamples, m)
    threshNum <- ceiling(thresh*m)
    logL <- 0
    ## prevSamp indicates whether resampling took place at the previous time point.
    prevSamp <- 1
    for(iNode in seq_along(bootStepFunctions)) {
      if(iNode == length(bootStepFunctions)){
        threshNum <- m  ## always resample at last time step so mvEWsamples is equally-weighted
      }
        out <- bootStepFunctions[[iNode]]$run(m, iterRun, storeModelValues, threshNum, prevSamp)
      logL <- logL + out[1]
      prevSamp <- out[2]
      #print(iNode)
      essVals[iNode] <<- bootStepFunctions[[iNode]]$returnESS()
      if(logL == -Inf) {lastLogLik <<- logL; return(logL)}
      if(is.nan(logL)) {lastLogLik <<- -Inf; return(-Inf)}
      if(logL == Inf)  {lastLogLik <<- -Inf; return(-Inf)}
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
