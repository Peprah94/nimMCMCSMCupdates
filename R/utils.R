#' @import nimble
#' @import methods

particleFilter_splitModelSteps <- function(model,
                                           nodes,
                                           iNode,
                                           notFirst) {
  ## if !notFirst (it is the first node), prevNode will not be used anyway, so
  ## giving its value is a dummy.
  prevNode <- nodes[if(notFirst) iNode-1 else iNode]
  thisNode <- nodes[iNode]
  thisNodeExpanded <- model$expandNodeNames(thisNode, sort = TRUE)

  ## Set up steps for calculations internal to thisNode and for downstream simulation from thisNode
  thisDeterm <- model$getDependencies(thisNode, determOnly = TRUE)
  if(length(thisDeterm) > 0) {
    thisDeterm_is_intermediate <- logical(length(thisDeterm))
    for(i in seq_along(thisDeterm)) {
      theseDeps <- model$getDependencies(thisDeterm[i], stochOnly = TRUE)
      thisDeterm_is_intermediate[i] <- any(theseDeps %in% thisNodeExpanded)
    }
    thisDeterm_self <- thisDeterm[ thisDeterm_is_intermediate ]
    thisDeterm <- thisDeterm[ !thisDeterm_is_intermediate ]
    calc_thisNode_self <-  model$expandNodeNames(c(thisNodeExpanded, thisDeterm_self), ## only for the sort
                                                 sort = TRUE)

  } else {
    calc_thisNode_self <- thisNodeExpanded
  }
  ## calc_thisNode_deps <- all_deps_from_thisNode[-(1:finalSelfIndex)] ## some of these could be unnecessary...

  ## Roughly we have prevNode -> prevDeterm -> thisNode -> thisDeterm -> thisData
  ## However, there is further sorting to do.
  ## prevDeterm may have (i) components not really needed for thisNode (leading to other parts of the graph).
  ##                     (ii) components that also depend on one part of thisNode, on which other parts of thisNode might depend.
  ## thisNode might have internal ordering and might comprise multiple nodes.
  ## thisDeterm might have (i) components not really needed for thisData (leading to other parts of the graph).
  ##                       (ii) components that also depend on one part of thisData, on which other parts of thisData might depend.
  ## thisData should be ordered and complete.
  ##
  ## The bootFstep only ever does thisDeterm and thisData in sequence, so they can be combined. Same for auxFstep
  ##
  ## In bootFstep, the prevDeterm and thisNode could be combined, but that is not the case in auxFstep,
  ## because it breaks apart thisNode for doing the lookahead
  prevDeterm <- model$getDependencies(prevNode, determOnly = TRUE, downstream = FALSE)
  ## Weed out any prevDeterm elements that do not lead to part of thisNode
  if(length(prevDeterm) > 0) {
    keep_prevDeterm <- logical(length(prevDeterm))
    for(i in seq_along(prevDeterm)) {
      theseDeps <- model$getDependencies(prevDeterm[i])
      keep_prevDeterm[i] <- any(theseDeps %in% calc_thisNode_self)
    }
    prevDeterm <- prevDeterm[ keep_prevDeterm ]
  }
  ## Weed out any prevDeterm elements that are redundant with deterministic parts of calc_thisNode_self
  ##
  if(length(prevDeterm) > 0) {
    keep_prevDeterm <- logical(length(prevDeterm))
    for(i in seq_along(prevDeterm)) {
      theseDeps <- model$getDependencies(prevDeterm[i], determOnly = TRUE)
      keep_prevDeterm[i] <- !any(theseDeps %in% calc_thisNode_self)
    }
    prevDeterm <- prevDeterm[ keep_prevDeterm ]
  }

  ## So far we have prevNode -> prevDeterm -> calc_thisNode_self -> calc_thisNode_deps
  ## calc_thisNode_deps combines the old thisDeterm and thisData.
  ## Separating out determ and data calcs is not necessary, but we do need to skip calcs
  ## that lead to the next latent nodes.
  ## A question is whether we condition only on data -- I think yes.

  ## Do similar weeding of calc_thisNode_deps in relation to thisData.
  ## calc_thisNode_deps includes determ+data dependencies of thisNode but excludes intermediate deterministic
  ##    nodes that need to be calculated between multiple nodes of thisNode.
  ## We will separately get determ and data dependencies,
  ##    filter determ by which ones lead to data (vs. other states)
  ## The list of determ will be only deterministic nodes that are in calc_thisNode_deps
  ## What if there is a data node in calc_thisNode_self?
  ## thisDeterm <- model$getDependencies(thisNode, determOnly = TRUE)
  ## from above, thisDeterm are determ deps on thisNode that are not needed
  ##    internally, i.e. between elements of thisNode

  ## Include self in this data to sort out the binary data assignment for the
  ## occupancy models
  thisData   <- model$getDependencies(thisNode, dataOnly = TRUE, self = FALSE)
  if(length(thisDeterm) > 0) {
    keep_thisDeterm <- logical(length(thisDeterm))
    for(i in seq_along(thisDeterm)) {
      theseDeps <- model$getDependencies(thisDeterm[i])
      keep_thisDeterm[i] <- any(theseDeps %in% thisData)
    }
    thisDeterm <- thisDeterm[keep_thisDeterm]
    calc_thisNode_deps <- model$expandNodeNames(c(thisDeterm, thisData), sort = TRUE) ## only for the sort
  } else {
    calc_thisNode_deps <- thisData
  }
  list(prevDeterm = prevDeterm,
       calc_thisNode_self = calc_thisNode_self,
       calc_thisNode_deps = calc_thisNode_deps)
}

fillIndices <- function(node, info, returnExpr = FALSE) {
  ## Fill missing indexes with full extent of that dimension.
  node <- parse(text = node)[[1]]
  if(info$nDim != length(node) - 2 || node[[1]] != '[')
    stop("findLatentNodes: invalid node expression: ", node, ".")
  for(i in seq_len(info$nDim)) {
    if(node[[i+2]] == "") {
      node[[i+2]] <- substitute(A:B,
                                list(A = info$mins[i], B = info$maxs[i]))
      if(info$mins[i] == info$maxs[i])  # avoid things like 3:3
        node[[i+2]] <- info$mins[i]
    }
  }
  if(!returnExpr)
    node <- deparse(node)
  return(node)
}

findLatentNodes <- function(model, nodes, timeIndex = NULL) {
  ## Determine set of latent 'nodes', one per time point.
  ## Note that each time point might have one node or a set of nodes.
  varName <- sapply(nodes, function(x) {
    return(model$getVarNames(nodes = x))
  })
  if(length(unique(varName)) > 1){
    stop("findLatentNodes: all latent nodes must come from same variable.")
  }
  varName <- varName[1]
  info <- model$getVarInfo(varName)

  if(length(nodes) > 1) {
    ## Check for and fill in empty dimensions if more than 1 dimension.
    ## Otherwise, assume user-provided indexing is valid.
    nodes <- sapply(nodes, fillIndices, info, returnExpr = FALSE,
                    USE.NAMES = FALSE)
  } else {
    if(nodes == varName)
      ## 'nodes' is a variable, so setup indexing
      nodes <- paste0(varName, '[', paste0(rep(',', info$nDim - 1), collapse = ''), ']')

    nodeExpr <- fillIndices(nodes, info, returnExpr = TRUE)
    indexLengths <- sapply(nodeExpr[3:length(nodeExpr)],
                           function(x) length(eval(x)))
    ## Determine time index as longest dimension if not provided
    if(is.null(timeIndex)){
      maxLength <- max(indexLengths)
      if(sum(indexLengths == maxLength) > 1)
        stop("findLatentNodes: unable to determine which dimension indexes time. Specify manually using the 'timeIndex' control list argument.")
      timeIndex <- which.max(indexLengths)
      timeLength <- maxLength
    } else{
      timeLength <- indexLengths[timeIndex]
    }

    timeIndices <- nodeExpr[[timeIndex+2]]  # timeIndices is unevaluated

    ## Expand so will have one 'node' per time point, overwriting time dimension.
    nodes <- rep('', timeLength)
    cnt <- 1
    for(i in eval(timeIndices)) {
      nodeExpr[[2+timeIndex]] <- as.numeric(i)
      nodes[cnt] <- deparse(nodeExpr)
      cnt <- cnt + 1
    }
  }
  return(nodes)
}

mySetAndCalculateUpdate <- nimbleFunction(
  name = 'mySetAndCalculateUpdate',
  setup = function(model, target, latents, mvSamplesEst, my_particleFilter, m, topParamsInter, mvSaved) {# postSamples) {
    targetNodesAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    calcNodes <- model$getDependencies(target)
    latentDep <- model$getDependencies(latents)
    particleMVold <- my_particleFilter$mvEWSamples
    #extras <- any(model$expandNodeNames(extraTargetVars) %in% model$expandNodeNames(latents))
    # print(topParamsInter)
    # print(targetNodesAsScalar)
    # print(target)
  },
  run = function(iterRan = integer()) {
    index <- ceiling(runif(1, 0, m))
    #nimCopy(from = mvSaved, to = model, nodes = topParamsInter, nodesTo = topParamsInter, row = 1)
    #nimCopy(from = mvSaved, to = model, nodes = extraTargetVars, row = 1)
    #print(values(model, topParamsInter))
    #for now, I am simulating the extraPars values, the idea is to use the previous values
    #model$calculate()
    #if(extras) nimCopy(from = mvSaved, to = model, nodes = extraTargetVars, nodesTo = extraTargetVars, row = 1)#model$simulate(extraTargetVars)
    #if(extras) nimCopy(from = mvSamplesEst, to = model, nodes = latents, row = iterRan, rowTo = 1)
    #model$calculate()
    #print(values(model, targetNodesAsScalar))

    #nimCopy(from = mvSamplesEst, to = model, nodes = calNodesStoch,row = iterRan)
    #model$simulate(calNodesStoch)
    # model$calculate()
    lp <-  my_particleFilter$run(m = m, iterRun = iterRan, storeModelValues = values(model, targetNodesAsScalar))
    #print(lp)
    nimCopy(from =  particleMVold, to = model, latents, latents, row = index, rowTo = 1)
    #print(values(model, targetNodesAsScalar))
    # #calculate(model, latentDep)
    # model$calculate()
    # copy(from = model, to = mvSaved, nodes = latentDep, row = 1, logProb = TRUE)
    # copy(from)
    #simulate(model)
    #calculate(model, latentDep)
    #nimCopy(from = model, to = mvSaved, nodes = latentDep, row = 1)
    #calculate(model, c(latentDep, targetNodes))
    #copy(model, mvSaved)
    #nimCopy(particleMVold, mvSaved, latents, latents, row = index)
    #values(model, targetNodesAsScalar) <<- postSamples[iterRan, targetNodesAsScalar]
    #lp <- model$calculate(calcNodes)
    returnType(double())
    return(lp)
  }
)

mySetAndCalculateUpdate11 <- nimbleFunction(
  name = 'mySetAndCalculateUpdate11',
  setup = function(model, target, latents, mvSamplesEst, my_particleFilter, m, topParamsInter, mvSaved, extraTargetVars, extras) {# postSamples) {
    targetNodesAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    calcNodes <- model$getDependencies(target)
    latentDep <- model$getDependencies(latents)
    particleMVold <- my_particleFilter$mvEWSamples
    #extras <- any(model$expandNodeNames(extraTargetVars) %in% model$expandNodeNames(latents))
  # print(topParamsInter)
  # print(targetNodesAsScalar)
  # print(target)
    },
  run = function(iterRan = integer()) {
    index <- ceiling(runif(1, 0, m))
    nimCopy(from = mvSamplesEst, to = model, nodes = topParamsInter, nodesTo = topParamsInter, row = iterRan, rowTo = 1)
    #nimCopy(from = mvSaved, to = model, nodes = extraTargetVars, row = 1)
    #print(values(model, topParamsInter))
    #for now, I am simulating the extraPars values, the idea is to use the previous values
    #model$calculate()
    if(extras) nimCopy(from = mvSaved, to = model, nodes = extraTargetVars, nodesTo = extraTargetVars, row = 1)#model$simulate(extraTargetVars)
    #if(extras) nimCopy(from = mvSamplesEst, to = model, nodes = latents, row = iterRan, rowTo = 1)
    #model$calculate()
    #print(values(model, targetNodesAsScalar))

    #nimCopy(from = mvSamplesEst, to = model, nodes = calNodesStoch,row = iterRan)
    #model$simulate(calNodesStoch)
    # model$calculate()
    lp <-  my_particleFilter$run(m = m, iterRun = iterRan, storeModelValues = values(model, targetNodesAsScalar))
    #print(lp)
     nimCopy(from =  particleMVold, to = model, latents, latents, row = index, rowTo = 1)
    #print(values(model, targetNodesAsScalar))
     # #calculate(model, latentDep)
    # model$calculate()
    # copy(from = model, to = mvSaved, nodes = latentDep, row = 1, logProb = TRUE)
    # copy(from)
    #simulate(model)
    #calculate(model, latentDep)
    #nimCopy(from = model, to = mvSaved, nodes = latentDep, row = 1)
    #calculate(model, c(latentDep, targetNodes))
    #copy(model, mvSaved)
    #nimCopy(particleMVold, mvSaved, latents, latents, row = index)
    #values(model, targetNodesAsScalar) <<- postSamples[iterRan, targetNodesAsScalar]
    #lp <- model$calculate(calcNodes)
    returnType(double())
    return(lp)
  }
)
#
# mygenerateProposalVector = function(model, target) {
#   propValueVector <- values(model,target) ## last argument specifies prec_param = FALSE
#   #returnType(double(1))
#   return(propValueVector)
# }
#
#
# mygenerateProposalVector <- nimbleFunction(
#   name = 'mygenerateProposalVector',
#   setup = function(model, target) {
#     targetNodesAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
#     #calcNodes <- model$getDependencies(targetNodes)
#   },
#   run = function(targetValues = double(1)) {
#     values(model, targetNodesAsScalar) <<- targetValues
#     #lp <- model$calculate(calcNodes)
#     returnType(double())
#     return(targetValues)
#   }
# )

#' Creates a nimbleFunction for executing the Metropolis-Hastings jumping decision,
#' and updating values in the model, or in a carbon copy modelValues object, accordingly.
#'
#' This nimbleFunction generator must be specialized to three required arguments: a model, a modelValues, and a character vector of node names.
#'
#' @param model An uncompiled or compiled NIMBLE model object.
#' @param mvSaved A modelValues object containing identical variables and logProb variables as the model. Can be created by \code{modelValues(model)}.
#' @param target A character vector providing the target node.
#' @param calcNodes A character vector representing a set of nodes in the model (and hence also the modelValues) object.
#' @author Daniel Turek
#' @export
#' @details
#' Calling decideAndJump(model, mvSaved, calcNodes) will generate a specialized nimbleFunction with four required numeric arguments:
#'
#' modelLP1: The model log-probability associated with the newly proposed value(s)
#'
#' modelLP0: The model log-probability associated with the original value(s)
#'
#' propLP1: The log-probability associated with the proposal forward-transition
#'
#' propLP0: The log-probability associated with the proposal reverse-tranisiton
#'
#' Executing this function has the following effects:
#' -- Calculate the (log) Metropolis-Hastings ratio, as logMHR = modelLP1 - modelLP0 - propLP1 + propLP0
#' -- Make the proposal acceptance decision based upon the (log) Metropolis-Hastings ratio
#' -- If the proposal is accepted, the values and associated logProbs of all calcNodes are copied from the model object into the mvSaved object
#' -- If the proposal is rejected, the values and associated logProbs of all calcNodes are copied from the mvSaved object into the model object
#' -- Return a logical value, indicating whether the proposal was accepted
myDecideAndJump <- nimbleFunction(
  name = 'myDecideAndJump',
  setup = function(model, mvSaved, target, latentAsScalar, mvSamplesEst) {
    ccList <- myMcmc_determineCalcAndCopyNodes(model, target)
    copyNodesDeterm <- ccList$copyNodesDeterm; copyNodesStoch <- ccList$copyNodesStoch  # not used: calcNodes, calcNodesNoSelf
    #calNodesStoch <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
    #calNodesStoch <- calNodesStoch[!calNodesStoch %in% model$expandNodeNames(c(target, latents))]
    #extraStoch <- copyNodesStoch[!copyNodesStoch %in% model$expandNodeNames(latents) ]
    #index <- ceiling(runif(1, 0, m))
  },
  run = function(modelLP1 = double(), modelLP0 = double(), propLP1 = double(), propLP0 = double()){#, iterRan = integer() ){#, pfModelValues = double(), predVals = double(), topParamsVals = double()) {
    logMHR <- modelLP1 - modelLP0 - propLP1 + propLP0
    jump <- decide(logMHR)
    #if(jump) {

    #} else {
    #values(model, latentAsScalar) <<- pfModelValues
    # values(model, predictivePars) <<- predVals
    #nimCopy(from = mvSamplesEst, to = model, row = iterRan, nodes = topParams)
    #  calculate(model)
    # nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
    #nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
    #nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
    #}
    returnType(logical())
    return(jump)
  }
)


## create the lists of calcNodes and copyNodes for use in MCMC samplers

myMcmc_determineCalcAndCopyNodes <- function(model, target) {
  targetExpanded <- model$expandNodeNames(target)
  modelPredictiveNodes <- model$getNodeNames(predictiveOnly = TRUE)
  targetExpandedPPbool <- targetExpanded %in% modelPredictiveNodes
  targetAllPP <- all(targetExpandedPPbool)
  targetAnyPP <- any(targetExpandedPPbool)
  ## if a particular sampler is assigned *jointly to PP and non-PP* nodes, then we're going to bail
  ## out and quit, if the option MCMCusePredictiveDependenciesInCalculations == FALSE.
  ## this is an extreme corner-case, which I think will lead to problems.
  if(targetAnyPP && !targetAllPP && !getNimbleOption('MCMCusePredictiveDependenciesInCalculations'))
    stop('cannot assign samplers jointly to posterior predictive (PP) nodes and non-PP nodes, when MCMCusePredictiveDependenciesInCalculations option is FALSE', call. = FALSE)
  ## if the sampler calling this, itself, is operating exclusively on posterior predictive nodes,
  ## then regardless of how the rest of the model is being sampled (w.r.t. inclusion of posterior predictive nodes),
  ## we'll include 'self' and all stochastic dependencies (the full markov blanket) in the calculations,
  ## which necessarily are taking place entirely within a posterior predictive network of nodes.
  ## this should lead to correct behaviour (consistent samples and joint posteriors) in all cases.
  if(targetAllPP) {
    ## when sampler is operating only on posterior predictive nodes,
    ## then always include all predictive dependencies:
    calcNodes <- model$getDependencies(target, includePredictive = TRUE)
    calcNodesNoSelf <- model$getDependencies(target, self = FALSE, includePredictive = TRUE)
    ##calcNodesPPomitted <- character()
  } else {
    ## usual case:
    calcNodes <- model$getDependencies(target)
    calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
    ##calcNodesPPomitted <- setdiff(model$getDependencies(target, includePredictive = TRUE), calcNodes)
  }
  ## copyNodes:
  copyNodes <- model$getDependencies(target, self = FALSE)
  isStochCopyNodes <- model$isStoch(copyNodes)
  copyNodesDeterm <- copyNodes[!isStochCopyNodes]
  copyNodesStoch <- copyNodes[isStochCopyNodes]
  ##
  ccList <- list(
    calcNodes = calcNodes,
    calcNodesNoSelf = calcNodesNoSelf,
    ##calcNodesPPomitted = calcNodesPPomitted,
    copyNodesDeterm = copyNodesDeterm,
    copyNodesStoch = copyNodesStoch
  )
  return(ccList)
}

sampleTopPars <- nimbleFunction(
  name = 'sampleTopPars',
  #contains = sampler_BASE,
  setup = function(model, mvSaved, topParams, mvSamplesEst, scale, latents, target, propCov1) {
    ccList <- myMcmc_determineCalcAndCopyNodes(model, topParams)
    calcNodes <- ccList$calcNodes
    copyNodesDeterm <- ccList$copyNodesDeterm; copyNodesStoch <- ccList$copyNodesStoch  # not used: calcNodes, calcNodesNoSelf
    finalTargetIndex <- max(match(model$expandNodeNames(topParams), calcNodes))
    if(!is.integer(finalTargetIndex) | length(finalTargetIndex) != 1 | is.na(finalTargetIndex[1]))   stop('problem with target node in RW_block sampler')
    calcNodesProposalStage <- calcNodes[1:finalTargetIndex]
    calcNodesDepStage <- calcNodes[-(1:finalTargetIndex)]


    propCovOriginal1 <- propCov1
    chol_propCov1 <- chol(propCov1)
    chol_propCov_scale1 <- scale * chol_propCov1
    my_setAndCalculate <- setAndCalculateDiff(model, topParams)

    #calNodesStoch <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
    #calNodesStoch <- calNodesStoch[!calNodesStoch %in% model$expandNodeNames(c(target, latents))]
    #extraStoch <- copyNodesStoch[!copyNodesStoch %in% model$expandNodeNames(latents) ]
    #index <- ceiling(runif(1, 0, m))
  },
  run = function(iterRan = integer()){#, pfModelValues = double(), predVals = double(), topParamsVals = double()) {
    nimCopy(from = mvSamplesEst, to = model, row = iterRan, rowTo = 1, nodes = target)
    nimCopy(from = mvSamplesEst, to = model, row = iterRan, rowTo = 1, nodes = latents, nodesTo = latents)
    #model$simulate(extraTargetVars)

    # Update model to take into accout the MCMC output before generating new samples
    model$calculate()
    propVector <- generateProposalVector()
    lpD <- my_setAndCalculate$run(propVector)
    #values(model, topParams) <<- exp(propVector)
    #lpD <- model$calculateDiff(calcNodesProposalStage)

    if(lpD == -Inf) {
      jump <- FALSE
      nimCopy(from = mvSamplesEst, to = model, row = iterRan, rowTo = 1, nodes = topParams)
      model$calculate()
      nimCopy(from = model, to = mvSaved,   row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
    } else {
      lpD <- lpD + model$calculateDiff(calcNodesDepStage)
      jump <- decide(lpD)
      if(jump) {
        ##model$calculate(calcNodesPPomitted)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
      } else {
        nimCopy(from = mvSamplesEst, to = model, nodes = target, nodesTo = target, row = iterRan, rowTo = 1)
        #model$calculate()
        #nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
        #nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
        #nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
      }
    }

    #logMHR <- modelLP1 - modelLP0 - propLP1 + propLP0
    #jump <- decide(logMHR)
    #if(jump) {

    #} else {
    #values(model, latentAsScalar) <<- pfModelValues
    # values(model, predictivePars) <<- predVals
    #nimCopy(from = mvSamplesEst, to = model, row = iterRan, nodes = topParams)
    #  calculate(model)
    # nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
    #nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
    #nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
    #}
    returnType(double())
    return(lpD)
  },
  methods = list(
    generateProposalVector = function() {
      propValueVector <- rmnorm_chol(1, values(model,topParams), chol_propCov_scale1, 0)  ## last argument specifies prec_param = FALSE
      returnType(double(1))
      return(propValueVector)
    }
  )
)


# mySetAndCalculateUpdate1 <- nimbleFunction(
#   name = 'mySetAndCalculateUpdate1',
#   setup = function(model, targetNodes, mvSamplesEst, myParticleFilter,m) {# postSamples) {
#     targetNodesAsScalar <- model$expandNodeNames(targetNodes, returnScalarComponents = TRUE)
#     calcNodes <- model$getDependencies(targetNodes)
#   },
#   run = function(iterRan = integer(0)) {
#     nimCopy(from = mvSamplesEst, to = model, nodes = targetNodes,row = iterRan)
#     my_particleFilter$run(m = m, iterRun = iterRan, storeModelValues = values(model, targetNodes))
#     #values(model, targetNodesAsScalar) <<- postSamples[iterRan, targetNodesAsScalar]
#     lp <- my_particleFilter$mvEWSamples
#     returnType(double())
#     return(lp)
#   }
# )

#' @rdname samplers
#' @export
mySampler_RW_block <- nimbleFunction(
  name = 'mySsampler_RW_block',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## control list extraction
    adaptive            <- extractControlElement(control, 'adaptive',            TRUE)
    adaptScaleOnly      <- extractControlElement(control, 'adaptScaleOnly',      FALSE)
    adaptInterval       <- extractControlElement(control, 'adaptInterval',       200)
    adaptFactorExponent <- extractControlElement(control, 'adaptFactorExponent', 0.8)
    scale               <- extractControlElement(control, 'scale',               1)
    propCov             <- extractControlElement(control, 'propCov',             'identity')
    tries               <- extractControlElement(control, 'tries',               1)
    maxDimCovHistory    <- extractControlElement(control, 'maxDimCovHistory',    10)
    ## node list generation
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    ccList <- myMcmc_determineCalcAndCopyNodes(model, target)
    calcNodes <- ccList$calcNodes; copyNodesDeterm <- ccList$copyNodesDeterm; copyNodesStoch <- ccList$copyNodesStoch   # not used: calcNodesNoSelf
    finalTargetIndex <- max(match(model$expandNodeNames(target), calcNodes))
    if(!is.integer(finalTargetIndex) | length(finalTargetIndex) != 1 | is.na(finalTargetIndex[1]))   stop('problem with target node in RW_block sampler')
    calcNodesProposalStage <- calcNodes[1:finalTargetIndex]
    calcNodesDepStage <- calcNodes[-(1:finalTargetIndex)]
    ## numeric value generation
    scaleOriginal <- scale
    timesRan      <- 0
    timesAccepted <- 0
    timesAdapted  <- 0
    d <- length(targetAsScalar)
    scaleHistory <- c(0, 0)                                                                  ## scaleHistory
    acceptanceHistory <- c(0, 0)                                                             ## scaleHistory
    propCovHistory <- if(d <= maxDimCovHistory) array(0, c(2,d,d)) else array(0, c(2,2,2))   ## scaleHistory
    saveMCMChistory <- if(getNimbleOption('MCMCsaveHistory')) TRUE else FALSE
    if(is.character(propCov) && propCov == 'identity')     propCov <- diag(d)
    propCovOriginal <- propCov
    chol_propCov <- chol(propCov)
    chol_propCov_scale <- scale * chol_propCov
    empirSamp <- matrix(0, nrow=adaptInterval, ncol=d)
    ## nested function and function list definitions
    targetNodesAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    my_calcAdaptationFactor <- calcAdaptationFactor(d, adaptFactorExponent)
    ## checks
    if(any(model$isDiscrete(target)))       warning('cannot use RW_block sampler on discrete-valued target')  # This will become an error once we fix the designation of distributions in nimbleSCR to not be discrete.
    if(!inherits(propCov, 'matrix'))        stop('propCov must be a matrix\n')
    if(!inherits(propCov[1,1], 'numeric'))  stop('propCov matrix must be numeric\n')
    if(!all(dim(propCov) == d))             stop('propCov matrix must have dimension ', d, 'x', d, '\n')
    if(!isSymmetric(propCov))               stop('propCov matrix must be symmetric')
  },
  run = function() {
    for(i in 1:tries) {
      propValueVector <- generateProposalVector()
      values(model, targetNodesAsScalar) <<- propValueVector
      lpD <- model$calculateDiff(calcNodesProposalStage)
      if(lpD == -Inf) {
        jump <- FALSE
        nimCopy(from = mvSaved, to = model,   row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
      } else {
        lpD <- lpD + model$calculateDiff(calcNodesDepStage)
        jump <- decide(lpD)
        if(jump) {
          ##model$calculate(calcNodesPPomitted)
          nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
          nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
          nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
        } else {
          nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
          nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
          nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
        }
      }
      if(adaptive)     adaptiveProcedure(jump)
    }
  },
  methods = list(
    generateProposalVector = function() {
      propValueVector <- rmnorm_chol(1, values(model,target), chol_propCov_scale, 0)  ## last argument specifies prec_param = FALSE
      returnType(double(1))
      return(propValueVector)
    },
    adaptiveProcedure = function(jump = logical()) {
      timesRan <<- timesRan + 1
      if(jump)     timesAccepted <<- timesAccepted + 1
      if(!adaptScaleOnly)     empirSamp[timesRan, 1:d] <<- values(model, target)
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

