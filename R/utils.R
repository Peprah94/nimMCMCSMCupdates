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
  thisData   <- model$getDependencies(thisNode, dataOnly = TRUE)
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
  },
  run = function(iterRan = integer()) {
    index <- ceiling(runif(1, 0, m))
    nimCopy(from = mvSamplesEst, to = model, nodes = target,row = iterRan)
    #nimCopy(from = mvSaved, to = model, nodes = extraTargetVars, row = 1)
    #for now, I am simulating the extraPars values, the idea is to use the previous values
    #model$simulate(extraTargetVars)
    model$calculate()
    #nimCopy(from = mvSamplesEst, to = model, nodes = calNodesStoch,row = iterRan)
    #model$simulate(calNodesStoch)
    # model$calculate()
    lp <-  my_particleFilter$run(m = m, iterRun = iterRan, storeModelValues = values(model, target))
    nimCopy(from =  particleMVold, to = model, latents, latents, row = index)
    calculate(model, latentDep)
    copy(from = model, to = mvSaved, nodes = latentDep, row = 1, logProb = TRUE)
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
  run = function(modelLP1 = double(), modelLP0 = double(), propLP1 = double(), propLP0 = double(), iterRan = integer() ){#, pfModelValues = double(), predVals = double(), topParamsVals = double()) {
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
  setup = function(model, mvSaved, topParams, mvSamplesEst, scale, latents, target) {
    ccList <- myMcmc_determineCalcAndCopyNodes(model, topParams)
    calcNodes <- ccList$calcNodes
    copyNodesDeterm <- ccList$copyNodesDeterm; copyNodesStoch <- ccList$copyNodesStoch  # not used: calcNodes, calcNodesNoSelf
    finalTargetIndex <- max(match(model$expandNodeNames(topParams), calcNodes))
    if(!is.integer(finalTargetIndex) | length(finalTargetIndex) != 1 | is.na(finalTargetIndex[1]))   stop('problem with target node in RW_block sampler')
    calcNodesProposalStage <- calcNodes[1:finalTargetIndex]
    calcNodesDepStage <- calcNodes[-(1:finalTargetIndex)]

    d1 <- length(topParams)
    propCov1 <- diag(d1)
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
        nimCopy(from = mvSamplesEst, to = model, row = iterRan, rowTo = 1, nodes = topParams)
        model$calculate()
        nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
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
