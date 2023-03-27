#' Compare models fitted with the proposed framework
#'
#' @description Returns a dataframe with the times the models run, the effective sample size of each model, the efficiency of each model (where efficiency = ESS/times run) and the Monte Carlo standard error
#'
#' @param models A list of MCMC results estimated from the sparta Updating functions in this paper
#' @param n.chains The number of chains for the models. The n.chains should be equal for all models.
#' @param modelNames The names to be assigned to each model in the models parameter.
#' The default is NULL and the model names are created with the names 'Model 1', Model 2, ... Model N (where N = n.chains)
#' @param method  The methods to be used in estimating the Monte Carlo standard error.
#' The choices of these should be the same as there is in the mcmcse package. The default chosen is batch means ("bm").
#' @param metrics list of logical responses for whether the following parameters should be returned: timeRun, ESS, efficiency, MCse.
#' @author  Kwaku Peprah Adjei
#' @export
#'
#' @family particle filtering methods
#' @references Flegal, J. M., Hughes, J., Vats, D., Dai, N., Gupta, K., Maji, U., ... & Rcpp, L. (2017). Package ‘mcmcse’.
#'
#' @examples

compareModelsPars <- function(models = list(),
                              n.chains = NULL,
                              modelNames = NULL,
                              method = "bm",
                              metrics = list(timeRun = TRUE,
                                             ESS = TRUE,
                                             efficiency = TRUE,
                                             MCse = TRUE)){
  #retrieving parameters needed for the function
  timeRun <- metrics[["timeRun"]]
  ESS <- metrics[["ESS"]]
  efficiency <- metrics[["efficiency"]]
  MCse <- metrics[["MCse"]]

  #function for g
  # interesting in estimating second moments
  #g needs to be defined
  g <- function(x){
  return(sum(x^2))
  }

  #assign names for models
  modelsLength <- length(models)
  if(is.null(modelNames)) modelNames = paste("Model", 1:modelsLength)

  #Define default number of chains
  if(is.null(n.chains)) n.chains = 1

  ######################
  # Times run
  ##########################

  if(timeRun) timesRet <- lapply(models, function(x){
    nDim <- length(x$timeRun)
    if(nDim == 1){
      x$timeRun
    }else{
      x$timeRun$all.chains
    }
  })%>%
    do.call('c', .)

  ######################
  # Effective sample size
  ##########################

  if(ESS) ESSret <- lapply(models, function(x){
    nDim <- length(x$samples[[1]])
    if(nDim > 2){
      ret <- lapply(as.list(1:n.chains), function(y){mcmcse::multiESS(as.matrix(x$samples),
                                                                      method = method,
                                                                      r = 1,
                       size = NULL, g = NULL, adjust = TRUE,
                       blather = TRUE)})%>%
        do.call('c', .)%>%
        mean(.)
    }else{
      ret <- lapply(as.list(1:n.chains), function(y){mcmcse::multiESS(as.matrix(x$samples[[y]][[2]]),
                                                                      method = method,
                                                                      r = 1,
                                                               size = NULL,
                                                               g = NULL,
                                                               adjust = TRUE,
                                                               blather = TRUE)})%>%
        do.call('c', .)%>%
        mean(.)
    }
  })%>%
    do.call('c', .)

  #setting names for efficiency
  colnames(ESSret) <- modelNames

  #####################
  # Monte Carlo Sample Error
  ######################
  if(MCse) MCseret <- lapply(models, function(x){
    nDim <- length(x$samples[[1]])
    if(nDim > 2){
      ret <- lapply(as.list(1:n.chains), function(y){
        N <- nrow(as.matrix(x$samples))
        mcse <- mcmcse::mcse.multi(as.matrix(x$samples),
                                   method = "bm",
                                   g = g,
                         blather = TRUE)
        mcseRet <- c(mcse$cov/N)
        return(mcseRet)
        }%>%
          do.call("c", .)%>%
          mean(.))
    }else{
      ret <- lapply(as.list(1:n.chains), function(y){
        N <- nrow(as.matrix(x$samples[[y]]))
        mcse <- mcmcse::mcse.multi(as.matrix(x$samples[[y]]),
                                   method = "bm",
                                   g = g,
                                   blather = TRUE)
        mcseRet <- c(mcse$cov/N)
        return(mcseRet)
      }%>%
        do.call("c", .)%>%
        mean(.))
    }
  })%>%
    do.call('rbind', .)

  #setting names for model
  rownames(MCseret) <- modelNames

#########################
# estimating efficiency
#######################

  if(efficiency) efficiencyRet <- ESSret/timesRet
  #setting names for efficiency
  colnames(efficiencyRet) <- modelNames

  retDataFrame <- data.frame(timesRun = timesRet,
                             ess = ESSret,
                             efficiency = efficiencyRet,
                             mcse = MCseret)
  return(retDataFrame)
}


#' Compare models fitted with the proposed framework for some specific parameters
#'
#' @description Returns a dataframe with the times the models run, the effective sample size of each model, the efficiency of each model (where efficiency = ESS/times run) and the Monte Carlo standard error
#'
#' @param models A list of MCMC results estimated from the sparta Updating functions in this paper
#' @param modelNames The names to be assigned to each model in the models parameter.
#' @param n.chains The number of chains for the models. The n.chains should be equal for all models.
#' The default is NULL and the model names are created with the names 'Model 1', Model 2, ... Model N (where N = n.chains)
#' @param nodes The parameters we are interested in retrieving the ESS, efficiency and Monte Carlo Standard error for.
#' @param method  The methods to be used in estimating the Monte Carlo standard error.
#' The choices of these should be the same as there is in the mcmcse package. The default chosen is batch means ("bm").
#' @param metrics list of logical responses for whether the following parameters should be returned: timeRun, ESS, efficiency, MCse.
#' The MCse is estimated with mcmcse package and ESS was estimated with ggmcmc package in R
#' @author  Kwaku Peprah Adjei
#' @export
#'
#' @family particle filtering methods
#' @references Fernández-i-Marín, X. (2013). Using the ggmcmc Package.
#' @references Flegal, J. M., Hughes, J., Vats, D., Dai, N., Gupta, K., Maji, U., ... & Rcpp, L. (2017). Package ‘mcmcse’.
#'
#' @examples

compareModelsIndividualPars <- function(models = list(),
                                        modelNames = NULL,
                                        n.chains = NULL,
                                        nodes = c(),
                                        method = "bm", #parameterisations for mcse.mat
                              metrics = list(timeRun = TRUE,
                                             ESS = TRUE,
                                             efficiency = TRUE,
                                             MCse = TRUE)){
  timeRun <- metrics[["timeRun"]]
  ESS <- metrics[["ESS"]]
  efficiency <- metrics[["efficiency"]]
  MCse = metrics[["MCse"]]

  #assign names for models
  modelsLength <- length(models)
if(is.null(modelNames)) modelNames = paste("Model", 1:modelsLength)

  #assign default number of chains
  if(is.null(n.chains)) n.chains = 1

  #assign chain names
  chainNames <- paste0("chain", 1:n.chains)

  ##################
  # Times run
  ##################
  timesRet <- lapply(models, function(x){
    nDim <- length(x$timeRun)
    if(nDim == 1){
      x$timeRun
    }else{
      x$timeRun$all.chains
    }
  })%>%
    do.call('c', .)
  ####################
  # estimate Monte Carlo Standard error
  ##################
MCseRet <- lapply(seq_along(models), function(i){
    x <- models[[i]]
    nDim <- length(x$samples[[1]])
    if(nDim >2 ){
      ret <- lapply(as.list(1:n.chains), function(y){seEst <- mcmcse::mcse.mat(as.matrix(x$samples[[y]][,nodes]),
                                                               method = "bm",
                                                               g = NULL)%>%
                                                              as.data.frame()%>%
                                                              dplyr::select(se)
                    colnames(seEst) <- modelNames[i]
                    return(seEst)
                    })

      names(ret) <- chainNames

      ret$all.chains <- do.call("cbind", ret)%>%
        rowMeans(.)
    }else{
      ret <- lapply(as.list(1:n.chains), function(y){seEst <- mcmcse::mcse.mat(as.matrix(x$samples[[y]][[2]][,nodes]),
                                                                               method = "bm",
                                                                               g = NULL)%>%
        as.data.frame()%>%
        dplyr::select(se)
      colnames(seEst) <- modelNames[i]
      return(seEst)
      })

      names(ret) <- chainNames

      ret$all.chains <- do.call("cbind", ret)%>%
        rowMeans(.)
    }

    return(ret)
  })
names(MCseRet) <- modelNames

  #############
  # Effective sample Size
  ##############
    ESSret <- lapply(models, function(x) {ggmcmc::ggs_effective(ggs(x$samples),
                        proportion = FALSE,
                        plot =  FALSE)%>%
   dplyr::filter(Parameter %in% nodes)
    })

  #############
  # Efficiency
  ##############
  efficiencyRet <- lapply(seq_along(models), function(i){
    ESSret[[i]]%>%
      dplyr::mutate(timeRan = as.numeric(timesRet[i]),
                    efficiency = Effective/ as.numeric(timesRet[i]),
                    mcse = MCseRet[[i]]$all.chains)
  })
  names(efficiencyRet) <- modelNames

  # Results to return
retlist <- list()
if(efficiency) retlist$efficiency <- efficiencyRet
if(timeRun) retlist$timeRun <- timesRet
if(MCse) retlist$mcse <- MCseRet
if(ESS) retlist$ess <- ESSret

  return(retlist)
}


#' Compare models plot for fitted SSMs with the proposed framework for some specific parameters
#'
#' @description Returns a dataframe with the times the models run, the effective sample size of each model, the efficiency of each model (where efficiency = ESS/times run) and the Monte Carlo standard error
#'
#' @param models A list of MCMC results estimated from the sparta Updating functions in this paper
#' @param modelNames The names to be assigned to each model in the models parameter.
#' @param n.chains The number of chains for the models. The n.chains should be equal for all models.
#' The default is NULL and the model names are created with the names 'Model 1', Model 2, ... Model N (where N = n.chains)
#' @param nodes The parameters we are interested in retrieving the ESS, efficiency and Monte Carlo Standard error for.
#' @param method  The methods to be used in estimating the Monte Carlo standard error.
#' The choices of these should be the same as there is in the mcmcse package. The default chosen is batch means ("bm").
#' @param metrics list of logical responses for whether the plots of the following parameters should be returned: ESS, efficiency, MCse, traceplot and  density plot.
#' The MCse is estimated with mcmcse package and ESS was estimated with ggmcmc package in R
#' @author  Kwaku Peprah Adjei
#' @export
#'
#' @family particle filtering methods
#' @references Fernández-i-Marín, X. (2013). Using the ggmcmc Package.
#' @references Flegal, J. M., Hughes, J., Vats, D., Dai, N., Gupta, K., Maji, U., ... & Rcpp, L. (2017). Package ‘mcmcse’.
#'
#' @examples

compareModelsPlots <- function(models = list(),
                                modelNames = NULL,
                                  n.chains = NULL,
                                  nodes = c(),
                                  method = "bm", #parameterisations for mcse.mat
                                  metrics = list(ESS = TRUE,
                                                       efficiency = TRUE,
                                                       MCse = TRUE,
                                                       traceplot = TRUE,
                                                        density = TRUE)){

  ESS <- metrics[["ESS"]]
  efficiency <- metrics[["efficiency"]]
  MCse = metrics[["MCse"]]
  traceplot = metrics[["traceplot"]]
  density = metrics[["density"]]

  #assign names for models
  modelsLength <- length(models)
  if(is.null(modelNames)) modelNames = paste("Model", 1:modelsLength)

  #assign default number of chains
  if(is.null(n.chains)) n.chains = 1

  #assign chain names
  chainNames <- paste0("chain", 1:n.chains)

  ##################
  # Times run
  ############
timesRet <- lapply(models, function(x){
    nDim <- length(x$timeRun)
    if(nDim == 1){
      x$timeRun
    }else{
      x$timeRun$all.chains
    }
  })%>%
    do.call('c', .)

  ####################
  # estimate Monte Carlo Standard error
  ##################

MCseRet <- lapply(seq_along(models), function(i){
    x <- models[[i]]
    nDim <- length(x$samples[[1]])
    if(nDim >2 ){
      ret <- lapply(as.list(1:n.chains), function(y){seEst <- mcmcse::mcse.mat(as.matrix(x$samples[[y]][,nodes]),
                                                                               method = "bm",
                                                                               g = NULL)%>%
        as.data.frame()%>%
        dplyr::select(se)
      colnames(seEst) <- modelNames[i]
      return(seEst)
      })

      names(ret) <- chainNames

      ret$all.chains <- do.call("cbind", ret)%>%
        rowMeans(.)
    }else{
      ret <- lapply(as.list(1:n.chains), function(y){seEst <- mcmcse::mcse.mat(as.matrix(x$samples[[y]][[2]][,nodes]),
                                                                               method = "bm",
                                                                               g = NULL)%>%
        as.data.frame()%>%
        dplyr::select(se)
      colnames(seEst) <- modelNames[i]
      return(seEst)
      })

      names(ret) <- chainNames

      ret$all.chains <- do.call("cbind", ret)%>%
        rowMeans(.)
    }

    return(ret)
  })
  names(MCseRet) <- modelNames


  #############
  # Effective sample Size
  ##############

  ESSret <- lapply(models, function(x) {ggmcmc::ggs_effective(ggmcmc::ggs(x$samples),
                                                              proportion = FALSE,
                                                              plot =  FALSE)%>%
      dplyr::filter(Parameter %in% nodes)
  })

  #############
  # Efficiency
  ##############
efficiencyRet <- lapply(seq_along(models), function(i){
    ESSret[[i]]%>%
      dplyr::mutate(timeRan = as.numeric(timesRet[i]),
                    efficiency = Effective/ as.numeric(timesRet[i]),
                    mcse = MCseRet[[i]]$all.chains)
  })
  names(efficiencyRet) <- modelNames


  ####################
  # Plot results
  ###################

  modelNamesPlots <- rep(modelNames, each = length(nodes))

  #efficiency plot
  efficiencyPlot <- efficiencyRet%>%
    do.call("rbind", .)%>%
    mutate(model = modelNamesPlots)%>%
    ggplot2::ggplot(mapping = aes(x = Parameter,
                                  y = efficiency,
                                  group = model,
                                  col = as.factor(model)))+
    geom_point()+
    geom_line()+
    theme_bw()+
    ylab("Efficiency = ESS/Run time")

  #effective sample size plot
  effectivePlot <- efficiencyRet%>%
    do.call("rbind", .)%>%
    mutate(model = modelNamesPlots)%>%
    ggplot2::ggplot(mapping = aes(x = Parameter,
                                  y = Effective,
                                  group = model,
                                  col = as.factor(model)))+
    geom_point()+
    geom_line()+
    theme_bw()+
    ylab("Effective Sample Size (ESS)")

  mcsePlot <- efficiencyRet%>%
    do.call("rbind", .)%>%
    mutate(model = modelNamesPlots)%>%
    ggplot2::ggplot(mapping = aes(x = Parameter,
                                  y = mcse,
                                  group = model,
                                  col = as.factor(model)))+
    geom_point()+
    geom_line()+
    theme_bw()+
    ylab("Monte Carlo standard error (MCSE)")

  #traceplot plot
  traceplotRet <- lapply(seq_along(models), function(i){
    x <- models[[i]]
    ggmcmc::ggs(x$samples)%>%
      dplyr::filter(Parameter %in% nodes)%>%
      ggs_traceplot()+
      facet_wrap( ~ Parameter, nrow = ceiling(length(nodes)/4), ncol = 4, scales = "free_y")+
      ggtitle(paste("Model ", i))+
      theme_bw()
  })%>%
    ggpubr::ggarrange(plotlist = .,
                      nrow = length(models),
                      common.legend = TRUE)

  #density plot
  densityplotRet <- lapply(seq_along(models), function(i){
    x <- models[[i]]
    ggmcmc::ggs(x$samples)%>%
      dplyr::filter(Parameter %in% nodes)%>%
      ggs_density()+
      ggtitle(paste("Model ", i))+
      facet_wrap( ~ Parameter, nrow = ceiling(length(nodes)/4), ncol = 4, scales = "free_y")+
      theme_bw()
  })%>%
    ggpubr::ggarrange(plotlist = .,
                      nrow = length(models),
                      common.legend = TRUE)

  # Results to return
  retlist <- list()
 if(efficiency) retlist$efficiencyPlot <- efficiencyPlot
 if(ESS) retlist$essPlot <- effectivePlot
 if(MCse) retlist$mcsePlot <- mcsePlot
  if(traceplot) retlist$tracePlot <- traceplotRet
  if(density) retlist$densityPlot <-   densityplotRet

  return(retlist)

}


#' Compare models plot for fitted SSMs with the proposed framework for some specific parameters
#'
#' @description Returns a dataframe with the times the models run, the effective sample size of each model, the efficiency of each model (where efficiency = ESS/times run) and the Monte Carlo standard error
#'
#' @param models A list of MCMC results estimated from the sparta Updating functions in this paper
#' @param modelNames The names to be assigned to each model in the models parameter.
#' @param fullModel The full NIMBLE model used for running the updated model.
#' @param nodes The parameters we are interested in retrieving the ESS, efficiency and Monte Carlo Standard error for.
#' @param method  The methods to be used in estimating the Monte Carlo standard error.
#' The choices of these should be the same as there is in the mcmcse package. The default chosen is batch means ("bm").
#' @param metrics list of logical responses for whether the plots of the following parameters should be returned: mean, se, confint.
#' If mean = TRUE, it plots the posterior mean of the nodes (parameters)
#' If se = TRUE, it plots the posterior mean plus or minus the standard error of the nodes (parameters).
#' If confint = TRUE, it plots the posterior mean of the nodes with the 95% credible intervals.
#' The MCse is estimated with mcmcse package and ESS was estimated with ggmcmc package in R
#' @author  Kwaku Peprah Adjei
#' @export
#'
#' @family particle filtering methods
#' @references Fernández-i-Marín, X. (2013). Using the ggmcmc Package.
#' @references Flegal, J. M., Hughes, J., Vats, D., Dai, N., Gupta, K., Maji, U., ... & Rcpp, L. (2017). Package ‘mcmcse’.
#'
#' @examples

compareModelsMeanPlots <- function(models = list(),
                               modelNames = NULL,
                               fullModel = stateSpaceModel,
                               nodes = c(), #parameterisations for mcse.mat
                               metrics = list(mean= TRUE,
                                              se = TRUE,
                                              confint = TRUE)){

  mean <- metrics[["mean"]]
  se <- metrics[["se"]]
  confint <- metrics[["confint"]]

  #Retrieve the MCMC output from all the models
 allData <-  lapply(models, function(x){
    x$summary$all.chains
  })

# Mean plot
 meanPlot <- lapply(seq_along(allData), function(i){
   x <- allData[[i]]
  lengthNodes <- length(sapply(nodes, function(r){rownames(x)[grepl(r, rownames(x))]}))
  nodeNames <- as.vector(sapply(nodes, function(r){c(rownames(x)[grepl(r, rownames(x))])}))

  output <-  sapply(nodes, function(r){x[grepl(r, rownames(x)),"Mean"]})%>%
  #outputDF <- data.frame(Parameters = names(grepl(nodes, rownames(x))),
                         #output = output)
   # t()%>%
    data.frame()%>%
     dplyr::mutate(Parameters = nodeNames,
                   model = rep(paste0("model", i), lengthNodes))%>%
     dplyr::full_join(data.frame(Parameters = fullModel$expandNodeNames(nodes),
                                 model = rep(paste0("model", i), length(fullModel$expandNodeNames(nodes)))),
                      by = c("Parameters", "model"))
  colnames(output)[1] <- "mean"
  return(output)
 })%>%
   do.call("rbind", .)%>%
   ggplot(., mapping = aes(x = as.factor(Parameters), y = mean, col = model, group = model))+
   geom_point()+
   geom_line()+
   theme_bw()+
   xlab("Parameters")

 # Mean plus/minus standard deviation
 sePlot <- lapply(seq_along(allData), function(i){
   x <- allData[[i]]
   lengthNodes <- length(sapply(nodes, function(r){rownames(x)[grepl(r, rownames(x))]}))
   nodeNames <- as.vector(sapply(nodes, function(r){c(rownames(x)[grepl(r, rownames(x))])}))

   output <-  sapply(nodes, function(r){x[grepl(r, rownames(x)),c(1,3)]})%>%
     #outputDF <- data.frame(Parameters = names(grepl(nodes, rownames(x))),
     #output = output)
     t()%>%
     data.frame()%>%
     dplyr::mutate(Parameters = nodeNames,
                   model = rep(paste0("model", i), lengthNodes))%>%
     dplyr::full_join(data.frame(Parameters = fullModel$expandNodeNames(nodes),
                                 model = rep(paste0("model", i), length(fullModel$expandNodeNames(nodes)))),
                      by = c("Parameters", "model"))
   colnames(output)[1:2] <- c("mean", "se")
   return(output)
 })%>%
   do.call("rbind", .)%>%
   ggplot(., mapping = aes(x = as.factor(Parameters), y = mean, col = model, group = model))+
   geom_point()+
   geom_line()+
   geom_ribbon(aes(ymin = mean - se, ymax = mean + se, fill = model), alpha = 0.1) +
   theme_bw()+
   xlab("Parameters")

# Mean with credible interval
 confintPlot <- lapply(seq_along(allData), function(i){
   x <- allData[[i]]
   lengthNodes <- length(sapply(nodes, function(r){rownames(x)[grepl(r, rownames(x))]}))
   nodeNames <- as.vector(sapply(nodes, function(r){c(rownames(x)[grepl(r, rownames(x))])}))

   output <-  sapply(nodes, function(r){x[grepl(r, rownames(x)),c(1,4,5)]})%>%
     #outputDF <- data.frame(Parameters = names(grepl(nodes, rownames(x))),
     #output = output)
     t()%>%
     data.frame()%>%
     dplyr::mutate(Parameters = nodeNames,
                   model = rep(paste0("model", i), lengthNodes))%>%
     dplyr::full_join(data.frame(Parameters = fullModel$expandNodeNames(nodes),
                                 model = rep(paste0("model", i), length(fullModel$expandNodeNames(nodes)))),
                      by = c("Parameters", "model"))
   colnames(output)[1:3] <- c("mean", "lower", "upper")
   return(output)
 })%>%
   do.call("rbind", .)%>%
   ggplot(., mapping = aes(x = as.factor(Parameters), y = mean, col = model, group = model))+
   geom_point()+
   geom_line()+
   geom_ribbon(aes(ymin = lower, ymax = upper, fill = model), alpha = 0.1) +
   theme_bw()+
   xlab("Parameters")



  # Results to return
  retlist <- list()

  if(mean) retlist$meanPlot <-  meanPlot
  if(se) retlist$sePlot <-  sePlot
  if(confint) retlist$confintPlot <- confintPlot

  return(retlist)

}
