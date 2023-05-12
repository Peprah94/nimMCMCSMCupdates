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
      as.numeric(x$timeRun, units = "secs")
    }else{
      as.numeric(x$timeRun$all.chains, units = "secs")
    }
  })%>%
    do.call('c', .)

  ######################
  # Effective sample size
  ##########################

  if(ESS) ESSret <- lapply(models, function(x){
    #nDim <- length(x$samples[[1]])
    #if(nDim > 2){
      ret <- mcmcse::multiESS(as.matrix(x$samples),
                               method = method,
                       r = 1,
                       size = NULL, g = NULL, adjust = TRUE,
                       blather = TRUE)
    # }else{
    #   ret <- lapply(as.list(1:n.chains), function(y){mcmcse::multiESS(as.matrix(x$samples[[y]][[2]]),
    #                                                                   method = method,
    #                                                                   r = 1,
    #                                                            size = NULL,
    #                                                            g = NULL,
    #                                                            adjust = TRUE,
    #                                                            blather = TRUE)})%>%
    #     do.call('c', .)%>%
    #     mean(.)
    # }
  })%>%
    do.call('c', .)

  #setting names for efficiency
  colnames(ESSret) <- modelNames

  #####################
  # Monte Carlo Sample Error
  ######################
  if(MCse) MCseret <- lapply(models, function(x){
    #nDim <- length(x$samples[[1]])
    #if(nDim > 2){
      #ret <- lapply(as.list(1:n.chains), function(y){
        #N <- nrow(as.matrix(x$samples))
        mcse <- mcmcse::mcse.mat(as.matrix(x$samples),
                                   method = "bm",
                                   g = NULL)
        #mcseRet <- c(mcse$cov/N)
        #return(mcse)
      #})
        # })%>%
        #   do.call("rbind", .)#%>%
        #   #mean(.))
    # }else{
    #   ret <- lapply(as.list(1:n.chains), function(y){
    #     N <- nrow(as.matrix(x$samples[[y]]))
    #     mcse <- mcmcse::mcse.multi(as.matrix(x$samples[[y]]),
    #                                method = "bm",
    #                                g = g,
    #                                blather = TRUE)
    #     mcseRet <- c(mcse$cov/N)
    #     return(mcseRet)
    #   }%>%
    #     do.call("c", .)%>%
    #     mean(.))
    # }
  #})%>%
  })
    #do.call('rbind', .)

  #setting names for model
  #rownames(MCseret) <- modelNames

#########################
# estimating efficiency
#######################

  if(efficiency) efficiencyRet <- ESSret/timesRet
  #setting names for efficiency
  colnames(efficiencyRet) <- modelNames



  retDataFrame <- lapply(1:length(models), function(x){
    data.frame(timesRun = timesRet[x],
               ess = ESSret[x],
               efficiency = efficiencyRet[x],
               mcse = MCseret[[x]])
  })%>%
    do.call("rbind", .)

  # retDataFrame <- data.frame(timesRun = timesRet,
  #                            ess = ESSret,
  #                            efficiency = efficiencyRet,
  #                            mcse = MCseret)
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
if(timeRun)  timesRet <-  lapply(models, function(x){
      nDim <- length(x$timeRun)
      if(nDim == 1){
        as.numeric(x$timeRun, units = "secs")
      }else{
        as.numeric(x$timeRun$all.chains, units = "secs")
      }
    })%>%
      do.call('c', .)
  ####################
  # estimate Monte Carlo Standard error
  ##################
if(MCse) {MCseRet <- lapply(seq_along(models), function(i){
    x <- models[[i]]
    #nDim <- length(x$samples[[1]])
    #if(nDim >2 ){
      seEst <- mcmcse::mcse.mat(as.matrix(x$samples[,nodes]),
                                                               method = "bm",
                                                               g = NULL)%>%
                                                              as.data.frame()%>%
                                                              dplyr::select(se)
                    #colnames(seEst) <- modelNames[i]

                    # seEst <- seEst%>%
                    #   do.call("rbind", .)

                    return(seEst)
      #
      #
      # names(ret) <- chainNames
      #
      # ret$all.chains <- do.call("cbind", ret)%>%
      #   rowMeans(.)
  #   }else{
  #     ret <- lapply(as.list(1:n.chains), function(y){seEst <- mcmcse::mcse.mat(as.matrix(x$samples[[y]][[2]][,nodes]),
  #                                                                              method = "bm",
  #                                                                              g = NULL)%>%
  #       as.data.frame()%>%
  #       dplyr::select(se)
  #     colnames(seEst) <- modelNames[i]
  #     return(seEst)
  #     })
  #
  #     names(ret) <- chainNames
  #
  #     ret$all.chains <- do.call("cbind", ret)%>%
  #       rowMeans(.)
  #   }
  #
  #   return(ret)
  # })
#names(MCseRet) <- modelNames
})%>%
  do.call("cbind", .)
}
  #############
  # Effective sample Size
  ##############
 if(ESS)   ESSret <- lapply(models, function(x) {
   ret <- mcmcse::ess(as.matrix(x$samples[,nodes]))%>%
     as.data.frame()

   colnames(ret) <- "Effective"

   return(ret)
   # ggmcmc::ggs_effective(ggs(x$samples),
   #                      proportion = FALSE,
   #                      plot =  FALSE)%>%
   # dplyr::filter(Parameter %in% nodes)
    })

  #############
  # Efficiency
  ##############
 if(efficiency) {efficiencyRet <- lapply(seq_along(models), function(i){
    ESSret[[i]]%>%
      dplyr::mutate(timeRan = as.numeric(timesRet[i]),
                    efficiency = Effective/ as.numeric(timesRet[i]),
                    mcse = MCseRet[[i]])
  })
  names(efficiencyRet) <- modelNames
}
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
                               line_type = NULL,
                                  nodes = c(),
                               bacthSize = 500,
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
      as.numeric(x$timeRun, units = "secs")
    }else{
      as.numeric(x$timeRun$all.chains, units = "secs")
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
                                                                               #size = bacthSize,
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
  line_type <- rep(line_type, each = length(nodes))

  #efficiency plot
  if(is.null(line_type)){
  efficiencyPlot <- efficiencyRet%>%
    do.call("rbind", .)%>%
    mutate(model = modelNamesPlots)%>%
    ggplot2::ggplot(mapping = aes(x = Parameter,
                                  y = efficiency,
                                  group = model,
                                  col = as.factor(model)))+
    geom_point(position = position_dodge(width = 0.1),
               aes(shape = as.factor(model)))+
    #geom_line(aes(linetype = as.factor(model)))+
    theme_bw()+
    xlab("")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = "bottom")+
    #ylab("Efficiency = ESS/Run time")
    ylab("")
  }else(
    efficiencyPlot <- efficiencyRet%>%
      do.call("rbind", .)%>%
      mutate(model = modelNamesPlots)%>%
      ggplot2::ggplot(mapping = aes(x = Parameter,
                                    y = efficiency,
                                    group = model,
                                    col = as.factor(model)))+
      geom_point(position = position_dodge(width = 0.1),
                 aes(shape = as.factor(model)))+
      #geom_line(linetype = line_type)+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            legend.position = "bottom")+
      #ylab("Efficiency = ESS/Run time")
      xlab("")+
      ylab("")
  )

  #effective sample size plot
  effectivePlot <- efficiencyRet%>%
    do.call("rbind", .)%>%
    mutate(model = modelNamesPlots)%>%
    ggplot2::ggplot(mapping = aes(x = Parameter,
                                  y = Effective,
                                  group = model,
                                  col = as.factor(model)))+
    geom_point(position = position_dodge(width = 0.1),
               aes(shape = as.factor(model)))+
    #geom_line(aes(linetype = as.factor(model)))+
    theme_bw()+
    #ylab("ESS")+
    ylab("")+
    xlab("")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = "bottom")+
    labs(color = "")+
    labs(shape = "")

  mcsePlot <- efficiencyRet%>%
    do.call("rbind", .)%>%
    mutate(model = modelNamesPlots)%>%
    ggplot2::ggplot(mapping = aes(x = Parameter,
                                  y = mcse,
                                  group = model,
                                  col = as.factor(model)))+
    geom_point(position = position_dodge(width = 0.1),
               aes(shape = as.factor(model)))+
    #geom_line(aes(linetype = as.factor(model)))+
    theme_bw()+
    #ylab("Monte Carlo standard error (MCSE)")
    #ylab("MCSE")+
    ylab("")+
    xlab("")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = "bottom")+
    labs(color = "")+
    labs(shape = "")

  #traceplot plot
  traceplotRet <- lapply(seq_along(models), function(i){
    x <- models[[i]]
    ggmcmc::ggs(x$samples)%>%
      dplyr::filter(Parameter %in% nodes)%>%
      ggs_traceplot()+
      facet_wrap( ~ Parameter, nrow = ceiling(length(nodes)/4), ncol = 4, scales = "free_y")+
     # ggtitle(paste("Model ", i))+
      ggtitle(modelNames[i])+
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
      facet_wrap( ~ Parameter, nrow = ceiling(length(nodes)/4), ncol = 4, scales = "free_y")+
      ggtitle(modelNames[i])+
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
                               updated = FALSE,
                               trueData,
                               tr = NULL,
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

 if(is.null(modelNames)) modelNames <- paste("Model", 1: length(models))

 nodeNames <- fullModel$expandNodeNames(nodes)

 expandedNodes <- length(nodeNames) == length(nodes)
 truth <- data.frame(Parameters = nodeNames,
                     mean = trueData,
                     model = "truth",
                     row = seq(1:length(nodeNames)))

# Mean plot
 meanPlot <- lapply(seq_along(allData), function(i){
   x <- allData[[i]]
  #lengthNodes <- length(sapply(nodes, function(r){rownames(x)[grepl(r, rownames(x))]}))
  #nodeNames <- as.vector(sapply(nodes, function(r){c(rownames(x)[grepl(r, rownames(x))])}))
 #true values


   lengthNodes <- sum(rownames(x) %in% nodeNames)

  if(updated[i] == TRUE){
    updatedParNames <- nodeNames[tr[[i]]:length(nodeNames)]
  }else{
    updatedParNames <- rownames(x)[rownames(x) %in% nodeNames]
  }

  output <-  x[rownames(x)[rownames(x) %in% nodeNames],"Mean"]%>%
  #outputDF <- data.frame(Parameters = names(grepl(nodes, rownames(x))),
                         #output = output)
   # t()%>%
    data.frame()%>%
     dplyr::mutate(Parameters = updatedParNames,
                   model = rep(modelNames[i], lengthNodes),
                   row = row_number())%>%
     dplyr::full_join(data.frame(Parameters = fullModel$expandNodeNames(nodes),
                                 model = rep(modelNames[i], length(fullModel$expandNodeNames(nodes)))),
                      by = c("Parameters", "model"))

  if(!expandedNodes) output$Parameters <- factor(output$Parameters,
                              levels =   paste0(nodes,"[", sort(as.numeric(gsub("\\D", "", output$Parameters))), "]"))

    #dplyr::mutate(Parameters = stringr::str_sort(Parameters, numeric = TRUE))%>%
    #dplyr::mutate(Parameters = factor(Parameters, levels = row))
  colnames(output)[1] <- "mean"
  return(output)
 })%>%
   do.call("rbind", .)%>%
   rbind(., truth)%>%
   #ggplot(., mapping = aes(x = reorder(as.factor(Parameters), row), y = mean, col = model, group = model))+
   ggplot(data = ., mapping = aes(x = Parameters, y = mean, col = model, group = model))+
   geom_point(aes(shape = model))+
   geom_line()+
   scale_color_brewer(palette = "Paired")+
   #geom_point(data = truth, mapping = aes(x = Parameters, y = mean))+
  # geom_line(data = truth, mapping = aes(x = Parameters, y = mean))+
   theme_bw()+
   xlab("Parameters")+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
         legend.position = "bottom")+
   labs(color = "")+
   labs(shape = "")

 # Mean plus/minus standard deviation
 # True values
 truth <- data.frame(Parameters = nodeNames,
                     mean = trueData,
                     se = NA,
                     model = "truth",
                     row = seq(1:length(nodeNames)))

 sePlotData <- lapply(seq_along(allData), function(i){
   x <- allData[[i]]
   #lengthNodes <- length(sapply(nodes, function(r){rownames(x)[grepl(r, rownames(x))]}))
   #nodeNames <- as.vector(sapply(nodes, function(r){c(rownames(x)[grepl(r, rownames(x))])}))
   lengthNodes <- sum(rownames(x) %in% nodeNames)

   if(updated[i] == TRUE){
     updatedParNames <- nodeNames[tr[[i]]:length(nodeNames)]
   }else{
     updatedParNames <- rownames(x)[rownames(x) %in% nodeNames]
   }

   output <-  x[rownames(x)[rownames(x) %in% nodeNames],c(1,3)]%>%
     #outputDF <- data.frame(Parameters = names(grepl(nodes, rownames(x))),
     #output = output)
     # t()%>%
     data.frame()%>%
     dplyr::mutate(Parameters = updatedParNames,
                   model = rep(modelNames[i], lengthNodes),
                   row = row_number())%>%
     dplyr::full_join(data.frame(Parameters = fullModel$expandNodeNames(nodes),
                                 model = rep(modelNames[i], length(fullModel$expandNodeNames(nodes)))),
                      by = c("Parameters", "model"))

   if(!expandedNodes) output$Parameters <- factor(output$Parameters,
                                                  levels =   paste0(nodes,"[", sort(as.numeric(gsub("\\D", "", output$Parameters))), "]"))

   colnames(output)[1:2] <- c("mean", "se")
   return(output)
 })%>%
   do.call("rbind", .)%>%
   rbind(., truth)

 if(!expandedNodes){
   sePlot  <- ggplot(sePlotData , mapping = aes(x = Parameters, y = mean, col = model, group = model))+
   geom_point(aes(col = model), position = position_dodge(width = 0.4))+
   geom_line(position = position_dodge(width = 0.4))+
   geom_errorbar(aes(ymin = mean - se, ymax = mean + se, col = model),
                 position = position_dodge(width = 0.4)) +
   scale_color_brewer(palette = "Paired")+
   theme_bw()+
   xlab("")+
     ylab("")+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
         legend.position = "bottom")+
     labs(color = "")+
     labs(shape = "")#+
     #ylim(c(-20,20))
 }else{
   sePlot  <-  ggplot(sePlotData , mapping = aes(x = Parameters, y = mean, col = model, group = model))+
     geom_point(aes(col = model), position = position_dodge(width = 0.4))+
     #geom_line()+
     geom_errorbar(aes(ymin = mean - se, ymax = mean + se, col = model),
                   position = position_dodge(width = 0.4)#,
                   #width = 0.05,
                  # size = 1
                   ) +
     #scale_color_brewer(palette = "Paired")+
     theme_bw()+
     xlab("")+
     ylab("")+
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
           legend.position = "bottom")+
     labs(color = "")+
     labs(shape = "")
 }

# Mean with credible interval
 truth <- data.frame(Parameters = nodeNames,
                     mean = trueData,
                     lower = NA,
                     upper = NA,
                     model = "truth",
                     row = seq(1:length(nodeNames)))

 confintPlot <- lapply(seq_along(allData), function(i){
   x <- allData[[i]]
   #lengthNodes <- length(sapply(nodes, function(r){rownames(x)[grepl(r, rownames(x))]}))
   #nodeNames <- as.vector(sapply(nodes, function(r){c(rownames(x)[grepl(r, rownames(x))])}))

   lengthNodes <- sum(rownames(x) %in% nodeNames)

   if(updated[i] == TRUE){
     updatedParNames <- nodeNames[tr[[i]]:length(nodeNames)]
   }else{
     updatedParNames <- rownames(x)[rownames(x) %in% nodeNames]
   }

   output <-  x[rownames(x)[rownames(x) %in% nodeNames],c(1,4,5)]%>%
     #outputDF <- data.frame(Parameters = names(grepl(nodes, rownames(x))),
     #output = output)
     # t()%>%
     data.frame()%>%
     dplyr::mutate(Parameters = updatedParNames,
                   model = rep(modelNames[i], lengthNodes),
                   row = row_number())%>%
     dplyr::full_join(data.frame(Parameters = fullModel$expandNodeNames(nodes),
                                 model = rep(modelNames[i], length(fullModel$expandNodeNames(nodes)))),
                      by = c("Parameters", "model"))

   if(!expandedNodes) output$Parameters <- factor(output$Parameters,
                                                  levels =   paste0(nodes,"[", sort(as.numeric(gsub("\\D", "", output$Parameters))), "]"))

   colnames(output)[1:3] <- c("mean", "lower", "upper")
   return(output)
 })%>%
   do.call("rbind", .)%>%
   rbind(., truth)%>%
   ggplot(., mapping = aes(x = Parameters, y = mean, col = model, group = model))+
   geom_point(aes(shape = model))+
   geom_line()+
   geom_ribbon(aes(ymin = lower, ymax = upper, fill = model), alpha = 0.1) +
   theme_bw()+
   scale_color_brewer(palette = "Paired")+
   xlab("Parameters")+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
         legend.position = "bottom")+
   labs(color = "")+
   labs(shape = "")



  # Results to return
  retlist <- list()

  if(mean) retlist$meanPlot <-  meanPlot
  if(se) retlist$sePlot <-  sePlot
  if(confint) retlist$confintPlot <- confintPlot

  return(retlist)

}
