#' permVimp
#'
#' @description A variable selection approach which creates a null model by
#' permuting the response, rebuilding the model, and calculating the inclusion proportion (IP) on the null model.
#' The final result displayed is the original model's IP minus the null IP.
#'
#' @param model Model created from either the BART, dbarts or bartMachine packages.
#' @param data A data frame containing variables in the model.
#' @param response The name of the response for the fit.
#' @param numTreesPerm The number of trees to be used in the null model.
#' As suggested by Chipman (2009), a small number of trees is recommended (~20) to force important
#' variables to used in the model. If NULL, then the number of trees from the true model is used.
#' @param plotType Either a bar plot ('barplot') or a point plot ('point')
#' @return A variable selection plot.
#'
#'
#' @importFrom dplyr tibble
#' @importFrom dplyr mutate
#' @importFrom dplyr arrange
#' @import ggplot2
#' @importFrom utils capture.output
#'
#' @examples
#'if (requireNamespace("dbarts", quietly = TRUE)) {
#'  # Load the dbarts package to access the bart function
#'  library(dbarts)
#'  # Get Data
#'  df <- na.omit(airquality)
#'  # Create Simple dbarts Model For Regression:
#'  set.seed(1701)
#'  dbartModel <- bart(df[2:6],
#'    df[, 1],
#'    ntree = 5,
#'    keeptrees = TRUE,
#'    nskip = 10,
#'    ndpost = 10
#'  )
#'  # Tree Data
#'  trees_data <- extractTreeData(model = dbartModel, data = df)
#'  permVimp(model = dbartModel, data = df, response = 'Ozone', numTreesPerm = 2, plotType = 'point')
#'}
#'
#' @export

permVimp <- function(model, data, response, numTreesPerm = NULL, plotType = 'barplot') {

  vimp <- perBart(
    model = model,
    data = data,
    response = response,
    numTreesPerm = numTreesPerm
  )

  vimpPlot <- permPlotFn(dat = vimp,
                         plotType = plotType)

  return(vimpPlot)
}



# -------------------------------------------------------------------------

# Main function:
perBart <- function(model, data, response, numTreesPerm = NULL) {
  UseMethod("perBart")
}



# BART --------------------------------------------------------------------


perBart.wbart <- function(model, data,  response, numTreesPerm = NULL) {

  if (!requireNamespace("BART", quietly = TRUE)) {
    stop("Package \"BART\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  # get model info
  modelTrees <- model$treedraws$trees
  modelInfo <- unlist(strsplit(modelTrees, " "))[1:3]
  modelInfo <- gsub("(^\\d+)([\a-zA-Z0-9]*)", "\\1", modelInfo)

  nMCMC <- as.integer(modelInfo[1])
  nTree <- as.integer(modelInfo[2])
  nVar <- as.integer(modelInfo[3])
  burnIn <- length(model$sigma) - nMCMC

  # get var inc props
  varProp <- model$varcount
  varPropAvg <- proportions(varProp, 1)

  # null model info
  responseIdx <- which((names(data) %in% response))
  if(is.null(numTreesPerm)){
    numTreesPerm <- nTree
  }

  # null model function
  permuteBARTFn <- function(data) {
    yPerm <- sample(data[, responseIdx], replace = FALSE)
    x <- data[, -responseIdx]


    # capture.output is used to suppress output of building model
    capture.output(
    bmodelPerm <- BART:: wbart(
      x.train = x,
      y.train = yPerm,
      nskip = burnIn,
      ndpost = nMCMC,
      nkeeptreedraws = nMCMC,
      ntree = numTreesPerm
    ),
    file = nullfile()
    )


    varPropsPerm <- bmodelPerm$varcount
    varPropsPerm <- proportions(varPropsPerm, 1)
    return(varPropsPerm)
  }

  perMats <- permuteBARTFn(data)
  finalMat <- varPropAvg - perMats

  return(finalMat)
}


# dbarts ------------------------------------------------------------------

perBart.bart <-  function(model, data,  response, numTreesPerm = NULL) {

  if (!requireNamespace("dbarts", quietly = TRUE)) {
    stop("Package \"dbarts\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  # get some information
  nTree <- model$call$ntree
  nMCMC  <- model$call$ndpost
  nVar  <- as.integer(length(colMeans((model$varcount))))
  varNames <- colnames(model$fit$data@x)
  burnIn <-  model$call$nskip

  # get var inc props
  varProp <- model$varcount
  varPropAvg <- proportions(varProp, 1)

  # null model info
  responseIdx <- which((names(data) %in% response))
  if(is.null(numTreesPerm)){
    numTreesPerm <- nTree
  }

  # null model function
  permuteDBART <- function(data){

    yPerm <- sample(data[, responseIdx], replace = FALSE)
    x <- data[, -responseIdx]


    bmodelPerm <- dbarts::bart(x.train = x,
                       y.train = yPerm,
                       ntree = numTreesPerm,
                       keeptrees = TRUE,
                       nskip = burnIn,
                       ndpost = nMCMC,
                       combinechains = F,
                       nchain = 1,
                       verbose = FALSE
    )


    varPropsPerm <- bmodelPerm$varcount
    varPropsPermAvg <- proportions(varPropsPerm, 1)
    return(varPropsPermAvg)
  }

  perMats <- permuteDBART(data)
  finalMat <- varPropAvg - perMats

  return(finalMat)


}


# bartMachine -------------------------------------------------------------

perBart.bartMachine <- function(model, data, response, numTreesPerm = NULL){

  if (!requireNamespace("bartMachine", quietly = TRUE)) {
    stop("Package \"bartMachine\" needed for this function to work. Please install it.",
         call. = FALSE)
  }


  # get some information
  nTree <-  model$num_trees
  nMCMC <-  model$num_iterations_after_burn_in
  nVar  <- model$p
  varNames <- colnames(model$X)
  burnIn <-  model$num_burn_in

  # get var inc props
  varProp <- bartMachine::get_var_counts_over_chain(model)
  varPropAvg <- proportions(varProp, 1)

  # null model info
  responseIdx <- which((names(data) %in% response))
  if(is.null(numTreesPerm)){
    numTreesPerm <- nTree
  }

  # null model fuunction
  permuteBMachine <- function(data){

    yPerm <- sample(data[, responseIdx], replace = FALSE)
    x <- data[, -responseIdx]

    bmodelPerm <- bartMachine::bartMachine(X = x,
                              y = yPerm,
                              num_trees = numTreesPerm,
                              flush_indices_to_save_RAM = FALSE,
                              num_burn_in = burnIn,
                              num_iterations_after_burn_in = nMCMC,
                              verbose = FALSE)


    varPropPerm <- bartMachine::get_var_counts_over_chain(bmodelPerm)
    varPropAvgPerm <- proportions(varPropPerm, 1)
    return(varPropAvgPerm)
  }

  perMats <- permuteBMachine(data)
  finalMat <- varPropAvg - perMats

  return(finalMat)
}


# Plotting function -------------------------------------------------------


permPlotFn <- function(dat, plotType = 'barplot'){

  points <- dplyr::tibble(
    variable = colnames(dat),
    sds = apply(dat, 2, sd),
    se = apply(dat, 2, function(a) sd(a) / sqrt(length(a))),
    mean = pmax(apply(dat, 2, mean), 0),
    low = pmax(mean - 2 * se, 0),
    high = pmax(mean + 2 * se, 0)
  )

  if (plotType == "barplot") {
    p <- points %>%
      arrange(mean) %>%
      mutate(Variable = factor(variable, unique(variable))) %>%
      ggplot() +
      aes(x = Variable, y = mean) +
      geom_bar(aes(x = Variable, y = mean), stat = "identity", fill = "steelblue", col = "black") +
      geom_segment(aes(x = Variable, xend = Variable, y = low, yend = high), color = "black") +
      theme_light() +
      coord_flip() +
      theme_bw() +
      xlab("Variable") +
      ylab("Importance") +
      theme(
        axis.title.y = element_text(angle = 90, vjust = 0.5),
        legend.key.size = unit(0.5, "cm")
      )
  } else if (plotType == "point") {

    if (!requireNamespace("ggforce", quietly = TRUE)) {
      stop("Package \"ggforce\" needed for this function to work. Please install it.",
           call. = FALSE)
    }

    p <- points %>%
      arrange(mean) %>%
      mutate(Variable = factor(variable, unique(variable))) %>%
      ggplot(aes(x = Variable, y = mean)) +
      ggforce::geom_link(aes(
        x = Variable, xend = Variable, yend = low,
        col = Variable, alpha = rev(after_stat(index))
      ),
      size = 5, n = 1000
      ) +
      ggforce::geom_link(aes(
        x = Variable, xend = Variable, yend = high,
        col = Variable, alpha = rev(after_stat(index))
      ),
      size = 5, n = 1000
      ) +
      geom_point(aes(x = Variable, y = mean), shape = 18, size = 2, color = "black") +
      coord_flip() +
      theme_bw() +
      labs(x = "Variable", y = "Importance") +
      theme(legend.position = "none")
  }

  return(p)
}
