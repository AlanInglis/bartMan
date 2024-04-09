#' permVint
#'
#' @description A variable interaction evaluation which creates a null model by
#' permuting the response, rebuilding the model, and calculating the inclusion proportion (IP)
#' of adjacent splits on the null model.
#' The final result displayed is the original model's IP minus the null IP.
#'
#' @param model Model created from either the BART, dbarts or bartMachine packages.
#' @param data A data frame containing variables in the model.
#' @param trees A data frame created by extractTreeData function.
#' @param response The name of the response for the fit.
#' @param numTreesPerm The number of trees to be used in the null model.
#' As suggested by Chipman (2009), a small number of trees is recommended (~20) to force important
#' variables to used in the model. If NULL, then the number of trees from the true model is used.
#' @param top Display only the top X interactions.
#' @return A variable interaction plot. Note that for a dbarts fit, due to the internal workings of
#' dbarts, the null model is hard-coded to 20 trees, a burn-in of 100, and 1000 iterations. Both a
#' BART and bartMachine null model will extract the identical parameters from the original model.
#'
#'
#' @importFrom dplyr tibble
#' @importFrom dplyr mutate
#' @importFrom dplyr arrange
#' @import ggplot2
#' @export
#'


permVint <- function(model,
                     data,
                     trees,
                     response,
                     numTreesPerm = NULL,
                     top = NULL) {

  # get og model vints
  actualVint <- viviBart(trees = trees, out = 'vint')

  # get null permutation vints
  permVint <- permBartVint(
  model = model,
  data = data,
  response = response,
  numTreesPerm = numTreesPerm
  )

  # final vints
  finalVintDF <- actualVint
  finalVintMat <- actualVint$propMean - permVint$propMean

  finalDf <- data.frame(
    var = actualVint$var,
    propMean = finalVintMat
  )

  if(!is.null(top)){
    finalDf <- finalDf |> arrange(-propMean) |> filter(row_number() %in% 1:top)
  }


  finalDf$propMean[finalDf$propMean<0] <- 0


    p <- finalDf %>%
      arrange(propMean) %>%
      mutate(Variable = factor(var, unique(var))) %>%
      ggplot() +
      aes(x = Variable, y = propMean) +
      geom_bar(aes(x = Variable, y = propMean), stat = "identity", fill = "steelblue", col = "black") +
      theme_light() +
      coord_flip() +
      theme_bw() +
      xlab("Variable") +
      ylab("Importance") +
      theme(
        axis.title.y = element_text(angle = 90, vjust = 0.5),
        legend.key.size = unit(0.5, "cm")
      )






  return(p)
}

# -------------------------------------------------------------------------

# Main function:
permBartVint <- function(model, data, response, numTreesPerm = NULL) {
  UseMethod("permBartVint")
}




# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

# BART --------------------------------------------------------------------
permBartVint.bart <- function(model, data,  response, numTreesPerm = NULL){

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

  permuteBARTFn <- function(data) {
    yPerm <- sample(data[, responseIdx], replace = FALSE)
    x <- data[, -responseIdx]

    bmodelPerm <- dbarts::bart(x.train = x,
                               y.train = yPerm,
                               ntree = 20,
                               keeptrees = TRUE,
                               nskip = 100,
                               ndpost = 1000,
                               combinechains = F,
                               nchain = 1,
                               verbose = FALSE
    )

    permDF <- extractTreeData(bmodelPerm, data)
    permVints <- viviBart(trees = permDF, out = 'vint')
    return(permVints)
  }

  perMats <- permuteBARTFn(data)

  return(perMats)

}

permBartVint.bartMachine <- function(model, data,  response, numTreesPerm = NULL) {

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

  # null model info
  responseIdx <- which((names(data) %in% response))
  if(is.null(numTreesPerm)){
    numTreesPerm <- nTree
  }

  # null model fuunction
  permuteBARTFn <- function(data){

    yPerm <- sample(data[, responseIdx], replace = FALSE)
    x <- data[, -responseIdx]

    bmodelPerm <- bartMachine::bartMachine(X = x,
                              y = yPerm,
                              num_trees = numTreesPerm,
                              flush_indices_to_save_RAM = FALSE,
                              num_burn_in = burnIn,
                              num_iterations_after_burn_in = nMCMC,
                              verbose = FALSE)




    permDF <- extractTreeData(bmodelPerm, data)
    permVints <- viviBart(trees = permDF,out = 'vint')
    return(permVints)
  }

  perMats <- permuteBARTFn(data)

  return(perMats)

}

permBartVint.wbart <- function(model, data,  response, numTreesPerm = NULL) {

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
    bmodelPerm <- BART::wbart(
      x.train = x,
      y.train = yPerm,
      nskip = burnIn,
      ndpost = nMCMC,
      nkeeptreedraws = nMCMC,
      ntree = numTreesPerm
    ),
    file = nullfile()
    )

    permDF <- extractTreeData(bmodelPerm, data)
    permVints <- viviBart(trees = permDF,  out = 'vint')
    return(permVints)
  }

  perMats <- permuteBARTFn(data)

  return(perMats)
}





