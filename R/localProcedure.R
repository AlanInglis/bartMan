#' localProcedure
#'
#' @description A variable selection approach performed by permuting the response.
#'
#' @param model Model created from either the BART, dbarts or bartMachine packages.
#' @param data A data frame containing variables in the model.
#' @param numRep The number of replicates to perform for the BART null model's variable inclusion proportions.
#' @param numTreesRep The number of trees to be used in the replicates.
#' As suggested by Chipman (2009), a small number of trees is recommended (~20) to force important
#' variables to used in the model. If NULL, then the number of trees from the true model is used.
#' @param alpha The cut-off level for the thresholds.
#' @param shift Whether to shift the incusion proportion points by the difference
#' in distance between the quantile and the value of the inclusion proportion point.
#'
#' @return A variable selection plot using the local procedure method.
#'
#'
#' @importFrom BART wbart
#' @importFrom dbarts bart
#' @importFrom bartMachine bartMachine
#' @importFrom bartMachine get_var_props_over_chain
#' @importFrom dplyr tibble
#' @import ggplot2
#'
#' @export

localProcedure <- function(model, data, numRep = 10, numTreesRep = NULL, alpha = 0.5, shift = FALSE){
  vimp <- lProd(model= model,
                data = data,
                numRep = numRep,
                numTreesRep = numTreesRep,
                alpha = alpha,
                shift = shift)
  return(vimp)
}


# -------------------------------------------------------------------------

# Main function:
lProd <- function(model, data, numRep = 10, numTreesRep = NULL, alpha = 0.5, shift = FALSE) {
  UseMethod("lProd")
}


# BART --------------------------------------------------------------------

lProd.wbart <- function(model, data, numRep = 10, numTreesRep = NULL, alpha = 0.5, shift = FALSE){

  # get some information
  modelTrees <- model$treedraws$trees
  modelInfo <- unlist(strsplit(modelTrees, " "))[1:3]
  modelInfo <- gsub("(^\\d+)([\a-zA-Z0-9]*)", "\\1", modelInfo)

  nMCMC <- as.integer(modelInfo[1])
  nTree <- as.integer(modelInfo[2])
  nVar  <- as.integer(modelInfo[3])
  burnIn <- length(model$sigma) - nMCMC

  # set up matrix
  permuteMat <- matrix(NA, nrow = numRep, ncol = nVar)
  colnames(permuteMat) <- colnames(model$varprob)
  varProp <- model$varcount
  varPropAvg <- colMeans(proportions(varProp, 1))
  varPropAvg <- sort(varPropAvg, decreasing = TRUE)


  responseIdx <- which(!(names(data) %in% colnames(model$varprob)))

  if(is.null(numTreesRep)){
    numTreesRep <- nTree
  }

  # null model fuunction
  permuteBART <- function(data){

    yPerm <- sample(data[, responseIdx], replace = FALSE)
    x <- data[, -responseIdx]

    bmodelPerm <-  wbart(x.train = x,
                         y.train = yPerm,
                         nskip = burnIn,
                         ndpost = nMCMC, # MCMC iters
                         nkeeptreedraws = nMCMC,
                         ntree = numTreesRep
    )

    varPropsPerm <- bmodelPerm$varcount
    varPropsPermAvg <- colMeans(proportions(varPropsPerm, 1))
    return(varPropsPermAvg)
  }

  for (i in 1:numRep) {
    permuteMat[i, ] = permuteBART(data)
  }


  permuteMat <- permuteMat[, names(varPropAvg)]
  Cutoffs <- apply(permuteMat, 2, quantile, probs = 1 - alpha)

  vimpName <- names(varPropAvg[varPropAvg > Cutoffs & varPropAvg > 0])

  vimpColNum <- sapply(1:length(vimpName), function(x){
    which(vimpName[x] == colnames(model$varprob))
    })

  # get metrics
  permSE = apply(permuteMat, 2, sd)/sqrt(nrow(permuteMat))
  permAvg = apply(permuteMat, 2, mean)
  maxCut = quantile(apply(permuteMat, 1, max), 1 - alpha)

  vimpIdx = which(varPropAvg > 0)[1:min(10, length(which(varPropAvg > 0)))]

  localThresholdsDF <- dplyr::tibble(
    Variable = names(Cutoffs),
    lThres = unname(Cutoffs)
  )

  incProp <- dplyr::tibble(
    Variable = names(varPropAvg),
    imp = unname(varPropAvg)
  )

  # reorder
  localThresholdsDF <- localThresholdsDF[ order(match(localThresholdsDF$Variable, incProp$Variable)), ]


  localThresholdsDF$Variable <- factor(localThresholdsDF$Variable, levels = names(varPropAvg))
  incProp$Variable <- factor(incProp$Variable, levels = names(varPropAvg))

  incProp$shape <- ifelse(incProp$imp > localThresholdsDF$lThres, 19, 4)
  incProp$threshold <- localThresholdsDF$lThres

  # add shift
  incProp$difference <- incProp$imp - incProp$threshold
  incProp$difference[incProp$difference <=  0] <- 0

  # for(i in seq_along(incProp$Variable)){
  #   incProp$zSc[i] <- (incProp$imp[i] - mean(incProp$imp))/sd(incProp$imp)
  # }

  incProp$Variable <- factor(incProp$Variable, levels = rev(incProp$Variable))


  if(shift){
    p <- ggplot(incProp, aes(x = Variable, y = difference)) +
      geom_point(size = 3) +
      theme_bw() + ylab('proportion included') + coord_flip()
  }else{
     p <-  ggplot(incProp, aes(x = Variable, y = threshold)) +
        geom_segment(aes(x=Variable, xend=Variable, y=0, yend=threshold), col = 'steelblue') +
        geom_point(aes(x = Variable, y = imp), shape = incProp$shape, size = 3) +
        theme_bw() + ylab('proportion included') + coord_flip()
  }


  return(p)
}


# dbarts ------------------------------------------------------------------

lProd.bart <- function(model, data, numRep = 10, numTreesRep = NULL, alpha = 0.5, shift = FALSE){

  # get some information
  nTree <- model$call$ntree
  nMCMC  <- model$call$ndpost
  nVar  <- as.integer(length(colMeans((model$varcount))))
  varNames <- colnames(model$fit$data@x)
  burnIn <-  model$call$nskip

  # set up matrix
  permuteMat <- matrix(NA, nrow = numRep, ncol = nVar)
  colnames(permuteMat) <- colnames(model$varcount)
  varProp <- model$varcount
  varPropAvg <- colMeans(proportions(varProp, 1))
  varPropAvg <- sort(varPropAvg, decreasing = TRUE)

  responseIdx <- which(!(names(data) %in% colnames(model$varcount)))

  if(is.null(numTreesRep)){
    numTreesRep <- nTree
  }

  # null model fuunction
  permuteDBART <- function(data){

    yPerm <- sample(data[, responseIdx], replace = FALSE)
    x <- data[, -responseIdx]

    bmodelPerm <- bart(x.train = x,
                       y.train = yPerm,
                       ntree = numTreesRep,
                       keeptrees = TRUE,
                       nskip = burnIn,
                       ndpost = nMCMC,
                       combinechains = F,
                       nchain = 1
    )


    varPropsPerm <- bmodelPerm$varcount
    varPropsPermAvg <- colMeans(proportions(varPropsPerm, 1))
    return(varPropsPermAvg)
  }

  for (i in 1:numRep) {
    permuteMat[i, ] = permuteDBART(data)
  }

  permuteMat <- permuteMat[, names(varPropAvg)]
  Cutoffs <- apply(permuteMat, 2, quantile, probs = 1 - alpha)

  vimpName <- names(varPropAvg[varPropAvg > Cutoffs & varPropAvg > 0])

  vimpColNum <- sapply(1:length(vimpName), function(x){
    which(vimpName[x] == colnames(model$varcount))
  })

  # get metrics
  permSE = apply(permuteMat, 2, sd)/sqrt(nrow(permuteMat))
  permAvg = apply(permuteMat, 2, mean)
  maxCut = quantile(apply(permuteMat, 1, max), 1 - alpha)

  vimpIdx = which(varPropAvg > 0)[1:min(10, length(which(varPropAvg > 0)))]

  localThresholdsDF <- dplyr::tibble(
    Variable = names(Cutoffs),
    lThres = unname(Cutoffs)
  )

  incProp <- dplyr::tibble(
    Variable = names(varPropAvg),
    imp = unname(varPropAvg)
  )

  # reorder
  localThresholdsDF <- localThresholdsDF[ order(match(localThresholdsDF$Variable, incProp$Variable)), ]


  localThresholdsDF$Variable <- factor(localThresholdsDF$Variable, levels = names(varPropAvg))
  incProp$Variable <- factor(incProp$Variable, levels = names(varPropAvg))

  incProp$shape <- ifelse(incProp$imp > localThresholdsDF$lThres, 19, 4)
  incProp$threshold <- localThresholdsDF$lThres

  # truncate difference to zero
  incProp$difference <- incProp$imp - incProp$threshold
  incProp$difference[incProp$difference <=  0] <- 0
  incProp$Variable <- factor(incProp$Variable, levels = rev(incProp$Variable))

  if(shift){
    p <- ggplot(incProp, aes(x = Variable, y = difference)) +
      geom_point(size = 3) +
      theme_bw() + ylab('proportion included')
  }else{
    p <-  ggplot(incProp, aes(x = Variable, y = threshold)) +
      geom_segment(aes(x=Variable, xend=Variable, y=0, yend=threshold), col = 'steelblue') +
      geom_point(aes(x = Variable, y = imp), shape = incProp$shape, size = 3) +
      theme_bw() + ylab('proportion included') + coord_flip()
  }

  return(p)
}


# bartMachine -------------------------------------------------------------


lProd.bartMachine <- function(model, data, numRep = 10, numTreesRep = NULL, alpha = 0.5, shift = FALSE){

  # get some information

  nTree <-  model$num_trees
  nMCMC <-  model$num_iterations_after_burn_in
  nVar  <- model$p
  varNames <- colnames(model$X)
  burnIn <-  model$num_burn_in

  # set up matrix
  permuteMat <- matrix(NA, nrow = numRep, ncol = nVar)
  colnames(permuteMat) = model$training_data_features_with_missing_features
  varPropAvg <- bartMachine::get_var_props_over_chain(model)
  varPropAvg <- sort(varPropAvg, decreasing = TRUE)

  responseIdx <- which(!(names(data) %in% colnames(permuteMat)))

  if(is.null(numTreesRep)){
    numTreesRep <- nTree
  }

  # null model fuunction
  permuteBMachine <- function(data){

    yPerm <- sample(data[, responseIdx], replace = FALSE)
    x <- data[, -responseIdx]

    bmodelPerm <- bartMachine(X = x,
                              y = yPerm,
                              num_trees = numTreesRep,
                              flush_indices_to_save_RAM = FALSE,
                              num_burn_in = burnIn,
                              num_iterations_after_burn_in = nMCMC)




    varPropsPermAvg <-  bartMachine::get_var_props_over_chain(bmodelPerm)
    return(varPropsPermAvg)
  }

  for (i in 1:numRep) {
    permuteMat[i, ] = permuteBMachine(data)
  }

  permuteMat <- permuteMat[, names(varPropAvg)]
  Cutoffs <- apply(permuteMat, 2, quantile, probs = 1 - alpha)

  vimpName <- names(varPropAvg[varPropAvg > Cutoffs & varPropAvg > 0])

  vimpColNum <- sapply(1:length(vimpName), function(x){
    which(vimpName[x] == colnames(model$training_data_features_with_missing_features))
  })

  # get metrics
  permSE = apply(permuteMat, 2, sd)/sqrt(nrow(permuteMat))
  permAvg = apply(permuteMat, 2, mean)
  maxCut = quantile(apply(permuteMat, 1, max), 1 - alpha)

  vimpIdx = which(varPropAvg > 0)[1:min(10, length(which(varPropAvg > 0)))]

  localThresholdsDF <- dplyr::tibble(
    Variable = names(Cutoffs),
    lThres = unname(Cutoffs)
  )

  incProp <- dplyr::tibble(
    Variable = names(varPropAvg),
    imp = unname(varPropAvg)
  )

  # reorder
  localThresholdsDF <- localThresholdsDF[ order(match(localThresholdsDF$Variable, incProp$Variable)), ]


  localThresholdsDF$Variable <- factor(localThresholdsDF$Variable, levels = names(varPropAvg))
  incProp$Variable <- factor(incProp$Variable, levels = names(varPropAvg))

  incProp$shape <- ifelse(incProp$imp > localThresholdsDF$lThres, 19, 4)
  incProp$threshold <- localThresholdsDF$lThres

  # truncate difference to zero
  incProp$difference <- incProp$imp - incProp$threshold
  incProp$difference[incProp$difference <=  0] <- 0
  incProp$Variable <- factor(incProp$Variable, levels = rev(incProp$Variable))

  if(shift){
    p <- ggplot(incProp, aes(x = Variable, y = difference)) +
      geom_point(size = 3) +
      theme_bw() + ylab('proportion included')
  }else{
    p <-  ggplot(incProp, aes(x = Variable, y = threshold)) +
      geom_segment(aes(x=Variable, xend=Variable, y=0, yend=threshold), col = 'steelblue') +
      geom_point(aes(x = Variable, y = imp), shape = incProp$shape, size = 3) +
      theme_bw() + ylab('proportion included') + coord_flip()# +
     # geom_hline(yintercept = maxCut, col = 'red')
  }

  return(p)
}




