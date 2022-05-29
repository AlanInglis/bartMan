#' bartClassifDiag
#'
#' @description Displays a selection of diagnostic plots for a BART model.
#'
#' @param model a model created from either the BART, dbarts, or bartMachine package.
#' @param data A dataframe
#' @param response The name of the response for the fit.
#' @param threshold A dashed line on some plots to indicate a chosen threshold value.
#' by default the Youden index is shown.
#' @param pNorm apply pnorm to the y-hat data
#' @param showInterval LOGICAL if TRUE then show 5\% and 95\% quantile intervals.
#'
#' @return A selection of diagnostic plots
#'
#' @import ggplot2
#' @importFrom patchwork area
#' @importFrom patchwork plot_layout
#' @importFrom ROCR prediction
#' @importFrom ROCR performance
#' @importFrom bartMachine investigate_var_importance
#' @importFrom caret confusionMatrix
#' @importFrom condvis2 CVpredict
#'
#' @export

bartClassifDiag <- function(model,
                            data,
                            response,
                            threshold = 'Youden',
                            pNorm = FALSE,
                            showInterval = FALSE){

  responseVals <- response


  if(class(model) == "bartMachine"){
    yhatTrain <- model$y_hat_train
    if(model$pred_type == 'classification'){
      yhatTrain <- as.integer(yhatTrain)-1
    }
  }else{
    yhatTrain <- colMeans(model$yhat.train)
  }

  if(pNorm){
    yhatTrain <- pnorm(yhatTrain)
  }

  # get prediction using ROCR package:
  pred <- prediction(yhatTrain, responseVals)

  # get auc value
  auc <- performance(pred, "auc")
  auc <- auc@y.values[[1]]
  aucLab <- print(paste0("AUC: ", round(auc, 5)))

  # get false/true positive rates
  perfTF <- performance(pred, "tpr", "fpr")

  # create dataframe for ROC plot
  dfROC <- data.frame(fpr = perfTF@x.values[[1]],
                      tpr = perfTF@y.values[[1]])

  # calculate Youden's Index
  youdenIndex <- function(pred) {

    sens <- performance(pred, measure = "sens")@y.values[[1]]
    spec <- performance(pred, measure = "spec")@y.values[[1]]
    #youdenVal <- mean(sens) + mean(spec) - 1
    youdenVal <- max(sens + spec -1)

    return(youdenVal)
  }

  yI <- youdenIndex(pred)

  # create dataframe for fitted vals plot
  dfFitClassBart <- data.frame(fitted = yhatTrain,
                               actual = responseVals)

  if(threshold == 'Youden'){
  threshold <- yI
  }else{
    threshold <- threshold
  }

  class <- as.numeric(yhatTrain > threshold)
  dfFitClassBart$class <- class


  # create data frame for histogram
  dfHist <- data.frame(vals = yhatTrain)
  dfHist$group = ifelse(dfHist$vals < threshold, "low", "high")

  # create data fro precision-recall plot
  predRec <- performance(pred, "prec", "rec")

  dfPR <- data.frame(Recall = predRec@x.values[[1]],
                     Precision = predRec@y.values[[1]])

  # get aucpr value
  aucpr <- performance(pred, "aucpr")
  aucpr <- aucpr@y.values[[1]]
  aucprLab <- print(paste0("AUCPR: ", round(aucpr, 5)))




  # -------------------------------------------------------------------------
  if(showInterval){
    ROC <- rocCI(model = model, response = response, data = data)
    PrecRec <- prCI(model = model, response = response, data = data)
  }else{
    ROC <- bartROC(dfROC, threshold = threshold, label = aucLab)
    PrecRec <- bartPrecRec(dfPR, label = aucprLab)
  }
  classF <- classFit(dfFitClassBart, threshold = threshold)
  histogram <- classHist(dfHist, threshold = threshold)
  vimp <- bartVimpClass(model)
  cM <- confMat(model, data, response)

  design <- c(
    area(1, 1, 3, 3),
    area(1, 5, 3, 7),
    area(5, 1, 7, 3),
    area(5, 5, 7, 7),
    area(9, 1, 11, 3),
    area(9, 5, 11, 7)
  )

  diagPlot <- ROC + PrecRec + classF + histogram + cM + vimp + plot_layout(design = design)

  return(diagPlot)



}
# ROC curve for dbarts ----------------------------------------------------


bartROC <- function(data, threshold, label){

  p <- ggplot(data, aes(x = fpr, y = tpr)) +
    geom_line() +
    xlab('False positive rate') +
    ylab('True positive rate') +
    ggtitle('ROC') +
    #geom_abline(intercept = 0, slope = 1, col = 'blue') +
    geom_segment(aes(x = threshold, xend = threshold, y = -Inf,  yend = 1), linetype = 2) +
    geom_segment(aes(x = -Inf, xend = threshold, y = 1, yend = 1),  linetype = 2) +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), col = 'blue') +
    xlim(0, 1) +
    annotate("text", x = 0.75, y = 0.1, label = label, size = 3) +
    #geom_vline(xintercept = threshold, col = 'black') +
    theme_bw()

  return(p)
}



# Precision-Recall --------------------------------------------------------


bartPrecRec <- function(data, label) {

  data <- na.omit(data)
  p <- ggplot(data, aes(x = Recall, y = Precision)) +
    geom_line() +
    xlab("Recall") +
    ylab("Precision") +
    ggtitle("Precision-Recall") +
    xlim(0, 1) +
    annotate("text", x = 0.75, y = 0.56, label = label, size = 3) +
    theme_bw()

  return(p)
}


# Class fitted ------------------------------------------------------------


classFit <- function(data, threshold){

  p <- ggplot(data, aes(x = fitted,
                        y = factor(actual),
                        fill = factor(class),
                        col = factor(class))
  ) +
    #geom_point(alpha = 0.2) +
    geom_jitter(height = 0.2, size = 1, alpha = 0.2,  na.rm=TRUE) +
    scale_color_manual(values = c("red",  "blue")) +
    xlab('Predicted probability') +
    ylab('') +
    xlim(0, 1) +
    ggtitle('Fitted values') +
    #xlim(0,1) +
    theme_bw() +
    theme(legend.position = 'none') +
    geom_vline(xintercept = threshold, col = 'black')

  return(p)

}

# Histogram fitted vals ---------------------------------------------------


classHist <- function(data, threshold){

  p <- ggplot(data, aes(vals, fill = group)) +
    geom_histogram(bins = 50, color = "black", alpha = 0.5) +
    ylab('') +
    xlab('Predicted probability') +
    expand_limits(x = 1) +
    ggtitle("Histogram") +
    geom_vline(xintercept = threshold, col = 'black') +
    scale_fill_manual(values = c("low" = 'red',
                                 "high" = "blue")) +
    theme_bw() +
    theme(legend.position = "none")

  return(p)

}



# VIMP --------------------------------------------------------------------

bartVimpClass <- function(model){


  if(class(model) == "bartMachine"){
    vImp <- bartMachine::get_var_counts_over_chain(model)
  } else {
    # get variable importance
    vImp <- model$varcount
  }

  vImpProps <- proportions(vImp, 1)
  vImp <- colMeans(vImpProps)


  # get quantiles of proportions
  vimp25 <- apply(vImpProps, 2, function(x) quantile(x, c(.25)))
  vimp75 <- apply(vImpProps, 2, function(x) quantile(x, c(.75)))

  vImp <- dplyr::tibble(
    Variable = names(vImp),
    imp = vImp,
    upperQ = vimp75,
    lowerQ = vimp25
  )

  p <- vImp %>%
    arrange(imp) %>%
    mutate(Variable = factor(Variable, unique(Variable))) %>%
    ggplot(aes(x = Variable, y = imp)) +
    ggforce::geom_link(aes(
      x = Variable, xend = Variable, yend = upperQ,
      col = Variable, alpha = rev(stat(index))
    ),
    size = 2, n = 1000
    ) +
    ggforce::geom_link(aes(
      x = Variable, xend = Variable, yend = lowerQ,
      col = Variable, alpha = rev(stat(index))
    ),
    size = 2, n = 1000
    ) +
    geom_point(aes(x = Variable, y = imp), shape = 18, size = 2, color = "black") +
    coord_flip() +
    theme_bw() +
    labs(x = "Variable", y = "Importance") +
    ggtitle('VImp') +
    theme(legend.position = "none")

  return(p)
}




# Confusion Matrix --------------------------------------------------------

confMat <- function(model, data, response){

  respIdx <- which(sapply(data, identical, y = response))
  pred <- round(as.numeric(condvis2::CVpredict(model, data[, -respIdx])), 0)

  if(class(model) == 'pbart' || class(model) == 'bart'){
    pred <-  ifelse(pred == 2, 1, 0)
  }

  pred <- as.factor(pred)
  if(class(model) == 'bartMachine'){
    response <- as.numeric(response)
    #response <-  ifelse(response == 2, 1, 0)
    #pred <-  ifelse(pred == 2, 1, 0)
  }else{
    response <- as.numeric(as.factor(response))
    response <-  ifelse(response == 2, 1, 0)
  }


  tab <- table(pred, response)
  acc <- 100 - (mean(pred != response) * 100) # accuracy
  acc <- round(acc, 2)


  confMatPlot <- function(mat){
    mytitle <- paste("Accuracy: ", acc, "%")

    p <- ggplot(data = as.data.frame(tab), aes(x = response, y = pred)) +
      geom_tile(aes(fill = Freq)) +
      scale_fill_gradient(low = "white", high = "steelblue") +
      geom_text(aes(x = response, y = pred, label = Freq)) +
      ylab('Prediciton') +
      xlab("Response") +
      theme_bw() +
      theme(legend.position = "none") +
      ggtitle(mytitle)
    return(p)
  }
  confMatPlot(cfm)
}


# ROC with CI -------------------------------------------------------------

rocCI <- function(model, response, data){

  # get model info
  if(class(model) == 'pbart' | class(model) == 'wbart'){

    modelTrees <- model$treedraws$trees
    modelInfo <- unlist(strsplit(modelTrees, " "))[1:3]
    modelInfo <- gsub("(^\\d+)([\a-zA-Z0-9]*)", "\\1", modelInfo)
    nMCMC <- as.integer(modelInfo[1])
    nTree <- as.integer(modelInfo[2])
    nVar <- as.integer(modelInfo[3])
    burnIn <- length(model$sigma) - nMCMC
    varNames <- names(model$varcount.mean)

  }else if(class(model) == 'bart'){

    nTree <- model$call$ntree
    nMCMC  <- model$call$ndpost
    nVar  <- as.integer(length(colMeans((model$varcount))))
    varNames <- colnames(model$fit$data@x)
    burnIn <-  model$call$nskip

  }else if(class(model) == 'bartMachine'){
    nTree <-  model$num_trees
    nMCMC <-  model$num_iterations_after_burn_in
    nVar  <- model$p
    varNames <- colnames(model$X)
    burnIn <-  model$num_burn_in
  }

  # get yhats
  responseIdx <- which(!(names(data) %in% varNames))
  if(class(model) == 'bartMachine'){
    yhatTrain = bartMachine::bart_machine_get_posterior(model, data[, -responseIdx])$y_hat_posterior_samples
  }else{
    yhatTrain = model$yhat.train
  }

  # small set up
  pred = NULL
  perfTF = NULL
  dfROC = NULL
  noIter = nMCMC


  # get roc predicitons
  for(i in 1:noIter){
    if(class(model) == 'bartMachine'){
      pred[[i]] <- ROCR::prediction(yhatTrain[,i], response)
    }else{
      pred[[i]] <- ROCR::prediction(yhatTrain[i,], response)
    }
    perfTF[[i]] <- ROCR::performance(pred[[i]], "tpr", "fpr")
    dfROC[[i]] <- data.frame(fpr = perfTF[[i]]@x.values[[1]],
                             tpr = perfTF[[i]]@y.values[[1]])
    dfROC[[i]]$n <- c(1:nrow(dfROC[[i]]))
  }

  # put together in single df
  newDF <- lapply(seq_along(dfROC),
                  function(x) {
                    dfROC[[x]] <- dfROC[[x]] %>%
                      mutate(id1 = 1:n(), id2 = x)
                  }) %>%
    bind_rows() %>%
    arrange(id1, id2) %>%
    select(-id1, -id2) %>%
    as.data.frame()


  # split into lists and remove var n
  dfList <- split(newDF, newDF$n)
  for(i in 1:length(dfList)){
    dfList[[i]]$n <- NULL
  }

  # get quantiles
  lo <- lapply(dfList, function(x) apply(x, 2, function(x) quantile(x, c(.05))))
  hi <- lapply(dfList, function(x) apply(x, 2, function(x) quantile(x, c(.95))))
  lowQROC <- bind_rows(lo)
  hiQROC <- bind_rows(hi)

  # get median
  me <-  lapply(dfList, function(x) apply(x, 2, function(x) median(x)))
  medianDF <- bind_rows(me)

  # unify into single df
  joinDF <- data.frame(
    tpr = medianDF$tpr,
    fpr = medianDF$fpr,
    hiTPR = hiQROC$tpr,
    hiFPR = hiQROC$fpr,
    lowTPR = lowQROC$tpr,
    lowFPR = lowQROC$fpr
  )

  # plot
  p <- ggplot(joinDF, aes(x = fpr, y = tpr)) +
    geom_line() +
    geom_ribbon(aes(ymax = hiTPR, ymin = lowTPR),
                fill="steelblue", alpha=.5) +
    xlab('False positive rate') +
    ylab('True positive rate') +
    ggtitle('ROC') +
    theme_bw()

  return(p)
}


prCI <- function(model, response, data){

  # get model info
  if(class(model) == 'pbart' | class(model) == 'wbart'){

    modelTrees <- model$treedraws$trees
    modelInfo <- unlist(strsplit(modelTrees, " "))[1:3]
    modelInfo <- gsub("(^\\d+)([\a-zA-Z0-9]*)", "\\1", modelInfo)
    nMCMC <- as.integer(modelInfo[1])
    nTree <- as.integer(modelInfo[2])
    nVar <- as.integer(modelInfo[3])
    burnIn <- length(model$sigma) - nMCMC
    varNames <- names(model$varcount.mean)

  }else if(class(model) == 'bart'){

    nTree <- model$call$ntree
    nMCMC  <- model$call$ndpost
    nVar  <- as.integer(length(colMeans((model$varcount))))
    varNames <- colnames(model$fit$data@x)
    burnIn <-  model$call$nskip

  }else if(class(model) == 'bartMachine'){
    nTree <-  model$num_trees
    nMCMC <-  model$num_iterations_after_burn_in
    nVar  <- model$p
    varNames <- colnames(model$X)
    burnIn <-  model$num_burn_in
  }

  # get yhats
  responseIdx <- which(!(names(data) %in% varNames))
  if(class(model) == 'bartMachine'){
    yhatTrain = bartMachine::bart_machine_get_posterior(model, data[, -responseIdx])$y_hat_posterior_samples
  }else{
    yhatTrain = model$yhat.train
  }

  # small set up
  pred = NULL
  predPR = NULL
  dfPR = NULL
  noIter = nMCMC

  # get roc predicitons
  for(i in 1:noIter){
    if(class(model) == 'bartMachine'){
      pred[[i]] <- ROCR::prediction(yhatTrain[,i], response)
    }else{
      pred[[i]] <- ROCR::prediction(yhatTrain[i,], response)
    }
    predPR[[i]] <- ROCR::performance(pred[[i]], "prec", "rec")
    dfPR[[i]] <- data.frame(rec = predPR[[i]]@x.values[[1]],
                            prec  = predPR[[i]]@y.values[[1]])
    dfPR[[i]]$n <- c(1:nrow(dfPR[[i]]))
  }

  # remove NaN
  # dfPR <- lapply(dfPR, function(dat) {
  #   dat[] <- lapply(dat, function(x) replace(x, is.nan(x), 0))
  # })

  # put together in single df
  newDF <- lapply(seq_along(dfPR),
                  function(x) {
                    dfPR[[x]] <- dfPR[[x]] %>%
                      mutate(id1 = 1:n(), id2 = x)
                  }) %>%
    bind_rows() %>%
    arrange(id1, id2) %>%
    select(-id1, -id2) %>%
    as.data.frame()

  # remove NaN
  newDF <- newDF[(nMCMC+1):(length(newDF[,1])),]

  # split into lists and remove var n
  dfList <- split(newDF, newDF$n)
  for(i in 1:length(dfList)){
    dfList[[i]]$n <- NULL
  }

  # get quantiles
  lo <- lapply(dfList, function(x) apply(x, 2, function(x) quantile(x, c(.05), na.rm = T)))
  hi <- lapply(dfList, function(x) apply(x, 2, function(x) quantile(x, c(.95), na.rm = T)))
  lowQROC <- bind_rows(lo)
  hiQROC <- bind_rows(hi)


  # get median
  me <-  lapply(dfList, function(x) apply(x, 2, function(x) median(x)))
  medianDF <- bind_rows(me)


  # unify into single df
  joinDF <- data.frame(
    prec = medianDF$prec,
    rec = medianDF$rec,
    hiprec = hiQROC$prec,
    hirec = hiQROC$rec,
    lowprec = lowQROC$prec,
    lowrec = lowQROC$rec
  )

  #joinDF <- joinDF[-1,]

  # plot
  p <- ggplot(joinDF, aes(x = rec, y = prec)) +
    geom_line() +
    geom_ribbon(aes(ymax = hiprec, ymin = lowprec),
                fill="steelblue", alpha=.5) +
    scale_x_continuous(limits = c(min(joinDF$rec), max(joinDF$rec))) +
    xlab("Recall") +
    ylab("Precision") +
    ggtitle("Precision-Recall") +
    theme_bw()


  return(p)

}



