#' bartDiag
#'
#' @description Displays a selection of diagnostic plots for a BART model.
#'
#' @param model a model created from either the BART, modelarts, or bartMachine package.
#' @param response The name of the response for the fit.
#' @param burnIn Trace plot will only show iterations above selected burn in value.
#' @param data A dataframe used to build the model.
#' @param threshold A dashed line on some plots to indicate a chosen threshold value (classification only).
#' by default the Youden index is shown.
#' @param pNorm apply pnorm to the y-hat data (classification only).
#' @param showInterval LOGICAL if TRUE then show 5\% and 95\% quantile intervals on ROC an PC curves (classification only).
#' @param combineFactors Whether or not to combine dummy variables (if present) in display.
#'
#'
#' @return A selection of diagnostic plots.
#'
#' @import tidytreatment
#' @importFrom patchwork area
#' @importFrom patchwork plot_layout
#' @importFrom tidybayes residual_draws
#' @importFrom tidybayes point_interval
#' @importFrom tidybayes geom_pointinterval
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr summarise
#' @importFrom dplyr tibble
#' @import ggplot2
#'
#' @examples
#' \donttest{
#' # For Regression
#' # Generate Friedman data
#' fData <- function(n = 200, sigma = 1.0, seed = 1701, nvar = 5) {
#'   set.seed(seed)
#'   x <- matrix(runif(n * nvar), n, nvar)
#'   colnames(x) <- paste0("x", 1:nvar)
#'   Ey <- 10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 + 10 * x[, 4] + 5 * x[, 5]
#'   y <- rnorm(n, Ey, sigma)
#'   data <- as.data.frame(cbind(x, y))
#'   return(data)
#' }
#' f_data <- fData(nvar = 10)
#' x <- f_data[, 1:10]
#' y <- f_data$y
#'
#' # Create dbarts model
#' library(dbarts)
#' set.seed(1701)
#' dbartModel <- bart(x, y, ntree = 5, keeptrees = TRUE, nskip = 10, ndpost = 10)
#'
#' bartDiag(model = dbartModel, response = "y", burnIn = 100, data = f_data)
#'
#'
#' # For Classification
#' data(iris)
#' iris2 <- iris[51:150, ]
#' iris2$Species <- factor(iris2$Species)
#'
#' # Create dbarts model
#' dbartModel <- bart(iris2[, 1:4], iris2[, 5], ntree = 5, keeptrees = TRUE, nskip = 10, ndpost = 10)
#'
#' bartDiag(model = dbartModel, data = iris2, response = iris2$Species)
#' }
#'
#' @export

bartDiag <- function(model,
                      data,
                      response,
                      burnIn = 0,
                      threshold = 'Youden',
                      pNorm = FALSE,
                      showInterval = TRUE,
                      combineFactors = FALSE){

  if(inherits(model, 'pbart')){

    output <- bartClassifDiag(model = model, data = data, response = response,
                              threshold = threshold, pNorm = pNorm, showInterval = showInterval,
                              combineFactors = combineFactors)

  }else if(inherits(model, 'wbart')){


    output <- bartRegrDiag(model = model, data = data, response = response,
                           burnIn = burnIn, combineFactors = combineFactors)

  }else if(inherits(model, 'bart')){

    if(model$fit$control@binary == TRUE){

      output <- bartClassifDiag(model = model, data = data, response = response,
                                threshold = threshold, pNorm = pNorm, showInterval = showInterval,
                                combineFactors = combineFactors)

    }else{
      output <- bartRegrDiag(model = model, data = data, response = response,
                             burnIn = burnIn, combineFactors = combineFactors)
    }

  }else if(inherits(model, 'bartMachine')){

    if(model$pred_type == 'classification'){

      output <- bartClassifDiag(model = model, data = data, response = response,
                                threshold = threshold, pNorm = pNorm, showInterval = showInterval,
                                combineFactors = combineFactors)
    }else{
      output <- bartRegrDiag(model = model, data = data, response = response,
                             burnIn = burnIn, combineFactors = combineFactors)
    }
  }
  return(output)
}


# -------------------------------------------------------------------------
# Diagnostics for regression ----------------------------------------------
# -------------------------------------------------------------------------

#' bartRegrDiag
#'
#' @description Displays a selection of diagnostic plots for a BART model.
#'
#' @param model a model created from either the BART, modelarts, or bartMachine package.
#' @param response The name of the response for the fit.
#' @param burnIn Trace plot will only show iterations above selected burn in value.
#' @param data A dataframe used to build the model.
#' @param combineFactors Whether or not to combine dummy variables (if present) in display.
#'
#'
#' @return A selection of diagnostic plots
#'

#' @import ggplot2
#' @import tidytreatment
#' @importFrom patchwork area
#' @importFrom patchwork plot_layout
#' @importFrom tidybayes residual_draws
#' @importFrom tidybayes point_interval
#' @importFrom tidybayes geom_pointinterval
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr summarise
#' @importFrom dplyr tibble
#' @export


bartRegrDiag <- function(model,
                         response,
                         burnIn = 0,
                         data,
                         combineFactors = FALSE) {

  qq <- bartQQ(model, response)
  trace <- bartTrace(model, burnIn = burnIn)
  residual <- bartResiduals(model, response = response)
  histogram <- bartHist(model, response)
  fitVSact <- bartFitted(model, response)
  vImp <- bartVimp(model, data = data, combineFactors = combineFactors)

  design <- c(
    area(1, 1, 3, 3),
    area(1, 5, 3, 7),
    area(5, 1, 7, 3),
    area(5, 5, 7, 7),
    area(9, 1, 11, 3),
    area(9, 5, 11, 7)
  )

  diagPlot <- qq + trace + residual + histogram +
    fitVSact + vImp + plot_layout(design = design)

  return(diagPlot)
}


# QQ plot -----------------------------------------------------------------

bartQQ <- function(model, response) {
  if (inherits(model, "wbart") || inherits(model, "bartMachine")){
    res <- tidybayes::residual_draws(model, response = response, include_newdata = FALSE)
    res <- res %>% summarise(.residual = mean(.residual))
  } else {
    res <- data.frame(.residual = residuals(model))
  }

  p <- res %>%
    ggplot(aes(sample = .residual)) +
    geom_qq(color = "blue", alpha = 0.2) +
    geom_qq_line() +
    xlab("Theoretical") +
    ylab("Sample") +
    theme_bw() +
    ggtitle("Q-Q plot")

  return(p)
}





# Trace plot --------------------------------------------------------------


bartTrace <- function(model, burnIn = 0) {
  if (inherits(model, "wbart") || inherits(model, "bartMachine")) {
    # get values
    varDraws <- tidytreatment::variance_draws(model, value = "siqsq")
    varDraws$sigma <- sqrt(varDraws$siqsq)
  } else {
    varDraws <- data.frame(
      sigma = model$sigma,
      .draw = c(1:model$call$ndpost)
    )
  }


  maxIter <- max(varDraws$.draw)
  if (burnIn > 0 && burnIn < maxIter) {
    # Plot includes all parts because there's data after burnIn
    p <- ggplot(varDraws[1:burnIn, ], aes(x = .draw, y = sigma)) +
      geom_vline(xintercept = burnIn, linetype = 5, alpha = 0.5) +
      geom_line(alpha = 0.5) +
      geom_line(data = varDraws[burnIn:maxIter, ], aes(x = .draw, y = sigma), colour = "blue", alpha = 0.5) +
      theme_bw() +
      xlab("Iteration") +
      ylab("Sigma") +
      ggtitle("Trace plot")
  } else {
    # Plot only the initial segment because burnIn is 0 or equals maxIter
    p <- ggplot(varDraws, aes(x = .draw, y = sigma)) +
      geom_line(alpha = 0.5) +
      theme_bw() +
      xlab("Iteration") +
      ylab("Sigma") +
      ggtitle("Trace plot")
  }

  return(p)
}





# Residuals vs Fitted --------------------------------------------------------------


bartResiduals <- function(model,
                          response) {
  if (inherits(model, "wbart") || inherits(model, "bartMachine")){
    res <- tidybayes::residual_draws(model, response = response, include_newdata = FALSE)
  } else {
    res <- data.frame(
      .residual = residuals(model),
      .fitted = fitted(model)
    )
  }


  if (inherits(model, "wbart") || inherits(model, "bartMachine")){
    p <- res %>%
      tidybayes::point_interval(.fitted, .residual, .width = c(0.95)) %>%
      select(-.fitted.lower, -.fitted.upper) %>%
      ggplot() +
      tidybayes::geom_pointinterval(aes(x = .fitted, y = .residual, ymin = .residual.lower, ymax = .residual.upper),
                                    alpha = 0.1,
                                    color = "blue"
      ) +
      theme_bw() +
      xlab("Fitted") +
      ylab("Residual") +
      ggtitle("Fitted vs Residuals")
  } else {
    p <- ggplot(res, aes(.fitted, .residual)) +
      geom_point(color = "blue", alpha = 0.1) +
      theme_bw() +
      xlab("Fitted") +
      ylab("Residual") +
      ggtitle("Fitted vs Residuals")
  }

  return(p)
}




# Histogram Residuals -----------------------------------------------------

bartHist <- function(model, response) {
  if (inherits(model, "wbart") || inherits(model, "bartMachine")) {
    res <- tidybayes::residual_draws(model, response = response, include_newdata = FALSE)
  } else {
    res <- data.frame(.residual = residuals(model))
  }

  p <- ggplot(data = res, aes(.residual)) +
    geom_histogram(bins = 50, color = "blue", fill = "white") +
    theme_bw() +
    xlab("Residual") +
    ylab("Frequency") +
    ggtitle("Histogram")

  return(p)
}



# Fitted Vs Actual --------------------------------------------------------


bartFitted <- function(model, response) {

  if (inherits(model, "wbart") || inherits(model, "bartMachine")) {
    res <- tidybayes::residual_draws(model, response = response, include_newdata = FALSE)
  } else {
    plquants = c(.05,.95)
    cols = c('blue', 'black')

    qLow <- apply(model$yhat.train, length(dim(model$yhat.train)), quantile, probs = plquants[1])
    qMid <- apply(model$yhat.train, length(dim(model$yhat.train)), quantile, probs = 0.5)
    qUpp <- apply(model$yhat.train, length(dim(model$yhat.train)), quantile, probs=plquants[2])

    res <- data.frame(y = response,
                      qMid = qMid,
                      qLow = qLow,
                      qUpp = qUpp)
  }


  if (inherits(model, "wbart") || inherits(model, "bartMachine")) {
    p <- res %>%
      tidybayes::point_interval(.fitted, y, .width = c(0.95)) %>%
      select(-y.lower, -y.upper) %>%
      ggplot() +
      tidybayes::geom_pointinterval(aes(x = y, y = .fitted, ymin = .fitted.lower, ymax = .fitted.upper),
                                    alpha = 0.1,
                                    color = "blue"
      ) +
      geom_smooth(aes(x = y, y = .fitted), method = "lm", color = "black", formula = y ~ x) +
      xlab("Actual") +
      ylab("Fitted") +
      theme_bw() +
      ggtitle("Actual vs Fitted")

  }else{
    p <-  ggplot(data = res, aes(y, qMid)) +
      geom_point(color = "blue", alpha = 0.1) +
      theme_bw() +
      xlab("Actual") +
      ylab("Fitted") +
      ggtitle("Actual vs Fitted") +
      geom_smooth(aes(x = y, y = qMid), method = "lm", color = "black", formula = y ~ x) +
      tidybayes::geom_pointinterval(aes(x = y, y = qMid, ymin = qLow, ymax = qUpp),
                                    alpha = 0.1,
                                    color = "blue"
      )
  }


  return(p)
}



# Variable Importance -----------------------------------------------------

bartVimp <- function(model,  data, combineFactors = FALSE) {

  if (inherits(model, "bartMachine")){
    vImp <- bartMachine::get_var_counts_over_chain(model)
  } else {
    # get variable importance
    vImp <- model$varcount
  }

  if(combineFactors){
    vImp <-  combineDummyDiag(data = data, vimp = vImp)
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
      colour = "gray50", alpha = rev(after_stat(index))
    ),
    linewidth = 2, n = 1000
    ) +
    ggforce::geom_link(aes(
      x = Variable, xend = Variable, yend = lowerQ,
      colour = "gray50", alpha = rev(after_stat(index))
    ),
    linewidth = 2, n = 1000
    ) +
    geom_point(aes(x = Variable, y = imp), shape = 18, size = 2, color = "black") +
    scale_colour_identity() +
    coord_flip() +
    theme_bw() +
    labs(x = "Variable", y = "Importance") +
    theme(legend.position = "none")

  return(p)

}


# -------------------------------------------------------------------------
# Diagnostic for classification -------------------------------------------
# -------------------------------------------------------------------------

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
#' @param combineFactors Whether or not to combine dummy variables (if present) in display.
#'
#' @return A selection of diagnostic plots
#'
#' @import ggplot2
#' @importFrom patchwork area
#' @importFrom patchwork plot_layout
#' @importFrom stats predict
#'
#' @export

bartClassifDiag <- function(model,
                            data,
                            response,
                            threshold = 'Youden',
                            pNorm = FALSE,
                            showInterval = TRUE,
                            combineFactors = FALSE){


  if (!requireNamespace("ROCR", quietly = TRUE)) {
    stop("Package \"ROCR\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  responseVals <- response

  if (inherits(model, "bartMachine")){
    yhatTrain <- model$y_hat_train
    if(model$pred_type == 'classification'){
      yhatTrain <- as.integer(yhatTrain)-1
    }
  }else{
    yhatTrain <- colMeans(model$yhat.train)
  }

  if(pNorm){
    yhatTrain <- stats::pnorm(yhatTrain)
  }

  # get prediction using ROCR package:
  pred <- ROCR::prediction(yhatTrain, responseVals)

  # get auc value
  auc <- ROCR::performance(pred, "auc")
  auc <- auc@y.values[[1]]
  aucLab <- message(paste0("AUC: ", round(auc, 5)))

  # get false/true positive rates
  perfTF <- ROCR::performance(pred, "tpr", "fpr")

  # create dataframe for ROC plot
  dfROC <- data.frame(fpr = perfTF@x.values[[1]],
                      tpr = perfTF@y.values[[1]])

  # calculate Youden's Index
  youdenIndex <- function(pred) {

    sens <- ROCR::performance(pred, measure = "sens")@y.values[[1]]
    spec <- ROCR::performance(pred, measure = "spec")@y.values[[1]]
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
  predRec <- ROCR::performance(pred, "prec", "rec")

  dfPR <- data.frame(Recall = predRec@x.values[[1]],
                     Precision = predRec@y.values[[1]])

  # get aucpr value
  aucpr <- ROCR::performance(pred, "aucpr")
  aucpr <- aucpr@y.values[[1]]
  aucprLab <- message(paste0("AUCPR: ", round(aucpr, 5)))




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
  vimp <- bartVimpClass(model, data = data)
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

  data <- stats::na.omit(data)
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

bartVimpClass <- function(model, data, combineFactors = FALSE){



  if (inherits(model, "bartMachine")){
    vImp <- bartMachine::get_var_counts_over_chain(model)
  } else {
    # get variable importance
    vImp <- model$varcount
  }

  if(combineFactors){
    vImp <-  combineDummyDiag(data = data, vimp = vImp)
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
      colour = "gray50", alpha = rev(after_stat(index))
    ),
    size = 2, n = 1000
    ) +
    ggforce::geom_link(aes(
      x = Variable, xend = Variable, yend = lowerQ,
      colour = "gray50", alpha = rev(after_stat(index))
    ),
    size = 2, n = 1000
    ) +
    geom_point(aes(x = Variable, y = imp), shape = 18, size = 2, color = "black") +
    scale_colour_identity() +
    coord_flip() +
    theme_bw() +
    labs(x = "Variable", y = "Importance") +
    theme(legend.position = "none")

  return(p)
}




# Confusion Matrix --------------------------------------------------------

confMat <- function(model, data, response){

  respIdx <- which(sapply(data, identical, y = response))
  if (inherits(model, 'bart') || inherits(model, 'pbart') || inherits(model, 'wbart')){
    pred <-  round(colMeans(stats::predict(model, data[, -respIdx])), 0)
  }else{
    pred <- round(stats::predict(model, data[, -respIdx]), 0)
  }


  if (inherits(model, 'pbart') || inherits(model, 'bart')){
    pred <-  ifelse(pred == 2, 1, 0)
  }

  pred <- as.factor(pred)
  if (inherits(model, 'bartMachine')){
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
  if (inherits(model, 'pbart') || inherits(model, 'wbart')){

    modelTrees <- model$treedraws$trees
    modelInfo <- unlist(strsplit(modelTrees, " "))[1:3]
    modelInfo <- gsub("(^\\d+)([\a-zA-Z0-9]*)", "\\1", modelInfo)
    nMCMC <- as.integer(modelInfo[1])
    nTree <- as.integer(modelInfo[2])
    nVar <- as.integer(modelInfo[3])
    burnIn <- length(model$sigma) - nMCMC
    varNames <- names(model$varcount.mean)

  }else if (inherits(model, 'bart')){

    nTree <- model$call$ntree
    nMCMC  <- model$call$ndpost
    if(is.null(nMCMC)){
      nMCMC <- 1000
    }
    nVar  <- as.integer(length(colMeans((model$varcount))))
    varNames <- colnames(model$fit$data@x)
    burnIn <-  model$call$nskip

  }else if (inherits(model, 'bartMachine')){
    nTree <-  model$num_trees
    nMCMC <-  model$num_iterations_after_burn_in
    nVar  <- model$p
    varNames <- colnames(model$X)
    burnIn <-  model$num_burn_in
  }

  # get yhats
  responseIdx <- which(!(names(data) %in% varNames))
  if (inherits(model, 'bartMachine')){
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
    if (inherits(model, 'bartMachine')){
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
  if (inherits(model, 'pbart') || inherits(model, 'wbart')){

    modelTrees <- model$treedraws$trees
    modelInfo <- unlist(strsplit(modelTrees, " "))[1:3]
    modelInfo <- gsub("(^\\d+)([\a-zA-Z0-9]*)", "\\1", modelInfo)
    nMCMC <- as.integer(modelInfo[1])
    nTree <- as.integer(modelInfo[2])
    nVar <- as.integer(modelInfo[3])
    burnIn <- length(model$sigma) - nMCMC
    varNames <- names(model$varcount.mean)

  }else if (inherits(model, 'bart')){

    nTree <- model$call$ntree
    nMCMC  <- model$call$ndpost
    if(is.null(nMCMC)){
      nMCMC <- 1000
    }
    nVar  <- as.integer(length(colMeans((model$varcount))))
    varNames <- colnames(model$fit$data@x)
    burnIn <-  model$call$nskip

  }else if (inherits(model, 'bartMachine')){
    nTree <-  model$num_trees
    nMCMC <-  model$num_iterations_after_burn_in
    nVar  <- model$p
    varNames <- colnames(model$X)
    burnIn <-  model$num_burn_in
  }

  # get yhats
  responseIdx <- which(!(names(data) %in% varNames))
  if (inherits(model, 'bartMachine')){
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
    if (inherits(model, 'bartMachine')){
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


# -------------------------------------------------------------------------


# combine factors ---------------------------------------------------------

combineDummyDiag <- function(data, vimp) {
  # Identify factor columns in trees$data
  dfOG <- data
  factorColNam <- names(which(!(sapply(dfOG[colnames(dfOG)], is.numeric))))
  factorCols <- which((colnames(dfOG) %in% factorColNam))

  # Create a list of the variables split into their factors
  facLevelsList <- list()
  for (i in 1:length(factorCols)) {
    facLevels <- unique(dfOG[, factorCols[i]])
    facLevelsList[[names(dfOG)[factorCols[i]]]] <- as.character(facLevels)
  }

  # Function to find original factor name for dummy variables
  find_original_factor_name <- function(varName, factorList) {
    for (factorName in names(factorList)) {
      if (grepl(factorName, varName)) {
        return(factorName)
      }
    }
    return(varName) # Return original if no match found
  }

  # Update entries in trees$structure$var
  colnames(vimp) <- sapply(colnames(vimp), find_original_factor_name, factorList = facLevelsList)
  names(colnames(vimp)) <- NULL

  vimp <- as.data.frame(vimp)
  vimp <- sapply(split.default(vimp, colnames(vimp)), rowSums, na.rm = TRUE)

  return(vimp)
}




