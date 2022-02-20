#' bartClassDiag
#'
#' @description Displays a selection of diagnostic plots for a BART model.
#'
#' @param model a model created from either the BART, dbarts, or bartMachine package.
#' @param data A dataframe
#' @param response The name of the response for the fit.
#' @param pNorm apply pnorm to the y-hat data
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

bartClassDiag <- function(model, data, response, pNorm = FALSE){

  responseVals <- response


  if(class(model) == "bartMachine"){
    yhatTrain <- model$y_hat_train
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

  threshold <- yI
  class <- as.numeric(yhatTrain > threshold)
  dfFitClassBart$class <- class


  # create data frame for histogram
  dfHist <- data.frame(vals = yhatTrain)
  dfHist$group = ifelse(dfHist$vals < yI, "low", "high")

  # create data fro precision-recall plot
  predRec <- performance(pred, "prec", "rec")

  dfPR <- data.frame(Recall = predRec@x.values[[1]],
                     Precision = predRec@y.values[[1]])




  # -------------------------------------------------------------------------

  ROC <- bartROC(dfROC, threshold = yI, label = aucLab)
  PrecRec <- bartPrecRec(dfPR)
  classF <- classFit(dfFitClassBart, threshold = yI)
  histogram <- classHist(dfHist, threshold = yI)
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


bartPrecRec <- function(data) {

  data <- na.omit(data)
  p <- ggplot(data, aes(x = Recall, y = Precision)) +
    geom_line() +
    xlab("Recall") +
    ylab("Precision") +
    ggtitle("Precision-Recall") +
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
    geom_jitter(height = 0.2, size = 1, alpha = 0.2) +
    scale_color_manual(values = c("red",  "blue")) +
    xlab('Predicted probability') +
    ylab('') +
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
    ggtitle("Histogram") +
    geom_vline(xintercept = threshold, col = 'black') +
    scale_fill_manual(values = c("low" = 'red',
                                 "high" = "blue")) +
    theme_bw() +
    theme(legend.position = "none")

  return(p)

}



# VIMP --------------------------------------------------------------------

bartVimpClass <- function(model) {

  if (class(model) == "pbart") {
    # get variable importance
    vImp <- model$varcount.mean
    vImp <- dplyr::tibble(Variable = names(vImp), Importance = vImp)
  } else if(class(model) == "bartMachine"){
    bmVimp <- bartMachine::investigate_var_importance(model, num_replicates_for_avg = 5, plot = FALSE)
    vimpVals <- bmVimp$avg_var_props
    vImp <- data.frame(Variable = names(vimpVals), Importance = vimpVals, row.names = NULL)
  }else{
    varCount <- NULL
    for(i in 1:length(model$varcount[,1])){
      varCount[[i]] <-  prop.table(model$varcount[i,])
    }
    vImp <- varCount %>%
      bind_rows() %>%
      colMeans()
    vImp <- tibble(Variable = names(vImp), Importance = vImp)
  }

  p <- vImp %>%
    arrange(Importance) %>%
    mutate(Variable = factor(Variable, unique(Variable))) %>%
    ggplot() +
    aes(x = Variable, y = Importance) +
    geom_segment(aes(x = Variable, xend = Variable, y = 0, yend = Importance), color = "blue") +
    geom_point(color = "blue") +
    theme_light() +
    coord_flip() +
    ggtitle(label = "Variable Importance") +
    theme_bw() +
    xlab("Variable") +
    ylab("Importance") +
    theme(
      axis.title.y = element_text(angle = 90, vjust = 0.5),
      legend.key.size = unit(0.5, "cm")
    )

  return(p)
}


# Confusion Matrix --------------------------------------------------------

confMat <- function(model, data, response){

  respIdx <- which(sapply(data, identical, y = response))
  pred <- as.numeric(condvis2::CVpredict(model, data[, -respIdx]))
  pred <-  ifelse(pred == 2, 1, 0)
  pred <- as.factor(pred)
  response <- as.factor(response)


  tab <- table(pred, response)
  acc <- 100 - (mean(pred != response) * 100) # accuracy


  confMatPlot <- function(mat){
    mytitle <- paste("Accuracy: ", acc, "%")

    p <- ggplot(data = as.data.frame(tab), aes(x = response, y = pred)) +
      geom_tile(aes(fill = log(Freq)), colour = "white") +
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



