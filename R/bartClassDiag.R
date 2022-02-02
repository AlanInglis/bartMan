#' bartClassDiag
#'
#' @description Displays a selection of diagnostic plots for a BART model.
#'
#' @param model a model created from either the BART, dbarts, or bartMachine package.
#' @param response The name of the response for the fit.
#'
#' @return A selection of diagnostic plots
#'
#' @import ggplot2
#' @importFrom patchwork area
#' @importFrom patchwork plot_layout
#' @importFrom ROCR prediction
#' @importFrom ROCR performance
#' @export

bartClassDiag <- function(model, response){

  responseVals <- response


  if(class(model) == "bartMachine"){
    yhatTrain <- model$y_hat_train
  }else{
    yhatTrain <- colMeans(model$yhat.train)
    yhatTrain <- pnorm(yhatTrain)
  }

  # get prediction using ROCR package:
  pred <- prediction(yhatTrain, responseVals)

  # get auc value
  auc <- performance(pred, "auc")
  auc <- auc@y.values[[1]]
  print(paste0("AUC: ", round(auc, 5)))

  # get false/true positive rates
  perfTF <- performance(pred, "tpr", "fpr")

  # create dataframe for ROC plot
  dfROC <- data.frame(fpr = perfTF@x.values[[1]],
                      tpr = perfTF@y.values[[1]])

  # create dataframe for fitted vals plot
  dfFitClassBart <- data.frame(fitted = yhatTrain,
                               actual = responseVals)

  threshold <- mean(dfFitClassBart$fitted)
  class <- as.numeric(yhatTrain > threshold)
  dfFitClassBart$class <- class


  # create data frame for histogram
  dfPnorm <- data.frame(vals = yhatTrain)


  # -------------------------------------------------------------------------

  ROC <- bartROC(dfROC)
  classF <- classFit(dfFitClassBart)
  histogram <- classHist(dfPnorm)
  vimp <- bartVimpClass(model)

  design <- c(
    area(1, 1, 3, 3),
    area(1, 5, 3, 7),
    area(5, 1, 7, 3),
    area(5, 5, 7, 7)
  )

  diagPlot <- ROC + classF + histogram + vimp + plot_layout(design = design)

  return(diagPlot)


}
# ROC curve for dbarts ----------------------------------------------------


bartROC <- function(data){

  p <- ggplot(data, aes(x = fpr, y = tpr)) +
    geom_line() +
    xlab('False positive rate') +
    ylab('True positive rate') +
    ggtitle('ROC') +
    geom_abline(intercept = 0, slope = 1, col = 'blue')+
    theme_bw()

  return(p)
}




# Class fitted ------------------------------------------------------------


classFit <- function(data){

  threshold <- mean(data$fitted)

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


classHist <- function(data){

  p <- ggplot(data, aes(vals)) +
    #geom_histogram(stat = 'bin', binwidth = 0.05) +
    geom_histogram(bins = 50, color = "blue", fill = "white") +
    ylab('') +
    xlab('Predicted probability') +
    ggtitle("Histogram") +
    theme_bw()

  return(p)

}



# VIMP --------------------------------------------------------------------

bartVimpClass <- function(model) {

  if (class(model) == "pbart") {
    # get variable importance
    vImp <- model$varcount.mean
    vImp <- dplyr::tibble(Variable = names(vImp), Importance = vImp)
  } else if(class(model) == "bartMachine"){
    bmVimp <- bartMachine::investigate_var_importance(model, num_replicates_for_avg = 5)
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
