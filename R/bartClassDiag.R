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
#' @importfrom ROCR prediction
#' @importFrom ROCR performance
#' @export


bartClassDiag <- function(model, response){

  responseVals <- response

  # get prediction using ROCR package:
  pred <- prediction(colMeans(pnorm(model$yhat.train)), responseVals)

  # get auc value
  auc <- performance(pred,"auc")@y.values[[1]]

  # get false/true positive rates
  perfTF <- performance(pred, "tpr", "fpr")

  # create dataframe for ROC plot
  dfROC <- data.frame(fpr=perfTF@x.values[[1]],
                      tpr=perfTF@y.values[[1]])

  # get sens and spec vals
  perfSS <- performance(pred, "sens", "spec")

  # extract values
  SS <- (perfSS@x.values[[1]] + perfSS@y.values[[1]] - 1)

  # turn into datframe
  dfSS <- data.frame(alpha = perfSS@alpha.values[[1]], trueSkill = SS)

  # get threshold value
  thresholdVal <- min(dfSS$alpha[which(dfSS$trueSkill == max(dfSS$trueSkill))])

  # create dataframe for fitted vals plot
  dfFitClassBart <- data.frame(fitted = pnorm(colMeans(model$yhat.train)),
                               class = as.numeric(pnorm(colMeans(model$yhat.train)) > thresholdVal),
                               actual = responseVals)

  # create data frame for histogram
  dfPnorm <- data.frame(pnorm = colMeans(pnorm(model$yhat.train)))


# -------------------------------------------------------------------------

  ROC <- bartROC(dfROC)
  classF <- classFit(dfFitClassBart, thresholdVal)
  histogram <- classHist(dfPnorm)
  thres <- thresPerf(dfSS, thresholdVal)

  design <- c(
    area(1, 1, 3, 3),
    area(1, 5, 3, 7),
    area(5, 1, 7, 3),
    area(5, 5, 7, 7)
  )

  diagPlot <- ROC + classF + histogram + thres + plot_layout(design = design)

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
  ylab('True classification') +
  ggtitle('Class of fitted values') +
  xlim(0,1) +
  theme_bw() +
  theme(legend.position = 'none') +
  geom_vline(xintercept = threshold, col = 'black')

return(p)

}

# Histogram fitted vals ---------------------------------------------------


classHist <- function(data){

  p <- ggplot(data, aes(pnorm)) +
    #geom_histogram(stat = 'bin', binwidth = 0.05) +
    geom_histogram(bins = 50, color = "blue", fill = "white") +
    ylab('Number of training data points') +
    xlab('Predicted probability') +
    ggtitle("Histogram") +
    theme_bw()

  return(p)

}



# Threshold performace ----------------------------------------------------

thresPerf <- function(data, threshold){

  p <- ggplot(data, aes(x = alpha, y = trueSkill)) +
    geom_line() +
    ggtitle('Threshold-performance curve') +
    xlab('Threshold') +
    ylab('True skill statistic') +
    geom_vline(xintercept = threshold, col = 'blue')+
    theme_bw()

  return(p)
}













