#' splitDensity
#'
#' @description Density plots of the split value for each variable.
#'
#' @param treeData A list of trees created using the treeData function.
#' @param data Data frame containing variables from the model.
#' @param colBy A parameter used to control the colour of the density plots.
#' @param display Choose how to display the plot. Either histogram, facet wrap or ridges.
#'
#' @return A faceted group of density plots
#'
#' @importFrom tidygraph activate
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr bind_rows
#' @importFrom ggridges geom_density_ridges
#' @import ggplot2
#'
#' @export


splitDensity <- function(treeData, data, colBy = NULL, display = "histogram") {

  if (!(display %in% c("histogram", "ridge", "density", 'both'))) {
    stop("display must be \"histogram\", \"ridge\", \"density\", or \"both\"")
  }

  # get just the variable and split value
  tt <- treeData$structure %>%
    ungroup() %>%
    select(var, splitValue) %>%
    na.omit()

  # create plotting order to match order of data
  nam <- treeData$varName
  tt$var <- factor(tt$var, levels = nam)

  varNames <- unique(tt$var)

  # create plot

  if (display == "density") {
    dPlot <- tt %>%
      ggplot(aes(x = splitValue)) +
      geom_density(aes(colour = var, fill = var)) +
      facet_wrap(~var) +
      theme_bw() +
      ylab("Density") +
      xlab("Split value") +
      theme(legend.position = "none")
  } else if (display == "ridge") {
    dPlot <- tt %>%
      ggplot(aes(x = splitValue, y = var, fill = stat(x))) +
      geom_density_ridges(aes(fill = var, alpha = 0.1)) +
      ylab("Variable") +
      xlab("Split value") +
      theme_bw() +
      theme(legend.position = "none")
  } else if(display == "histogram") {
    dPlot <- tt %>%
      ggplot(aes(x = splitValue)) +
      geom_histogram(aes(colour = var, fill = var), bins = 30) +
      facet_wrap(~var) +
      theme_bw() +
      ylab("Density") +
      xlab("Split value") +
      theme(legend.position = "none")
  }else if(display == 'both'){

    dataIdx <- which((names(data) %in% varNames))
    dat <- data[, dataIdx]

    meltDat <- melt(dat)
    names(tt) <- c('variable', 'value')

    dataList <- list(meltDat, tt)
    names(dataList)  <- c('dat', 'sv')
    dfList <- plyr::ldply(dataList)

   dPlot <- ggplot(dfList) +
      geom_density(aes(x = value, fill = .id), alpha = 0.5) +
      facet_wrap(~variable) +
      scale_fill_discrete(name = "", labels = c("Data", "Split Value")) +
      ylab('Density') +
      xlab("Split value") +
      theme_bw()
  }

  suppressMessages(print(dPlot))
  #return(dPlot)
}
