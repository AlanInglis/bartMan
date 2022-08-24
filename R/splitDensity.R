#' splitDensity
#'
#' @description Density plots of the split value for each variable.
#'
#' @param treeData A list of trees created using the treeData function.
#' @param data Data frame containing variables from the model.
#' @param colBy A parameter used to control the colour of the density plots.
#' @param display Choose how to display the plot. Either histogram, facet wrap, ridges
#' or display both the split value and density of the predictor by using either both1 or both2.
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

  if (!(display %in% c("histogram", "ridge", "density", 'both1', 'both2'))) {
    stop("display must be \"histogram\", \"ridge\", \"density\", \"both1\", or \"both2\"")
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
  }else if(display == 'both1'){

    dataIdx <- which((names(data) %in% varNames))
    dat <- data[, dataIdx]

    meltDat <- melt(dat)
    names(tt) <- c('variable', 'value')

    dataList <- list(meltDat, tt)
    names(dataList)  <- c('data', 'split  \nvalue')
    dfList <- plyr::ldply(dataList)

   dPlot <- ggplot(dfList) +
      geom_density(aes(x = value, fill = .id), alpha = 0.5) +
      facet_wrap(~variable) +
      scale_fill_discrete(name = "", labels = c("Data", "Split Value")) +
      ylab('Density') +
      xlab("Split value") +
      theme_bw()
  }else if(display == 'both2'){

    dataIdx <- which((names(data) %in% varNames))
    dat <- data[, dataIdx]

    meltDat <- reshape2::melt(dat)
    names(tt) <- c('variable', 'value')

    dataList <- list(meltDat, tt)
    names(dataList)  <- c('data', 'split  \nvalue')
    dfList <- plyr::ldply(dataList)

    dPlot <- dfList %>%
      ggplot(aes(x = value, y = .id, fill = stat(x))) +
      geom_density_ridges(aes(fill = .id, alpha = 0.1), panel_scaling = F) +
      facet_wrap(~variable)+
      ylab("") +
      xlab("Value") +
      theme_bw() +
      theme(legend.position = "none")
  }

  suppressMessages(print(dPlot))
  #return(dPlot)
}
