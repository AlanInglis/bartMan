#' splitDensity
#'
#' @description Density plots of the split value for each variable.
#'
#' @param model a model created from either the BART dbarts, or bartMachine packages
#' @param treeList A list of trees created using the treeList function.
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


splitDensity <- function(model, treeList, colBy = NULL, display = "histogram") {

  # split into dataframe of trees
  dfTrees <- NULL
  for (i in 1:(length(treeList))) {
    dfTrees[[i]] <- treeList[[i]] %>%
      activate(nodes) %>%
      data.frame() %>%
      select(var, value) %>%
      na.omit()
  }

  # turn into one dataframe
  tt <- bind_rows(dfTrees, .id = "column_label")

  # create plotting order to match order of data
  nam <- colnames(model$varcount)
  tt$var <- factor(tt$var, levels = nam)

  # create plot

  if (display == "density") {
    dPlot <- tt %>%
      ggplot(aes(x = value)) +
      geom_density(aes(colour = var, fill = var)) +
      facet_wrap(~var) +
      theme_bw() +
      ylab("Density") +
      xlab("Split value") +
      theme(legend.position = "none")
  } else if (display == "ridge") {
    dPlot <- tt %>%
      ggplot(aes(x = value, y = var, fill = stat(x))) +
      geom_density_ridges(aes(fill = var, alpha = 0.1)) +
      theme_bw() +
      theme(legend.position = "none")
  } else if(display == "histogram") {
    dPlot <- tt %>%
      ggplot(aes(x = value)) +
      geom_histogram(aes(colour = var, fill = var), bins = 30) +
      facet_wrap(~var) +
      theme_bw() +
      ylab("Density") +
      xlab("Split value") +
      theme(legend.position = "none")
  }

  return(dPlot)
}
