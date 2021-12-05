#' splitDensity
#'
#' @description Density plots of the split value for each variable.
#'
#' @param model a model created from either the BART dbarts, or bartMachine packages
#' @param treeList A list of trees created using the treeList function.
#' @param colBy A parameter used to control the colour of the density plots.
#' @param display Choose how to display the plot. Either facet wrap of overlayed.
#'
#' @return A faceted group of density plots
#'
#' @importFrom tidygraph activate
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr bind_rows
#' @import ggplot2
#'
#' @export


splitDensity <- function(model, treeList, colBy = NULL, display = "ridge") {

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

  intPal <- rev(colorspace::sequential_hcl(palette = "Purples 3", n = 100))
  # create plot

  if (display == "facet") {
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
  }

  return(dPlot)
}
