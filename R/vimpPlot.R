#' vimpPlot
#'
#' @description Plot the variable importance for a BART model with the 25% and 75%
#' quantile.
#'
#' @param treeData A data frame created by treeData function.
#' @param type What value to return. Either the raw count 'count'
#' or the proportions 'prop' averaged over iterations.
#' @param plotType Which type of plot to return. Either a barplot 'barplot' with the
#' quantiles shown as a line, a point plot with the quantiles shown as a gradient 'pointGrad', or a
#' letter-value plot 'lvp'.
#' @param metric Whether to show the 'mean' or 'median' importance values. Note, this has
#' no effect when using plotType = 'lvp'.
#'
#' @return A plot of variable importance.
#'
#' @import ggplot2
#' @importFrom dplyr tibble
#' @importFrom dplyr %>%
#' @importFrom dplyr arrange
#' @importFrom dplyr mutate
#' @importFrom ggforce geom_link
#' @importFrom reshape melt
#' @importFrom lvplot geom_lv
#'
#'
#' @export
#'

vimpPlot <- function(treeData, type = "prop", plotType = "barplot", metric = "median") {

  # warning
  if (!(type %in% c("val", "prop"))) {
    stop("type must be \"val\", or \"prop\"")
  }


  if (!(plotType %in% c("barplot", "pointGrad", "lvp"))) {
    stop("type must be \"barplot\", \"pointGrad\"  or \"lvp\"")
  }

  vimp <- vimpBart(treeData, type = type)
  vimp <- as.data.frame(vimp)

  # get quantiles of proportions
  vimp25 <- apply(vimp, 2, function(x) quantile(x, c(.25)))
  vimp75 <- apply(vimp, 2, function(x) quantile(x, c(.75)))

  # get median and mean value
  vimpMed <- apply(vimp, 2, function(x) median(x))
  vimpMean <- colMeans(vimp)

  if (metric == "mean") {
    importance <- vimpMean
  } else if (metric == "median") {
    importance <- vimpMed
  }

  # create df of values
  vImp <- dplyr::tibble(
    Variable = names(vimpMean),
    imp = importance,
    upperQ = vimp75,
    lowerQ = vimp25
  )

  if (plotType == "barplot") {
    p <- vImp %>%
      arrange(imp) %>%
      mutate(Variable = factor(Variable, unique(Variable))) %>%
      ggplot() +
      aes(x = Variable, y = vimpMedian) +
      geom_bar(aes(x = Variable, y = imp), stat = "identity", fill = "steelblue", col = "black") +
      geom_segment(aes(x = Variable, xend = Variable, y = lowerQ, yend = upperQ), color = "black") +
      theme_light() +
      coord_flip() +
      theme_bw() +
      xlab("Variable") +
      ylab("Importance") +
      theme(
        axis.title.y = element_text(angle = 90, vjust = 0.5),
        legend.key.size = unit(0.5, "cm")
      )
  } else if (plotType == "pointGrad") {
    p <- vImp %>%
      arrange(imp) %>%
      mutate(Variable = factor(Variable, unique(Variable))) %>%
      ggplot(aes(x = Variable, y = imp)) +
      ggforce::geom_link(aes(
        x = Variable, xend = Variable, yend = upperQ,
        col = Variable, alpha = rev(stat(index))
      ),
      size = 5, n = 1000
      ) +
      ggforce::geom_link(aes(
        x = Variable, xend = Variable, yend = lowerQ,
        col = Variable, alpha = rev(stat(index))
      ),
      size = 5, n = 1000
      ) +
      geom_point(aes(x = Variable, y = imp), shape = 18, size = 2, color = "black") +
      coord_flip() +
      theme_bw() +
      labs(x = "Variable", y = "Importance") +
      theme(legend.position = "none")
  } else if (plotType == "lvp") {
    suppressMessages(
      dfvimp <- reshape::melt(vimp)
    )

    p <- ggplot(dfvimp, aes(reorder(variable, value), value)) +
               # aes(x = variable, y = value)) +
      lvplot::geom_lv(aes(fill = ..LV..),
        conf = 0.5,
        outlier.colour = "blue",
        outlier.shape = 5,
        varwidth = TRUE
      ) +
      # geom_jitter(width = 0.2, alpha = 0.08) +
      scale_fill_brewer(palette = "Blues", direction = -1) +
      theme_bw() +
      theme(legend.position = "none")
  }

  return(p)
}






