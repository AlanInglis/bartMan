#' vimpPlot
#'
#' @description Plot the variable importance for a BART model with the 25% and 75%
#' quantile.
#'
#' @param trees A data frame created by `extractTreeData` function.
#' @param type What value to return. Either the raw count 'count'
#' or the proportions 'prop' averaged over iterations.
#' @param plotType Which type of plot to return. Either a barplot 'barplot' with the
#' quantiles shown as a line, a point plot with the quantiles shown as a gradient 'point', or a
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
#'
#' @examples
#' if(requireNamespace("dbarts", quietly = TRUE)){
#'  # Load the dbarts package to access the bart function
#'  library(dbarts)
#'  # Get Data
#'  df <- na.omit(airquality)
#'  # Create Simple dbarts Model For Regression:
#'  set.seed(1701)
#'  dbartModel <- bart(df[2:6], df[, 1], ntree = 5, keeptrees = TRUE, nskip = 10, ndpost = 10)
#'
#'  # Tree Data
#'  trees_data <- extractTreeData(model = dbartModel, data = df)
#'  vimpPlot(trees = trees_data, plotType = 'point')
#' }
#'
#' @export
#'

vimpPlot <- function(trees,
                     type = "prop",
                     plotType = "barplot",
                     metric = "median") {

  # warning
  if (!(type %in% c("val", "prop"))) {
    stop("type must be \"val\", or \"prop\"")
  }


  if (!(plotType %in% c("barplot", "point", "lvp"))) {
    stop("type must be \"barplot\", \"point\"  or \"lvp\"")
  }

  if (!(metric %in% c("mean", "median"))) {
    stop("metric must be \"mean\"  or \"median\"")
  }

  vimp <- vimpBart(trees, type = type)
  vimp <- as.data.frame(vimp)


  # get quantiles of proportions
  vimp25 <- apply(vimp, 2, function(x) quantile(x, c(.25)))
  vimp75 <- apply(vimp, 2, function(x) quantile(x, c(.75)))

  # get median and mean value
  vimpMed <- apply(vimp, 2, function(x) quantile(x, c(.5)))
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
      xlab("") +
      ylab("Importance") +
      theme(
        axis.title.y = element_text(angle = 90, vjust = 0.5),
        legend.key.size = unit(0.5, "cm")
      )
  } else if (plotType == "point") {

    if (!requireNamespace("ggforce", quietly = TRUE)) {
      stop("Package \"ggforce\" needed for this function to work. Please install it.",
           call. = FALSE)
    }

    p <- vImp %>%
      arrange(imp) %>%
      mutate(Variable = factor(Variable, unique(Variable))) %>%
      ggplot(aes(x = Variable, y = imp)) +
      ggforce::geom_link(aes(
        x = Variable, xend = Variable, yend = upperQ,
        colour = "gray50", alpha = rev(after_stat(index))
      ),
      size = 5, n = 1000
      ) +
      ggforce::geom_link(aes(
        x = Variable, xend = Variable, yend = lowerQ,
        colour = 'gray50', alpha = rev(after_stat(index))
      ),
      size = 5, n = 1000
      ) +
      geom_point(aes(x = Variable, y = imp), shape = 18, size = 2, color = "black") +
      scale_colour_identity() +
      coord_flip() +
      theme_bw() +
      labs(x = "", y = "Importance") +
      theme(legend.position = "none")
  } else if (plotType == "lvp") {

    if (!requireNamespace("lvplot", quietly = TRUE)) {
      stop("Package \"lvplot\" needed for this function to work. Please install it.",
           call. = FALSE)
    }


    dfvimp <- utils::stack(as.data.frame(vimp))
    colnames(dfvimp) <- c('value', 'variable')

   # pal <- rev(colorRampPalette(RColorBrewer::brewer.pal(9,name = 'Blues'))(10))
   pal <- c("#08306B", "#084D96", "#1B69AF", "#3787C0", "#58A1CE",
            "#81BADA", "#ABCFE5", "#CBDEF0","#E0ECF7", "#F7FBFF")

     p <- ggplot(dfvimp, aes(stats::reorder(variable, value), value)) +
               # aes(x = variable, y = value)) +
    #  lvplot::geom_lv(aes(fill = ..LV..),
       lvplot::geom_lv(aes(fill = after_stat(LV)),
        conf = 0.5,
        outlier.colour = "blue",
        outlier.shape = 5,
        varwidth = TRUE,
        col = 'black'
      ) +
      # geom_jitter(width = 0.2, alpha = 0.08) +
       scale_fill_manual(values = pal) +
      #scale_fill_brewer(palette = "Blues", direction = -1) +
       labs(x = "", y = "Importance") +
      theme_bw() +
      theme(legend.position = "none") + coord_flip()
  }

  return(p)
}






