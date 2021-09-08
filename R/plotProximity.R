#' plotProximity
#'
#' @description Plot a proximity matrix
#'
#' @param matrix A matrix of proximities created by the proximityMatrix function
#' @param pal A vector of colours to show proximity scores, for use with scale_fill_gradientn.
#' @param limit Specifies the fit range for the color map for proximity scores.
#'
#' @return A plot of proximity values.
#'
#' @importFrom reshape melt
#' @importFrom ggplot2 ggplot
#' @export


plotProximity <- function(matrix,
                          pal = rev(colorspace::sequential_hcl(palette = "Blues 2", n = 100)),
                          limit = NULL){

  # Set the limits
  if (is.null(limit)) {
    limit <- range(as.dist(matrix))
    limit <- range(labeling::rpretty(limit[1], limit[2]))
  }

  # melt the matrox
  suppressWarnings(
    df <- reshape::melt(matrix)
  )

  # order axis names
  varNames <- unique(df$values)
  df$values   <- factor(df$values,   levels = varNames)
  df$values.1 <- factor(df$values.1, levels = varNames)

  p <- ggplot(df, aes(values, values.1)) +
    geom_point(aes(size = value), show.legend = F) +
    geom_tile(aes(fill = value)) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(limits = rev(levels(df$values.1))) +
    scale_fill_gradientn(
      colors = pal, limits = limit, name = "Proximity",
      guide = guide_colorbar(
        order = 1,
        frame.colour = "black",
        ticks.colour = "black"
      ), oob = scales::squish
    ) +
    xlab("") +
    ylab("") +
    theme_bw() +
    theme(aspect.ratio = 1)

  return(p)

}
