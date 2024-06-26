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
#' @import ggplot2
#' @examples
#' if(requireNamespace("dbarts", quietly = TRUE)){
#'  # Load the dbarts package to access the bart function
#'  library(dbarts)
#'  # Get Data
#'  df <- na.omit(airquality)
#'  # Create Simple dbarts Model For Regression:
#'  set.seed(1701)
#'  dbartModel <- bart(df[2:6],
#'    df[, 1],
#'    ntree = 5,
#'    keeptrees = TRUE,
#'    nskip = 10,
#'    ndpost = 10)
#'  # Tree Data
#'  trees_data <- extractTreeData(model = dbartModel, data = df)
#'  # Create Proximity Matrix
#'  mProx <- proximityMatrix(trees = trees_data, reorder = TRUE, normalize = TRUE, iter = 1)
#'  # Plot
#'  plotProximity(matrix = mProx)
#' }
#' @export


plotProximity <- function(matrix,
                          pal = rev(colorspace::sequential_hcl(palette = "Blues 2", n = 100)),
                          limit = NULL){

  # Set the limits
  if (is.null(limit)) {
    limit <- range(stats::as.dist(matrix))
    limit <- range(pretty(c(limit[1], limit[2])))
      #range(labeling::rpretty(limit[1], limit[2]))
  }

  # melt the matrix
  # suppressWarnings(
  #   df <- reshape::melt(matrix)
  # )
  df <- utils::stack(as.data.frame(matrix))
  colnames(df) <- c('value', 'values.1')
  df$values.1 <- as.integer(as.character(df$values.1))
  df$values <- colnames(matrix)
  df$values <- as.integer(df$values)

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
