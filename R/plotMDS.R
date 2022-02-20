#' plotMDS
#'
#' @description Multi-dimensional Scaling Plot of proximity matrix from a BART model.
#'
#' @param matrix A matrix of proximities created by the proximityMatrix function
#' @param pal A vector of colours used to colour the facrtor levels..
#' @param factorLevel A factorized response used to train BART model.
#'
#' @return A MDS plot.
#'
#' @importFrom stats cmdscale
#' @importFrom ggplot2 ggplot
#' @export
#'

plotMDS <- function(matrix, pal = NULL, factorLevel){

  fit <- cmdscale(1 - matrix, eig = TRUE, k = 2)

  nlevs <- nlevels(factorLevel)

  suppressWarnings(
  if (is.null(pal)) {
    pal <- if (requireNamespace("RColorBrewer", quietly = TRUE) &&
                   nlevs < 12)
      RColorBrewer::brewer.pal(nlevs, "Set1")
    else rainbow(nlevs)
  }
  )

  pal <- pal[as.numeric(factorLevel)]

  mdsDF <- data.frame(x = fit$points[, 1],
                      y = fit$points[, 2]
  )


  p <- ggplot(mdsDF, aes(x, y)) +
    geom_point(size = 2, col = pal) +
    xlab("Dimension 1") +
    ylab("Dimension 2") +
    theme_bw()

  return(p)

}
