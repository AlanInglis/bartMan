#' plotAllTrees
#'
#' @description Plots all the trees created by the treeList function.
#'
#' @param treeList A list of trees created by treeList function.
#' @param sampleSize Sample the tree list.
#'
#' @return A list containing vectors of the indices of observations from leaf nodes.
#'
#' @importFrom gridExtra grid.arrange
#' @importFrom ggraph ggraph
#' @importFrom dplyr %>%
#' @importFrom dplyr pull
#' @importFrom purrr keep
#' @importFrom tibble as_tibble
#'
#' @export


plotAllTrees <- function(treeList, sampleSize = 0) {

  # plot a sample of trees
  if (sampleSize > 0) {
    treeList <- sample(treeList, sampleSize, replace = FALSE)
  }

  # remove any of the lists that are only a single node
  treeList <- keep(treeList, ~ .x %>%
                   as_tibble %>%
                   pull(var) %>%
                   is.na %>%
                   `!` %>%
                   any)

  allPlots <- lapply(treeList, plotFun, n = length(treeList))
  n <- length(allPlots)
  nRow <- floor(sqrt(n))
  do.call("grid.arrange", c(allPlots, nrow = nRow))
}


plotFun <- function(t, n) {

  if (n > 10) {
    g <- ggraph(t, "partition") +
      geom_node_tile(aes(fill = var), size = 0.25) +
      scale_y_reverse() +
      theme_void() +
      theme(legend.position = "none")
  } else {
    g <- ggraph(t, "partition") +
      geom_node_tile(aes(fill = var), size = 0.25) +
      geom_node_label(aes(label = label, color = var)) +
      scale_y_reverse() +
      theme_void() +
      theme(legend.position = "none")
  }
  return(g)
}
