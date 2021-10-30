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
#' @importFrom igraph gsize
#'
#' @export


plotAllTrees <- function(treeList, sampleSize = 0) {

  # plot a sample of trees
  if (sampleSize > 0) {
    treeList <- sample(treeList, sampleSize, replace = FALSE)
  }

 # remove stumps
 treeList <- Filter(function(x) igraph::gsize(x) > 0, treeList)

  # set node colours
 nodenames <- unique(na.omit(unlist(lapply(treeList, .%>%activate(nodes) %>% pull(var) ))))
 nodecolors <- setNames(scales::hue_pal(c(0,360)+15, 100, 64, 0, 1)(length(nodenames)), nodenames)


  allPlots <- lapply(treeList, plotFun, n = length(treeList), color = nodecolors)
  n <- length(allPlots)
  nRow <- floor(sqrt(n))
  do.call("grid.arrange", c(allPlots, nrow = nRow))
}

plotFun <- function(List, colors = NULL, n) {
  if (n > 10) {
    plot <- ggraph(List, "partition") +
      geom_node_tile(aes(fill = var), size = 0.25) +
      scale_y_reverse() +
      theme_void() +
      theme(legend.position = "none")
    if (!is.null(colors)) {
      plot <- plot + scale_fill_manual(values = colors) +
        scale_color_manual(values = colors, na.value = "grey")
    }
  } else {
    plot <- ggraph(List, "partition") +
      geom_node_tile(aes(fill = var), size = 0.25) +
      geom_node_label(aes(label = label, color = label)) +
      scale_y_reverse() +
      theme_void() +
      theme(legend.position = "none")
    if (!is.null(colors)) {
      plot <- plot + scale_fill_manual(values = colors) +
        scale_color_manual(values = colors, na.value = "grey")
    }
  }
  return(plot)
}
