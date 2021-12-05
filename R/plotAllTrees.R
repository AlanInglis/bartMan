#' plotAllTrees
#'
#' @description Plots all the trees created by the treeList function.
#'
#' @param treeList A list of trees created by treeList function.
#' @param sampleSize Sample the tree list.
#' @param cluster LOGICAL. If TRUE, then cluster by tree structures.
#'
#' @return A list containing vectors of the indices of observations from leaf nodes.
#'
#' @importFrom gridExtra grid.arrange
#' @importFrom ggraph ggraph
#' @importFrom igraph gsize
#' @importFrom cowplot get_legend
#' @importFrom cowplot plot_grid
#'
#' @export


plotAllTrees <- function(treeList, sampleSize = 0, cluster = FALSE) {

  # plot a sample of trees
  if (sampleSize > 0) {
    treeList <- sample(treeList, sampleSize, replace = FALSE)
  }

 # remove stumps
 treeList <- Filter(function(x) igraph::gsize(x) > 0, treeList)

  # set node colours
 nodenames <- unique(na.omit(unlist(lapply(treeList, .%>%activate(nodes) %>% pull(var) ))))
 nodenames <- sort(nodenames)
 nodecolors <- setNames(scales::hue_pal(c(0,360)+15, 100, 64, 0, 1)(length(nodenames)), nodenames)

 # cluster by tree structure
 if(cluster){
 indIDS <- map(treeList, function(x){
   x %>%
     pull(var) %>%
     replace_na("a") %>%
     paste0(collapse = "")
 }) %>%
   unlist(use.names = F) %>%
   as_tibble() %>%
   mutate(ids = 1:n()) %>%
   group_by(value) %>%
   mutate(count = n():1) %>%
   arrange(value)

 ind <- indIDS  %>%
   group_by(value) %>%
   mutate(valrank = max(count)) %>%
   ungroup() %>%
   arrange(-valrank, value, -count) %>%
   pull(ids)


 treeList <- treeList[ind]
 }

  allPlots <- lapply(treeList, plotFun, n = length(treeList), color = nodecolors)
  # get legend
  legend <- cowplot::get_legend(allPlots[[1]])
  # remove legends from individual plots
  allPlots <- lapply(allPlots, function(x) x + theme(legend.position = "none"))
  n <- length(allPlots)
  nRow <- floor(sqrt(n))
  allTreesPlot <- do.call("grid.arrange", c(allPlots, nrow = nRow))

  cowplot::plot_grid(allTreesPlot, legend, rel_widths = c(1, .1), ncol = 2)

}


plotFun <- function(List, colors = NULL, n) {

  plot <- ggraph(List, "partition") +
    geom_node_tile(aes(fill = var), size = 0.25) +
    geom_node_text(aes(label = ""), size = 4) +
    scale_y_reverse() +
    theme_void()
  if (!is.null(colors)) {
    plot <- plot + scale_fill_manual(values = colors, name = "Variable") +
      scale_color_manual(values = colors, na.value = "grey")
  }
  return(plot)
}


