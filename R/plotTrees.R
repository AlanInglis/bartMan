#' plotTrees
#'
#' @description plots individual trees from the R-packages BART, dbarts, and
#' bartMachine
#'
#' @param treeData A data frame created by treeList function
#' @param treeNum The tree number to plot.
#' @param iteration The MCMC iteration or chain to plot.
#' @param plotType What type of plot to display. either dendrogram or icicle.
#' @return A plot of an individual tree
#'
#' @import ggraph
#' @import ggplot2
#' @importFrom tidygraph activate
#' @importFrom dplyr %>%
#' @importFrom tibble as_tibble
#' @export

plotTrees <- function(treeData,
                      treeNum = 1,
                      iteration = 1,
                      plotType = c("dendrogram", "icicle")) {

  getTreeListNumber <- function(treeData, iter, tNum){
    which(sapply(treeData, function(x) {
      nodes <- tidygraph::activate(x, nodes) %>% tibble::as_tibble()
      all(nodes$iteration == iteration & nodes$treeNum == treeNum)
    }))
  }


  listIndex <- getTreeListNumber(treeData, iter = iteration, tNum = treeNum)

  if (plotType == "dendrogram") {
    gp <- ggraph(treeData[[listIndex]], "dendrogram") +
      geom_edge_elbow() +
      geom_node_label(aes(label = label, color = label)) +
      theme_graph() +
      theme(legend.position = "none")
  } else if (plotType == "icicle") {
    gp <- ggraph(treeData[[listIndex]], "partition") +
      geom_node_tile(aes(fill = var), size = 0.25) +
      geom_node_label(aes(label = label, color = var)) +
      scale_y_reverse() +
      theme(legend.position = "none")
  }

  return(gp)
}
