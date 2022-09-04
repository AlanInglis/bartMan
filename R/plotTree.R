#' plotTree
#'
#' @description plots individual trees from the R-packages BART, dbarts, and
#' bartMachine
#'
#' @param treeData A data frame created by treeData function
#' @param treeNo The tree number to plot.
#' @param iter The MCMC iteration or chain to plot.
#' @param plotType What type of plot to display. either dendrogram or icicle.
#' @return A plot of an individual tree
#'
#' @import ggraph
#' @import ggplot2
#' @importFrom tidygraph activate
#' @importFrom tidygraph tbl_graph
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr one_of
#' @importFrom dplyr group_split
#' @importFrom dplyr transmute
#' @importFrom dplyr row_number
#' @importFrom dplyr pull
#' @importFrom dplyr filter
#' @importFrom dplyr n
#' @importFrom stats setNames
#' @importFrom tidygraph bfs_dist
#'
#' @export

plotTree <- function(treeData,
                     iter = 1,
                     treeNo = 1,
                     plotType =  "icicle") {

  p <- plotAll(treeData, iter = 1, treeNo = 1)

  if (plotType == "dendrogram") {
    gp <- ggraph(p[[1]], "dendrogram") +
      geom_edge_elbow() +
      geom_node_label(aes(label = label, color = label)) +
      theme_graph() +
      theme(legend.position = "none")
  } else if (plotType == "icicle") {
    gp <- ggraph(p[[1]], "partition") +
      geom_node_tile(aes(fill = var), size = 0.25) +
      geom_node_label(aes(label = label, color = var)) +
      scale_y_reverse() +
      theme_void() +
      theme(legend.position = "none")
  }

  return(gp)
}



# -------------------------------------------------------------------------
