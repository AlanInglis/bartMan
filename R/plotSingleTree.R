#' plotSingleTree
#'
#' @description Plots individual trees.
#'
#' @param treeData A data frame created by \code{extractTreeData} function
#' @param iter The MCMC iteration or chain to plot.
#' @param treeNo The tree number to plot.
#' @param plotType What type of plot to display. either dendrogram or icicle.
#' @return A plot of an individual tree
#'
#' @import ggraph
#' @import ggplot2
#'
#' @export

plotSingleTree <- function(treeData,
                           iter = 1,
                           treeNo = 1,
                           plotType =  "icicle") {

  # Error checks for 'iter' and 'treeNo' to ensure they are not NULL
  if(is.null(iter)) {
    stop("Error: 'iter' cannot be NULL")
  }
  if(is.null(treeNo)) {
    stop("Error: 'treeNo' cannot be NULL")
  }

  p <- treeList(treeData, iter = iter, treeNo = treeNo)

  if (plotType == "dendrogram") {
    gp <- ggraph(p[[1]], "dendrogram") +
      geom_edge_elbow() +
      geom_node_label(aes(label = label, color = label)) +
      theme_graph() +
      theme(legend.position = "none")
  } else if (plotType == "icicle") {
    gp <- ggraph(p[[1]], "partition") +
      geom_node_tile(aes(fill = var), linewidth = 0.25) +
      geom_node_label(aes(label = label, color = var)) +
      scale_y_reverse() +
      theme_void() +
      theme(legend.position = "none")
  }

  return(gp)
}



# -------------------------------------------------------------------------
