#' plotTrees
#'
#' @description plots individual trees from the R-packages BART, dbarts, and
#' bartMachine
#'
#' @param treeData A data frame created by XXX function
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

  list_obj <- lapply(treeData, function(x){
    edges <- tidygraph::activate(x, edges) %>% tibble::as_tibble()
    nodes <- tidygraph::activate(x, nodes) %>% tibble::as_tibble()
    return(list(edges = edges, nodes = nodes))
  } )


  getTreeListNumber <- function(treeData, iter, tNum){

    res <- 0
    listNumber <- NA

    for(i in 1:length(treeData)){
      res <- iter %in% treeData[[i]]$nodes$iteration && tNum %in% treeData[[i]]$nodes$treeNum
      if(res == TRUE){
        listNumber <- i
      }
    }
    if(is.na(listNumber)){
      stop("ERROR: Either iteration or treeNum subscript out of bounds")
    }

    return(listNumber)
  }

  listIndex <- getTreeListNumber(list_obj, iter = iteration, tNum = treeNum)

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
