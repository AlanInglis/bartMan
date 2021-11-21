#' treeBarPlotOLD
#'
#' @description Creates a barplot displaying the frequency of different tree structures.
#'
#' @param treeList A list of trees created using the treeList function.
#'
#' @return A barplot plot.
#'
#'
#' @import ggplot2
#' @importFrom purrr map_chr
#' @importFrom purrr map2
#' @importFrom tidygraph activate
#' @importFrom tidygraph tbl_graph
#' @importFrom tidyr separate
#' @importFrom dplyr as_tibble
#' @importFrom dplyr arrange
#' @importFrom dplyr mutate
#' @importFrom igraph gsize
#' @importFrom ggraph ggraph
#' @importFrom ggraph geom_node_tile
#' @importFrom ggraph geom_node_text
#' @importFrom patchwork area
#' @importFrom patchwork plot_layout
#' @importFrom treemapify geom_treemap
#' @importFrom treemapify geom_treemap_text
#'
#'
#' @export

treeBarPlotOLD <- function(treeList){


# Get Edge Frequency ------------------------------------------------------

  edgeFreq <- bmTreesList %>%
    map_chr(~as_tibble(activate(.x, edges)) %>%
              map_chr(str_c, collapse = " ") %>%
              toString())%>%
    table() %>%
    as_tibble() %>%
    setNames(c("data", "Frequency")) %>%
    separate(data, c("from", "to"), ", ") %>%
    filter(. != "") %>%
    arrange(-Frequency) %>%
    mutate(treeNum = row_number())


# Create barplot of frequency ---------------------------------------------

  bp <- edgeFreq %>%
    ggplot() +
    geom_bar(aes(x = from, y = Frequency, fill = Frequency), stat = "identity") +
    ggtitle("Tree Frequency") +
    ylab("Frequency") +
    xlab("") +
    theme_bw() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    guides(fill = guide_colourbar(
      frame.colour = "black",
      ticks.colour = "black"))


# Create trees for x-axis -------------------------------------------------

  # extract all edges
  edgeList <- NULL
  for(i in 1:length(treeList)){
    edgeList[[i]] <- treeList[[i]] %>%
      activate(edges) %>%
      data.frame()
  }


  # get list of unique plots
  edgeListKeep <- edgeList[!duplicated(lapply(edgeList, function(x) x[,c("from","to")]))]

  # remove stumps
  edgeListKeep <- Filter(NROW, edgeListKeep)

  # turn into table graph objects
  edgeListTBL <- NULL
  for(i in 1:length(edgeListKeep)){
    edgeListTBL[i] <- list(tbl_graph(edges = edgeListKeep[[i]]))
  }

  # add plot name as number
  for(i in 1:(length(edgeListTBL))){
    edgeListTBL[[i]] <- edgeListTBL[[i]] %>%
      activate(nodes) %>%
      mutate(name = c(i, rep("", length.out = igraph::gsize(edgeListTBL[[i]]))))
  }


  # plotting function
  plotFun <- function(List, n) {

    plot <- ggraph(List, "partition") +
      geom_node_tile(aes(fill = "red"), size = 0.25) +
      geom_node_text(aes(label = name), size = 4) +
      scale_y_reverse() +
      theme_void() +
      theme(legend.position = "none")
  }

  allPlots <- lapply(edgeListTBL, plotFun, n = length(edgeListTBL))


# Create final barplot ----------------------------------------------------


  width <- .9 # set default width of bars
  p_axis <- ggplot(edgeFreq) +
    geom_blank(aes(x = from)) +
    purrr::map2(allPlots, seq_along(allPlots), ~ annotation_custom(ggplotGrob(.x), xmin = .y - width / 2, xmax = .y + width / 2)) +
    theme_void()


  bpFinal <- bp / p_axis + plot_layout(heights = c(4, 1))

  return(bpFinal)


}
















