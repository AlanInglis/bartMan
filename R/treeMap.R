#' treeMap
#'
#' @description Creates a treemap displaying the frequency of different tree structures.
#'
#' @param treeList A list of trees created using the treeList function.
#'
#' @return A treemap plot.
#'
#'
#' @import ggplot2
#' @importFrom purrr map_chr
#' @importFrom purrr reduce
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


treeMap <- function(treeList){

  # get edge frequency
  edgeFreq <- treeList %>%
    map_chr(~as_tibble(activate(.x, edges)) %>%
              map_chr(str_c, collapse = " ") %>%
              toString())%>%
    table() %>%
    as_tibble() %>%
    setNames(c("data", "frequency")) %>%
    separate(data, c("from", "to"), ", ") %>%
    filter(. != "") %>%
    arrange(-frequency) %>%
    mutate(treeNum = row_number())

  # create treemap of tree structure frequency
  tMap <- ggplot(edgeFreq,
         aes(fill = frequency,
             area = frequency,
             label = treeNum
             )) +
    theme(legend.position = "none") +
    treemapify::geom_treemap() +
   treemapify::geom_treemap_text(colour = "white",
                     place = "centre")


# Create Legend -----------------------------------------------------------

  # extract all edges from tree list
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


# Create areas for legend -------------------------------------------------

  vals = seq(1:length(allPlots))
  maxVal <- max(vals)

  # create list of all areas
  areaList <- lapply(vals, function(x) area(x, maxVal+1))
  # turn into a df
  areaReduced <- purrr::reduce(areaList, c)
  # put all areas together
  design <- c(area(1, 1, maxVal, maxVal),
              areaReduced)


# Plot treeMap ------------------------------------------------------------

 tMapPlot <- tMap + allPlots + plot_layout(design = design)


  return(tMapPlot)



}
