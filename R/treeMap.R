#' treeMap
#'
#' @description Creates a treemap displaying the frequency of different tree structures.
#'
#' @param treeList A list of trees created using the treeList function.
#' @param topTrees Integer of the number of plots to display. Starting from the most frequent and then decreasing.
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
#' @importFrom dplyr slice
#' @importFrom dplyr pull
#' @importFrom dplyr rename
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



treeMap <- function(treeList, topTrees = NULL){

  # remove stumps
  treeList <- Filter(function(x) igraph::gsize(x) > 0, treeList)

  # get the frequency off similar trees:
  freqs <- map(treeList, function(x){
    x %>%
      pull(var) %>%
      replace_na("..") %>%
      paste0(collapse = "")
  }) %>%
    unlist(use.names = F) %>%
    as_tibble() %>%
    mutate(ids = 1:n()) %>%
    group_by(value) %>%
    mutate(val = n():1)


  freqDf <-  freqs %>% slice(1) %>%  arrange(-val) %>%  rename(frequency = val)  # frequency tibble
  freqDf$treeNum <- seq(1:nrow(freqDf)) # a~dd tree number

  if(!is.null(topTrees)){
    freqDf <- freqDf[1:topTrees,]
  }

  ids <- freqs %>% slice(1) %>% pull(ids) # remove duplicates
  freqs <- freqs %>% pull(val) # get frequencies

  treeListNew <- purrr::imap(treeList, ~.x %>%
                           mutate(frequency = freqs[.y]) %>%
                           select(var, frequency))

  # return new list of trees
  treeList <- treeListNew[sort(ids)]

  # add plot name as number
  for(i in 1:(length(treeList))){
    treeList[[i]] <- treeList[[i]] %>%
      activate(nodes) %>%
      mutate(name = c(i, rep("", length.out = igraph::gsize(treeList[[i]]))))
  }

# Create treemap ----------------------------------------------------------


  # create treemap of tree structure frequency
  tMap <- ggplot(freqDf,
                 aes(fill = frequency,
                     area = frequency,
                     label = treeNum
                 )) +
    theme(legend.position = "none") +
    treemapify::geom_treemap() +
    treemapify::geom_treemap_text(colour = "white",
                                  place = "centre")



# ggraph plotting function ------------------------------------------------

  # set node colours
  nodenames <- unique(na.omit(unlist(lapply(treeList, .%>%activate(nodes) %>% pull(var)))))
  nodenames <- sort(nodenames)
  nodecolors <- setNames(scales::hue_pal(c(0,360)+15, 100, 64, 0, 1)(length(nodenames)), nodenames)

  # plotting function
  plotFun <- function(List, colors = NULL, n) {


    plot <- ggraph(List, "partition") +
      geom_node_tile(aes(fill = var), size = 0.25) +
     # geom_node_label(aes(label = var, color = var)) +
      geom_node_text(aes(label = name), size = 4) +
      scale_y_reverse() +
      theme_void()
    if (!is.null(colors)) {
      plot <- plot + scale_fill_manual(values = colors) +
        scale_color_manual(values = colors, na.value = "grey")
    }
  }

  allPlots <- lapply(treeList, plotFun, n = length(treeList), color = nodecolors)

  # get legend
  legend <- cowplot::get_legend(allPlots[[1]])

  # remove legends from individual plots
  allPlots <- lapply(allPlots, function(x) x + theme(legend.position = "none"))


# filter top X% of plots
  if(!is.null(topTrees)){
  allPlots <- allPlots[1:topTrees]
  }


# Create areas for legend -------------------------------------------------

  vals = seq(1:length(allPlots))
  maxVal <- max(vals)

  # create list of all areas
  areaList <- lapply(vals, function(x) area(x, maxVal+1))
  # turn into a df
  areaReduced <- purrr::reduce(areaList, c)
  # add legend
  legendArea <- area(floor(median(vals)), maxVal+2)
  # put all areas together
  design <- c(area(1, 1, maxVal, maxVal),
              areaReduced,
              legendArea)



# Plot treeMap ------------------------------------------------------------

  tMapPlot <- tMap + allPlots + legend +plot_layout(design = design)


  return(tMapPlot)





}
