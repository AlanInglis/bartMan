#' treeBarPlot
#'
#' @description Creates a barplot displaying the frequency of different tree structures.
#'
#' @param treeList A list of trees created using the treeList function.
#' @param topTrees integer value to show the top x variables.
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
#' @importFrom tidyr replace_na
#' @importFrom dplyr as_tibble
#' @importFrom dplyr arrange
#' @importFrom dplyr mutate
#' @importFrom igraph gsize
#' @importFrom ggraph ggraph
#' @importFrom ggraph geom_node_tile
#' @importFrom ggraph geom_node_text
#' @importFrom patchwork area
#' @importFrom patchwork plot_layout
#' @importFrom cowplot get_legend
#' @importFrom cowplot plot_grid
#' @importFrom treemapify geom_treemap
#' @importFrom treemapify geom_treemap_text
#'
#'
#' @export

treeBarPlot <- function(treeList, topTrees = NULL){

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
  freqDf$treeNum <- seq(1:nrow(freqDf)) # add tree number


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


# Create barplot of frequency ---------------------------------------------
  names <- freqDf$value
  bp <- freqDf %>%
    ggplot() +
    geom_bar(aes(x = value, y = frequency, fill = frequency), stat = "identity") +
    scale_x_discrete(limits = names) +
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



# ggraph plotting funtion -------------------------------------------------


  # set node colours
  nodenames <- unique(na.omit(unlist(lapply(treeList, .%>%activate(nodes) %>% pull(var)))))
  nodenames <- sort(nodenames)
  nodecolors <- setNames(scales::hue_pal(c(0,360)+15, 100, 64, 0, 1)(length(nodenames)), nodenames)

  # plotting function
  plotFun <- function(List, colors = NULL, n) {

    plot <- ggraph(List, "partition") +
      geom_node_tile(aes(fill = var), size = 0.25) +
      geom_node_text(aes(label = ''), size = 4) +
      theme(legend.position = "bottom") +
      scale_y_reverse() +
      theme_void()
    if (!is.null(colors)) {
      plot <- plot + scale_fill_manual(values = colors, name = "Variable") +
        scale_color_manual(values = colors, na.value = "grey")  +
        theme(legend.position = "bottom")
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

# Create final barplot ----------------------------------------------------


  width <- 1 # set default width of bars
  p_axis <- ggplot(freqDf) +
    geom_blank(aes(x = value)) +
    purrr::map2(allPlots, seq_along(allPlots), ~ annotation_custom(ggplotGrob(.x),  .y - width / 2, xmax = .y + width / 2)) +
    theme_void()

  bpAxis <- bp / p_axis + plot_layout(heights = c(4, 1))
  bpFinal <- cowplot::plot_grid(bpAxis, legend, rel_heights = c(1, .1), ncol = 1)

  return(bpFinal)

}

