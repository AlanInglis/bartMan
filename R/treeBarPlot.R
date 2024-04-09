#' Plot Frequency of Tree Structures
#'
#' Generates a bar plot showing the frequency of different tree structures
#' represented in a list of tree graphs. Optionally, it can filter to show only the top N trees
#' and handle stump trees specially.
#'
#' @param trees A list of tree graphs to display
#' @param iter Optional; specifies the iteration to display.
#' @param topTrees Optional; the number of top tree structures to display. If NULL, displays all.
#' @param removeStump Logical; if TRUE, trees with no edges (stumps) are excluded from the display
#'
#' @details
#' This function processes a list of tree structures to compute the frequency of each unique structure,
#' represented by a bar plot. It has options to exclude stump trees (trees with no edges) and to limit
#' the plot to the top N most frequent structures.
#'
#' @return A `ggplot` object representing the bar plot of tree frequencies.
#'
#' @examples
#' if(requireNamespace("dbarts", quietly = TRUE)){
#'  # Load the dbarts package to access the bart function
#'  library(dbarts)
#'  # Get Data
#'  df <- na.omit(airquality)
#'  # Create Simple dbarts Model For Regression:
#'  set.seed(1701)
#'  dbartModel <- bart(df[2:6], df[, 1], ntree = 5, keeptrees = TRUE, nskip = 10, ndpost = 10)
#'
#'  # Tree Data
#'  trees_data <- extractTreeData(model = dbartModel, data = df)
#'  plot <- treeBarPlot(trees = trees_data, topTrees = 3, removeStump = TRUE)
#' }
#'
#' @importFrom dplyr mutate pull group_by arrange filter slice n
#' @importFrom purrr map imap
#' @importFrom tidyr replace_na
#' @import ggplot2
#' @import ggraph
#' @import patchwork
#' @importFrom cowplot plot_grid get_legend
#' @importFrom tidygraph activate bind_nodes bind_edges
#' @export

treeBarPlot <- function(trees, iter = NULL, topTrees = NULL, removeStump = FALSE) {
  # Create list of trees
  treeList <- treeList(trees = trees, iter = iter, treeNo = NULL)

  # Optionally remove stumps
  if (removeStump) {
    treeList <- Filter(function(x) igraph::gsize(x) > 0, treeList)
  }

  # Get frequencies of similar trees
  freqs <- purrr::map(treeList, function(x) {
    x |>
      dplyr::pull(var) |>
      tidyr::replace_na("..") |>
      paste0(collapse = "")
  }) |>
    unlist(use.names = F) |>
    dplyr::as_tibble() |>
    dplyr::mutate(ids = 1:dplyr::n()) |>
    dplyr::group_by(value) |>
    dplyr::mutate(val = dplyr::n():1)

  # Sort frequency data frame and rename columns
  freqDf <- freqs |>
    dplyr::slice(1) |>
    dplyr::arrange(-val) |>
    dplyr::rename(frequency = val)
  freqDf$treeNum <- seq(1:nrow(freqDf))

  # Limit to top trees if specified
  if (!is.null(topTrees)) {
    if (length(freqDf$ids) < topTrees) {
      stop("Number of trees chosen to display is greater than available trees. Alter the topTrees argument.")
    }
    freqDf <- freqDf[1:topTrees, ]
  }

  # Update frequencies and filter tree list
  ids <- freqDf$ids
  freqs <- freqs[ids, ] |> dplyr::pull(val)
  treeList <- treeList[ids]
  treeList <- purrr::imap(treeList, ~ .x |> dplyr::mutate(frequency = freqs[.y]) |> dplyr::select(var, frequency))

  # Generate barplot of tree frequencies
  names <- factor(freqDf$value, levels = freqDf$value)
  bp <- freqDf |>
    ggplot() +
    geom_bar(aes(x = value, y = frequency), fill = "steelblue", stat = "identity") +
    scale_x_discrete(limits = rev(levels(names))) +
    ylab("Count") +
    xlab("") +
    theme_bw() +
    theme(legend.position = "none") +
    coord_flip()

  # find if any stumps
  is_stump <- which(sapply(treeList, function(tree) igraph::gsize(tree) == 0))
  # Optional stump processing
  if(length(is_stump) >= 1){
    if (!removeStump) {
      tree_stumps <- treeList[is_stump]
      frequencies <- tree_stumps[[1]] |> tidygraph::activate(nodes) |> dplyr::pull(frequency)
      tree_stumps[[1]] <- tree_stumps[[1]] |> tidygraph::activate(nodes) |> dplyr::mutate(var = "Stump") |>
        tidygraph::bind_nodes(data.frame(var = "Stump", frequency = frequencies)) |>
        tidygraph::bind_edges(data.frame(from = c(1, 1), to = c(1, 1)))
      treeList[[is_stump]] <- tree_stumps[[1]]
      }
    }

  # Set node colors
  nodenames <- unique(stats::na.omit(unlist(lapply(treeList, function(tree) {
    tree |>
      tidygraph::activate(nodes) |>
      dplyr::pull(var)
  }))))
  nodenames <- sort(nodenames)
  nodecolors <- setNames(scales::hue_pal(c(0, 360) + 15, 100, 64, 0, 1)(length(nodenames)), nodenames)

  # Optional stump color setting
  if (!removeStump && length(is_stump) >= 1) {
    nodecolors[["Stump"]] <- '#808080'
  }

  # Define plot function
  plotFun <- function(List, colors = NULL, n) {
    ggraph(List, "partition") +
      geom_node_tile(aes(fill = var), linewidth = 0.25) +
      geom_node_text(aes(label = ""), size = 4) +
      theme(legend.position = "bottom") +
      scale_y_reverse() +
      theme_void() +
      scale_fill_manual(values = colors, name = "Variable", na.value = "#808080") +
      scale_color_manual(values = colors, na.value = "grey") +
      theme(aspect.ratio = 1)
  }

  # Generate plots for each tree
  allPlots <- lapply(treeList, plotFun, n = length(treeList), colors = nodecolors)

  # Create common legend
  ggdf <- data.frame(x = names(nodecolors), y = 1:length(nodecolors))
  ggLegend <- ggplot(ggdf, aes(x = x, y = y)) +
    geom_point(aes(color = x), shape = 15, size = 8) +
    scale_color_manual(values = unname(nodecolors), labels = names(nodecolors), name = "Variable") +
    theme_bw() +
    theme(legend.position = "bottom")
  legend <- cowplot::get_legend(ggLegend)

  # Remove legends from individual plots
  allPlots <- lapply(allPlots, function(x) x + theme(legend.position = "none"))

  # Remove vertical line from stump if not removed
  if(length(is_stump) >= 1){
    if (!removeStump) {
      allPlots[[is_stump]]$data <- allPlots[[is_stump]]$data[-2, ]
    }
  }

  # Filter top plots if specified
  if (!is.null(topTrees)) {
    allPlots <- allPlots[1:topTrees]
  }

  # Combine barplot and tree plots
  width <- 1
  p_axis <- ggplot(freqDf) +
    geom_blank(aes(y = value)) +
    purrr::map2(allPlots, rev(seq_along(allPlots)), ~ annotation_custom(ggplotGrob(.x), ymin = .y - width / 2, ymax = .y + width / 2)) +
    theme_void()

  bp1 <- bp + theme(axis.text.y = element_blank())
  ppp <- p_axis + theme(aspect.ratio = 10)
  px <- ppp | bp1
  bpFinal <- cowplot::plot_grid(px, legend, rel_heights = c(1, .1), ncol = 1)

  return(bpFinal)
}
