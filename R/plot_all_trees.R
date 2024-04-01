#' Plot Trees with Customizations
#'
#' This function plots trees from a list of tidygraph objects. It allows for various
#' customisations such as fill colour based on node response or value, node size adjustments,
#' and color palettes.
#'
#' @param trees A data frame of trees.
#' @param iter An integer specifying the iteration number of trees to be included in the output.
#'             If NULL, trees from all iterations are included.
#' @param treeNo An integer specifying the number of the tree to include in the output.
#'               If NULL, all trees are included.
#' @param fillBy A character string specifying the attribute to color nodes by.
#'               Options are 'response' for coloring nodes based on their mean response values or
#'               'mu' for coloring nodes based on their predicted value, or NULL for no
#'               specific fill attribute.
#' @param sizeNodes A logical value indicating whether to adjust node sizes.
#'                  If TRUE, node sizes are adjusted; if FALSE, all nodes are given the same size.
#' @param removeStump A logical value. If TRUE, then stumps are removed from plot.
#' @param selectedVars A vector of selected variables to display. Either a character vector of names
#'                  or the variables column number.
#' @param pal A colour palette for node colouring. Palette is used when 'fillBy' is specified for gradient colouring.
#' @param center_Mu A logical value indicating whether to center the color scale for the 'mu'
#'                  attribute around zero. Applicable only when 'fillBy' is set to "mu".
#' @param cluster A character string that specifies the criterion for reordering trees in the output.
#'                Currently supports "depth" for ordering by the maximum depth of nodes, and "var" for a
#'                clustering based on variables. If NULL, no reordering is performed.
#'
#' @return A ggplot object representing the plotted trees with the specified customizations.
#'
#' @importFrom purrr map
#' @importFrom tidygraph activate pull tbl_graph
#' @importFrom dplyr as_tibble mutate case_when
#' @importFrom tidyr replace_na
#' @importFrom scales hue_pal squish
#' @importFrom  ggnewscale new_scale_fill new_scale_color
#' @importFrom igraph gsize
#' @import ggplot2
#' @import ggraph
#'
#' @examples
#' plotFun(my_trees, fillBy = 'response', sizeNodes = TRUE)
#'
#' @export

plotTrees <- function(trees,
                      iter = NULL,
                      treeNo = NULL,
                      fillBy = NULL,
                      sizeNodes = FALSE,
                      removeStump = FALSE,
                      selectedVars = NULL,
                      pal = rev(colorRampPalette(c('steelblue', '#f7fcfd', 'orange'))(5)),
                      center_Mu = TRUE,
                      cluster = NULL){


  # get variable names
  variable_names <- trees$varName

  # make list of trees
  trees <- treeList(trees = trees,
                    iter = iter,
                    treeNo = treeNo)

  # sort graph
  if(is.null(cluster)){
    facet_name = "iteration"
  } else{
    if (cluster == "var") {
      trees <- clusterTrees(trees)
      facet_name <- "varID"
    } else if (cluster == "depth") {
      facet_name <- "-depthMax"
    }
  }

  # remove stump
  if (removeStump) {
    # Remove trees with no edges
    trees <- Filter(function(x) igraph::gsize(x) > 0, trees)
  }

  # get unique 'var' values from each tree
  all_vars <- purrr::map(trees, function(tree) {
    tree %>%
      activate(nodes) |>
      as_tibble() |>
      pull(var) |>
      unique()
  }) |>
    unlist() |>
    unique()



  # Select variables to display
  if(!is.null(selectedVars)){
    if(is.numeric(selectedVars)){
      selected_names <- variable_names[selectedVars]
      not_selected_names <- variable_names[-selectedVars]
    }else if(is.character(selectedVars)){
      selected_names <- variable_names[variable_names %in% selectedVars]
      not_selected_names <- variable_names[!(variable_names %in% selectedVars)]
    }
  }


  # set stump name
  if (!is.null(fillBy)) {
    stump_name <- "Stump"
  } else {
    stump_name <- "Stump/Leaf"
  }


  # Initialize empty data frames for nodes and edges
  all_nodes <- data.frame()
  all_edges <- data.frame()

  # Extract and combine nodes and edges from each tree
  for (i in seq_along(trees)) {
    nodes <- trees[[i]] |>  activate(nodes) |>  as_tibble()
    edges <- trees[[i]] |>  activate(edges) |>  as_tibble()

    # Adjust node IDs in edges to make them unique across combined graph
    edge_offset <- ifelse(nrow(all_nodes) == 0, 0, max(all_nodes$node))
    edges$from <- edges$from + edge_offset
    edges$to <- edges$to + edge_offset

    # Combine
    all_nodes <- rbind(all_nodes, mutate(nodes, node = node + edge_offset))
    all_nodes$var[is.na(all_nodes$var)] <- stump_name

    all_edges <- rbind(all_edges, edges)
  }


  # rename other variables if selected vars
  if(!is.null(selectedVars)){
    all_nodes$var <- ifelse(all_nodes$var %in% not_selected_names, "Others", all_nodes$var)
  }


  # add varID columns for use when clustering by var
  all_nodes <- all_nodes |>
    mutate(varID = cumsum(c(1, diff(treeNum) != 0)))

  # get the limits
  leaf_stumps <- all_nodes |>
    filter(var == "Stump")

  if(is.null(fillBy)) {
    fill_value <- NULL
    lims <- NULL
    legend_name <- 'Variable'
  } else if(fillBy == 'response') {
    lims <- range(leaf_stumps$respNode, na.rm = TRUE)
    lims <- pretty(c(lims[1], lims[2]))
    lims <- c(min(lims), max(lims))
    legend_name <- 'Mean \nResponse'
    all_nodes <- all_nodes |>
      mutate(fill_value = respNode)
  } else if(fillBy == "mu") {
    lims <- range(leaf_stumps$value, na.rm = TRUE)
    if(center_Mu){
      lims <- c(-max(abs(lims)), max(abs(lims)))
    }else{
      lims <- pretty(c(lims[1], lims[2]))
      lims <- c(min(lims), max(lims))
    }
    legend_name <- 'Mu'
    all_nodes <- all_nodes |>
      mutate(fill_value = value)
  }



  # set node colours
  nodeNames <-  unique(all_nodes$var) #sort(all_vars)
  nodecolors <- setNames(scales::hue_pal(c(0, 360) + 15, 100, 64, 0, 1)(length(nodeNames)), nodeNames)



  # set stump colour
  if(is.null(fillBy)){
    nodecolors[[stump_name]] <- '#808080'
  } else{
    mean_value <- ifelse(fillBy == 'response', all_nodes$respNode[1],
                         if(fillBy == 'mu') mean(leaf_stumps$value))
    # make sure stumps are coloured by mean value
    all_nodes$fill_value <- ifelse(all_nodes$isStump, mean_value, all_nodes$fill_value)

    max_colours <- 1000000  # threshold for the maximum number of colors
    # error handling if too many colours are selected
    if(lims[2] > max_colours) {
      warning(paste("lims[2] is too large (", lims[2], "). Using max_colors =", max_colours))
      pal_stump <- rev(colorRampPalette(c('steelblue', '#f7fcfd', 'orange'))(max_colours))
    } else {
      pal_stump <- rev(colorRampPalette(c('steelblue', '#f7fcfd', 'orange'))(length(nodecolors)))
    }

    nodecolors[[stump_name]] <- get_stump_colour_for_legend(lims = lims,
                                                            mean_value = mean_value,
                                                            palette = pal_stump)
  }


  # set "Others" node colour if selected
  if(!is.null(selectedVars)){
    nodecolors[['Others']] <- '#e6e6e6'
  }

  if(removeStump){
    nodecolors <- nodecolors[setdiff(names(nodecolors), stump_name)]
  }

  # set node size
  if(sizeNodes){
    all_nodes <- all_nodes |>
      mutate(dynamic_weight = noObs)
  }else{
    all_nodes <- all_nodes |>
      mutate(dynamic_weight = 1)
  }


  # Create the combined graph
  combined_graph <- tbl_graph(nodes = all_nodes, edges = all_edges, directed = TRUE)

  # reconfigure var column for plotting
  combined_graph <- combined_graph |>
    mutate(var = dplyr::case_when(
      var == 'Stump' & isStump == FALSE ~ NA_character_,
      TRUE ~ var
    ))



  # plot set up
  num_plots <- length(trees)

  # Dynamic settings based on the number of plots
  if (num_plots <= 250) {
    panel_spacing_x <- unit(1, "lines")
  } else if (num_plots <= 500){
    panel_spacing_x <- unit(0.5, "lines")
  } else {
    panel_spacing_x <- unit(0.25, "lines")
  }

  if(!is.null(cluster)){
    facet_formula <- as.formula(paste0("~ ", facet_name, " + iteration + treeNum"))
  } else{
    facet_formula <- as.formula(paste0("~ iteration + treeNum"))
  }

  suppressMessages(
  p <- ggraph(combined_graph, layout = "partition", weight = dynamic_weight) +
    geom_node_tile(aes(fill = var), linewidth = 0.25) +
    scale_y_reverse() +
    scale_fill_manual(values = nodecolors,
                      name = "Variable",
                      na.value = "#808080") +
    facet_nodes(facets = facet_formula, scales = "free") +
    theme_void() +
    theme(aspect.ratio = 1,
          legend.position = "right",
          panel.spacing.x = unit(panel_spacing_x, "lines"),
          strip.text.x = element_text(size = 0)) +
    ggnewscale::new_scale_fill() +
    ggnewscale::new_scale_color() +
    geom_node_tile(linewidth = 0.15,
                   data = . %>% filter(is.na(var)),
                   aes(fill = fill_value)) +
    scale_fill_gradientn(
      colours = pal,
      limits = lims,
      name = legend_name,
      guide = guide_colorbar(
        frame.colour = "black",
        ticks.colour = "black",
        order = 2
      ),
    )
  )

  return(p)
}
# END





