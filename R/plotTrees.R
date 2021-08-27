#' plotTrees
#'
#' @description plots individual trees from the R-packages BART, dbarts, and
#' bartMachine
#'
#' @param model A bart model created from he R-packages BART, dbarts, or
#' bartMachine.
#' @param data Data frame used for fit.
#' @param response The name of the response for the fit.
#' @param treeNum The tree number to plot.
#' @param iter The MCMC iteration or chain to plot.
#' @param plotType What type of plot to display. either dendrogram or icicle.
#' @return A plot of an individual tree
#'
#' @importFrom ggraph ggraph
#' @export


# Main Function -----------------------------------------------------------


plotTrees <- function(model,
                      data,
                      response,
                      treeNum = 1,
                      iteration = 1,
                      plotType = c("dendrogram", "icicle")
                      ) {

  treeStructure <-  plotPackageTrees(model = model)



  if(plotType == "dendrogram"){
  gp <- ggraph(treeStructure[[treeNum]], "dendrogram", ) +
    geom_edge_elbow() +
    geom_node_label(aes(label = label, color = label)) +
    theme_graph() +
    theme(legend.position = "none")
  }else if(plotType == "icicle"){
   gp <-  ggraph(treeStructure[[treeNum]], "partition") +
    geom_node_tile(aes(fill = var), size = 0.25) +
    geom_node_label(aes(label = label, color = var)) +
    scale_y_reverse() +
    theme(legend.position = "none")
  }

  return(gp)
}


# Individual package methods: -------------------------------------------

# Main method:
plotPackageTrees <- function(model,
                             ...) {
  UseMethod("plotPackageTrees")
}


# BART package ------------------------------------------------------------

plotPackageTrees.wbart <- function(model) {
  bartTrees(model)
}

bartTrees <- function(model) {


  # variable names:
  varNames <- names(model$varcount.mean)

  # extracting tree structure
  trees <- list()
  trees$structure <- suppressWarnings(
    readr::read_table(
      file = model$treedraws$trees,
      col_names = c("node", "var", "cut", "leafValue"),
      col_types =
        readr::cols(
          node = readr::col_integer(),
          var = readr::col_integer(),
          cut = readr::col_integer(),
          leaf = readr::col_double()
        ),
      skip = 1,
      na = c("")
    )
  )

  trees$structure <- dplyr::mutate(
    trees$structure,
    tier = as.integer(floor(log2(node))),
    cut_id = cut + 1L,
    var = varNames[var + 1L]
  )

  # getting cut points
  cutPoints <- purrr::map_df(
    .x = model$treedraws$cutpoints,
    .f = ~ dplyr::tibble(cut = ., cut_id = 1:length(.)),
    .id = "var"
  )

  trees$structure <- dplyr::left_join(
    dplyr::select(trees$structure, -cut),
    cutPoints,
    by = c("var", "cut_id")
  )

  # define tree id and mcmc iteration number
  # first line contains mcmc draws
  firstLine <- strsplit(
    readr::read_lines(
      file = model$treedraws$trees,
      n_max = 1
    ),
    " "
  )[[1]]
  trees$nMCMC <- as.integer(firstLine[1])
  trees$nTree <- as.integer(firstLine[2])
  trees$nVar <- as.integer(firstLine[3])


  trees$structure <- dplyr::mutate(
    trees$structure,
    uniqueTreeID = cumsum(is.na(var) & is.na(cut) & is.na(leafValue)),
    iter = (uniqueTreeID - 1L) %/% trees$nTree + 1L,
    tree_id = (uniqueTreeID - 1L) %% trees$nTree + 1L,
    uniqueTreeID = NULL
  )

  # remove information about tree groups (was stored as missing lines)
  trees$structure <- dplyr::filter(trees$structure, stats::complete.cases(trees$structure))

  # children
  child_left <- function(nodes) {

    # must be grouped by iter and tree to apply
    pot_child <- nodes * 2L
    pot_child[!pot_child %in% nodes] <- NA_integer_

    return(pot_child)
  }

  child_right <- function(nodes) {

    # must be grouped by iter and tree to apply
    pot_child <- nodes * 2L + 1L
    pot_child[!pot_child %in% nodes] <- NA_integer_

    return(pot_child)
  }

  parent <- function(nodes) {
    parents <- nodes %/% 2L
    parents[parents == 0L] <- NA_integer_

    return(parents)
  }

  trees$structure <- dplyr::group_by(trees$structure, iter, tree_id)
  trees$structure <- dplyr::mutate(
    trees$structure,
    child_left = child_left(node),
    child_right = child_right(node)
  )


  # remove leaf info if no children
  trees$structure <- dplyr::mutate(
    dplyr::ungroup(trees$structure),
    isLeaf = is.na(child_left) & is.na(child_right),
    leafValue = ifelse(isLeaf, leafValue, NA_real_),
    cut = ifelse(isLeaf, NA_real_, cut), # is leaf, then no cut for stem
    var = ifelse(isLeaf, NA_character_, var), # is leaf, then no var for cut
    label = ifelse(
      isLeaf,
      as.character(round(leafValue, digits = 2)),
      paste(var, "<", round(cut, digits = 2))
    ),
    parent = parent(node)
  )

  # regroup
  trees$structure <- dplyr::select(
    dplyr::group_by(trees$structure, iter, tree_id),
    iter,
    tree_id,
    node,
    parent,
    label,
    tier,
    var,
    cut,
    isLeaf,
    leafValue,
    child_left,
    child_right
  )

  keepCols <- c("iter", "tree_id", "node", "parent", "label", "var")

  res <- dplyr::select(
    trees$structure,
    dplyr::one_of(keepCols)
  )

  # reorder
  res <- dplyr::select(res, -iter, -tree_id, dplyr::everything())

  # create edgle and node list
  res <- dplyr::mutate(res,
    new_node = seq_along(node),
    new_parent = new_node[match(parent, node)],
    node = new_node,
    parent = new_parent
  )
  res <- dplyr::select(res, -new_node, -new_parent)

  node_list <- dplyr::group_split(dplyr::select(res, -parent), keep = T)
  edge_list <- purrr::map(dplyr::group_split(dplyr::select(res, iter, tree_id, parent, node), keep = F), ~ dplyr::filter(., !is.na(parent)))

  tbl_graphList <- purrr::map2(
    .x = node_list,
    .y = edge_list,
    .f = ~ tidygraph::tbl_graph(
      nodes = .x,
      edges = .y,
      directed = T
    )
  )

  return(tbl_graphList)
}

# dbarts package ------------------------------------------------------------

plotPackageTrees.bart <- function(model, treeNum = NULL, iteration = NULL) {
  dbartsTrees(model)
}

dbartsTrees <- function(model, treeNum = NULL, iteration = NULL) {

  if(is.null(treeNum)){
    treeNum <- 1
  }
  if(is.null(iteration)){
    iteration <- 1
  }

    # Get trees
    trees <- model$fit$getTrees(treeNums = treeNum, chainNums = iteration)

    # Get variable names
    varNames <- colnames(model$fit$data@x)

    # set up data for plot
    trees$node <- c(1:(nrow(trees)))

    trees$value <- round(trees$value, 4)
    trees <- transform(trees, isLeaf = ifelse(var < 0, TRUE, FALSE))
    trees <- transform(trees, leafValue = ifelse(isLeaf == TRUE, trees$value, NA))
    trees <- transform(trees, cut = ifelse(isLeaf == FALSE, value, NA))
    trees <- transform(trees, var = ifelse(var < 0, NA, var))
    trees$var <- varNames[trees$var]
    trees <- transform(trees, label = ifelse(is.na(var), value, paste(var, value, sep = ' ≤ ')))

    # Create edges
    edges <- data.frame(
      parent = c("NEED TO FIND METHOD HERE!!!"),
      node   = c(2:(nrow(trees)))
    )

    # Turn into tbl_graph
    tbl_graphList <- tidygraph::tbl_graph(trees, edges)

    return(tbl_graphList)

    # Need to turn into a list of all trees!!
}

