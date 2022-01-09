#' plotAllTrees
#'
#' @description Plots all the trees. By default, it shows the last iteration. If number of
#' trees is greater than 200, a sample of trees from that iteration will be shown.
#'
#' @param treeData A list of tree attributes created by treeDatafunction.
#' @param iter The selected iteration
#' @param treeNo The selected tree number.
#' @param sampleSize Sample the tree list.
#' @param cluster LOGICAL. If TRUE, then cluster by tree structures.
#'
#' @return A list containing vectors of the indices of observations from leaf nodes.
#'
#' @importFrom gridExtra grid.arrange
#' @importFrom ggraph ggraph
#' @importFrom igraph gsize
#' @importFrom cowplot get_legend
#' @importFrom cowplot plot_grid
#' @importFrom scales hue_pal
#' @importFrom dplyr select
#' @importFrom dplyr one_of
#' @importFrom dplyr mutate
#' @importFrom dplyr group_split
#' @importFrom dplyr filter
#' @importFrom purrr map
#' @importFrom purrr map2
#' @importFrom tidygraph tbl_graph
#' @importFrom dplyr transmute
#' @importFrom dplyr row_number
#' @importFrom dplyr pull
#' @importFrom dplyr n
#' @importFrom stats setNames
#'
#' @export

plotAllTrees <- function(treeData, iter = NULL, treeNo = NULL, cluster = FALSE) {
  allTrees <- plotAll(treeData, iter = iter, treeNo = treeNo)
  suppressWarnings(
    p <- plotAllTreesPlotFn(allTrees, cluster = cluster)
  )
  return(p)
}



# -------------------------------------------------------------------------


# Main plot function:
plotAll <- function(treeData, iter = NULL, treeNo = NULL) {
  UseMethod("plotAll")
}



# BART --------------------------------------------------------------------


plotAll.bart <- function(treeData, iter = NULL, treeNo = NULL) {
  maxIter <- treeData$nMCMC

  if (is.null(iter) & is.null(treeNo)) {
    df <- treeData$structure %>%
      filter(iteration == maxIter)
  } else if (is.null(iter) & !is.null(treeNo)) {
    df <- treeData$structure %>%
      filter(treeNum == treeNo)
  } else if (!is.null(iter) & is.null(treeNo)) {
    df <- treeData$structure %>%
      filter(iteration == iter)
  } else {
    df <- treeData$structure %>%
      filter(iteration == iter, treeNum == treeNo)
  }


  # Which columns to display
  keeps <- c("var", "node", "parent", "iteration", "treeNum", "label", "value", "depth")

  res <- dplyr::select(
    df,
    dplyr::one_of(keeps)
  )

  # Create edge and node list
  res <- dplyr::mutate(res,
    newNode = seq_along(node),
    newParent = newNode[match(parent, node)],
    node = newNode,
    parent = newParent
  )
  res <- dplyr::select(res, -newNode, -newParent)

  nodeList <- dplyr::group_split(dplyr::select(res, -parent), .keep = TRUE)
  edgeList <- purrr::map(
    dplyr::group_split(dplyr::select(
      res,
      iteration,
      treeNum,
      parent,
      node
    ), .keep = FALSE),
    ~ dplyr::filter(., !is.na(parent))
  )

  # Turn into data structure for tidy graph manipulation
  tblgList <- purrr::map2(
    .x = nodeList,
    .y = edgeList,
    .f = ~ tidygraph::tbl_graph(
      nodes = .x,
      edges = .y,
      directed = TRUE
    )
  )

  return(tblgList)
}


# dbarts ------------------------------------------------------------------


plotAll.dbarts <- function(treeData, iter = NULL, treeNo = NULL) {
  maxIter <- treeData$nMCMC

  if (is.null(iter) & is.null(treeNo)) {
    treeData$structure <- treeData$structure %>%
      filter(iteration == maxIter)
  } else if (is.null(iter) & !is.null(treeNo)) {
    treeData$structure <- treeData$structure %>%
      filter(treeNum == treeNo)
  } else if (!is.null(iter) & is.null(treeNo)) {
    treeData$structure <- treeData$structure %>%
      filter(iteration == iter)
  } else {
    treeData$structure <- treeData$structure %>%
      filter(iteration == iter, treeNum == treeNo)
  }

  noObservations <- max(treeData$structure$n)

  # Which columns to display
  keeps <- c("var", "node", "isLeaf", "iteration", "treeNum", "label", "n", "value", "depth")

  treeData$structure <- dplyr::select(
    treeData$structure,
    dplyr::one_of(keeps)
  )

  treeData$structure <- transform(treeData$structure, varValue = ifelse(is.na(var), -1, var))


  # split into individual trees
  treesSplit <- treeData$structure %>%
    group_split(cumsum(n == noObservations), .keep = FALSE)

  # add the depth of the tree
  treeDepth <- function(trees) {
    if (trees$varValue[1] == -1) {
      return(c(depth = 1, size = 1))
    }

    left <- treeDepth(trees[-1, , drop = FALSE])
    right <- treeDepth(trees[seq.int(2 + left[["size"]], nrow(trees)), , drop = FALSE])

    depthSize <- c(
      depth = 1 + max(left[["depth"]], right[["depth"]]),
      size = 1 + left[["size"]] + right[["size"]]
    )
    depthSize
  }

  # get every tree max depth and add to trees list
  dS <- as.vector(sapply(treesSplit, treeDepth)[1, ])
  treesSplit <- mapply(cbind, treesSplit, "depthNEW" = dS, SIMPLIFY = F)



  # Get Nodes and Edges -----------------------------------------------------

  # Get all the edges
  getEdges <- function(trees) {
    maxDepth <- trees$depthNEW[1]

    if (maxDepth > 2) {
      # Set the edge value:
      # this gets the 1st group of edges
      edgeSet1 <- do.call(rbind, lapply(
        split(trees, ~ cumsum(!isLeaf)),
        function(x) {
          with(x, expand.grid(from = node[!isLeaf], to = node[isLeaf]))
        }
      ))

      # get second group of edges
      edgeSet2 <- setNames(rev(data.frame(embed(trees$node[!trees$isLeaf], 2))), c("from", "to"))

      # bind them together
      edges <- rbind(edgeSet1, edgeSet2)

      # If any number in the to column appears more than twice then
      # replace with n-(number of times it appears-2)
      newFrom <- edges %>%
        group_by(from) %>%
        transmute(from := from - c(rep(0, 2), row_number())[row_number()]) %>%
        ungroup()

      # add to edges df
      edges$from <- pull(newFrom)

      # Change to right side of tree
      backToRoot <- edges %>%
        group_by(from) %>%
        filter(n() == 1)

      backToRoot$from <- 1
      edges <- rbind(edges, backToRoot)

      # remove duplicated rows and single row entries
      edges <- edges[!duplicated(edges), ]
      dups <- edges$from[duplicated(edges$from)]
      edges$unique <- edges$from %in% dups
      edgeIndex <- which(edges$unique == F)
      if (length(edgeIndex == 0)) {
        edges <- edges[-edgeIndex, ]
      }
      drops <- c("unique")
      edges <- edges[, !(names(edges) %in% drops)]

      # reorder
      "edges" <- `row.names<-`(edges[with(edges, order(from, to)), ], NULL)
      edges <- edges[order(edges$to), ]
    } else {
      edges <- do.call(rbind, lapply(
        split(trees, ~ cumsum(!isLeaf)),
        function(x) {
          with(x, expand.grid(from = node[!isLeaf], to = node[isLeaf]))
        }
      ))
    }
  }

  allEdges <- sapply(treesSplit, getEdges)

  # remove unnecessary columns
  treesSplit <- lapply(treesSplit, function(x) {
    x["varValue"] <- x["isLeaf"] <- NULL
    x
  })


  # Turn into a table graph object
  eachTree <- list()
  for (i in 1:length(treesSplit)) {
    eachTree[[i]] <- tidygraph::tbl_graph(treesSplit[[i]], as.data.frame(allEdges[, i]))
  }

  return(eachTree)
}

plotAll.bartMach <- function(treeData, iter = NULL, treeNo = NULL) {
  df <- treeData$structure

  maxIter <- treeData$MCMC

  if (is.null(iter) & is.null(treeNo)) {
    df <- df %>%
      filter(iteration == maxIter)
  } else if (is.null(iter) & !is.null(treeNo)) {
    df <- df %>%
      filter(treeNum == treeNo)
  } else if (!is.null(iter) & is.null(treeNo)) {
    df <- df %>%
      filter(iteration == iter)
  } else {
    df <- df %>%
      filter(iteration == iter, treeNum == treeNo)
  }


  # split the dataframe into a list of dfs, one for each tree
  list_edges <- split(df, df$treeNumID)

  # remove unnecessary columns
  treesSplit <- lapply(list_edges, function(x) {
    x["isStump"] <- x["to"] <- x["from"] <- x["node"] <- x["parentNode"] <- x["treeNumID"] <- NULL
    x
  })

  # create dataframe of edges
  dfOfEdges <- lapply(list_edges, function(df_tree) {
    res <- data.frame(
      from = df_tree$from,
      to = df_tree$to
    )
    # delete NAs from result
    res <- na.omit(res)
    return(res)
  })


  # Turn into a table graph object
  eachTree <- list()
  for (i in 1:length(treesSplit)) {
    eachTree[[i]] <- tidygraph::tbl_graph(treesSplit[[i]], dfOfEdges[[i]])
  }

  return(eachTree)
}


#' plotAllTreesplotFn
#'
#' @description This function is used to creat the plots.
#'
#' @param treeList A list of trees created by treeList function.
#' @param sampleSize Sample the tree list.
#' @param cluster LOGICAL. If TRUE, then cluster by tree structures.
#'
#' @return A list containing vectors of the indices of observations from leaf nodes.
#'
#'

plotAllTreesPlotFn <- function(treeList, sampleSize, cluster = FALSE) {

  # plot a sample of trees
  if (length(treeList) > 200) {
    sampleSize <- 200
    treeList <- sample(treeList, sampleSize, replace = FALSE)
  }

  # remove stumps
  treeList <- Filter(function(x) igraph::gsize(x) > 0, treeList)

  # set node colours
  nodenames <- unique(na.omit(unlist(lapply(treeList, . %>% activate(nodes) %>% pull(var)))))
  nodenames <- sort(nodenames)
  nodecolors <- setNames(scales::hue_pal(c(0, 360) + 15, 100, 64, 0, 1)(length(nodenames)), nodenames)

  # cluster by tree structure
  if (cluster) {
    indIDS <- map(treeList, function(x) {
      x %>%
        pull(var) %>%
        replace_na("a") %>%
        paste0(collapse = "")
    }) %>%
      unlist(use.names = F) %>%
      as_tibble() %>%
      mutate(ids = 1:n()) %>%
      group_by(value) %>%
      mutate(count = n():1) %>%
      arrange(value)

    ind <- indIDS %>%
      group_by(value) %>%
      mutate(valrank = max(count)) %>%
      ungroup() %>%
      arrange(-valrank, value, -count) %>%
      pull(ids)


    treeList <- treeList[ind]
  }

  allPlots <- lapply(treeList, plotFun, n = length(treeList), color = nodecolors)
  # get legend
  legend <- cowplot::get_legend(allPlots[[1]])
  # remove legends from individual plots
  allPlots <- lapply(allPlots, function(x) x + theme(legend.position = "none"))
  n <- length(allPlots)
  nRow <- floor(sqrt(n))
  allTreesPlot <- do.call("grid.arrange", c(allPlots, nrow = nRow))

  cowplot::plot_grid(allTreesPlot, legend, rel_widths = c(1, .1), ncol = 2)
}


plotFun <- function(List, colors = NULL, n) {
  plot <- ggraph(List, "partition") +
    geom_node_tile(aes(fill = var), size = 0.25) +
    geom_node_text(aes(label = ""), size = 4) +
    scale_y_reverse() +
    theme_void()
  if (!is.null(colors)) {
    plot <- plot + scale_fill_manual(values = colors, name = "Variable") +
      scale_color_manual(values = colors, na.value = "grey")
  }
  return(plot)
}
