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
#' @param sizeNode Whether to size node width by the number of observations that fall into that node.
#'
#' @return A list containing vectors of the indices of observations from leaf nodes.
#'
#' @importFrom gridExtra arrangeGrob
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
#' @importFrom dplyr add_tally
#' @importFrom stats setNames
#'
#' @export

plotAllTrees <- function(treeData,
                         iter = NULL,
                         treeNo = NULL,
                         sampleSize = NULL,
                         cluster = NULL,
                         sizeNode = FALSE) {

  allTrees <- plotAll(treeData, iter = iter, treeNo = treeNo, cluster = cluster)

  suppressWarnings(
    p <- plotAllTreesPlotFn(allTrees, sampleSize = sampleSize, sizeNode = sizeNode)
  )
  return(p)
}



# -------------------------------------------------------------------------

#' @export
# Main plot function:
plotAll <- function(treeData, iter = NULL, treeNo = NULL, cluster = NULL) {
  UseMethod("plotAll")
}



# BART --------------------------------------------------------------------

#' @export
plotAll.bart <- function(treeData, iter = NULL, treeNo = NULL, cluster = NULL) {
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

  if (!is.null(cluster)) {
    if (cluster == "depth") {
      df <- df[with(df, order(-depth)), ]
    }
  }

  # Which columns to display
  keeps <- c("var", "node", "parent", "iteration", "treeNum", "label", "value", "depth", 'noObs')

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

  res$helper <- cumsum(is.na(res$parent))

  res <- res %>%
    ungroup() %>%
    group_by(helper)


  nodeList <- dplyr::group_split(dplyr::select(res, -parent), .keep = TRUE)

  nL <- map(nodeList, function(x) {
    x %>%
      select(-helper) %>%
      group_by(iteration, treeNum)
  })

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

  edgeList <- map(edgeList, function(x) {
    x %>%
      select(-iteration, -treeNum)
  })

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

  tblgList <- map(tblgList, function(x) {
    x %>%
      select(-helper)
  })

  if (!is.null(cluster)) {
    if (cluster == "var") {
      tblgList <- clusterTrees(tblgList)
    }
  }

  return(tblgList)
}


# dbarts ------------------------------------------------------------------
# -------------------------------------------------------------------------
#' @export
plotAll.dbarts <- function(treeData, iter = NULL, treeNo = NULL, cluster = NULL) {

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

  noObservations <- max(treeData$structure$noObs)

  treeData$structure <- treeData$structure %>%
    group_by(iteration, treeNum)

  # -------------------------------------------------------------------------

  # cluster trees
  if (!is.null(cluster)) {
    if (cluster == "depth") {
      df <- df[with(df, order(-depth)), ]
    }
  }

  # Which columns to display
  keeps <- c("var", "node", "isLeaf", "iteration", "treeNum", "label", "noObs", "value", "depth")

  treeData$structure <- dplyr::select(
    treeData$structure,
    dplyr::one_of(keeps)
  )


  treeData$structure <- transform(treeData$structure, varValue = ifelse(is.na(var), -1, var))
  treeData$structure <- treeData$structure %>%
    tibble() %>%
    group_by(iteration, treeNum)



  treesSplit <- treeData$structure %>%
    group_split(cumsum(noObs == noObservations), .keep = FALSE)

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

  if (!is.null(cluster)) {
    if (cluster == "var") {
      eachTree <- clusterTrees(eachTree)
    }
  }

  return(eachTree)
}



# bartMachine -------------------------------------------------------------
# -------------------------------------------------------------------------
#' @export
plotAll.bartMach <- function(treeData, iter = NULL, treeNo = NULL, cluster = NULL) {

  df <- treeData$structure
  maxIter <- treeData$MCMC
  noObservations <- max(treeData$structure$noObs)

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

  df <- df %>%
    group_by(iteration, treeNum)

  # -------------------------------------------------------------------------

  # cluster trees
  if (!is.null(cluster)) {
    if (cluster == "depth") {
      df <- df[with(df, order(-depth)), ]
    }
  }

  # split the dataframe into a list of dfs, one for each tree
  list_edges <- df %>%
    ungroup() %>%
    group_split(cumsum(noObs == noObservations), .keep = FALSE)


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

  if (!is.null(cluster)) {
    if (cluster == "var") {
      eachTree <- clusterTrees(eachTree)
    }
  }

  return(eachTree)
}


# -------------------------------------------------------------------------
# cluster function by variable --------------------------------------------
# -------------------------------------------------------------------------

clusterTrees <- function(treeList) {

  #df <- data
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

  return(treeList)
}

# -------------------------------------------------------------------------

#' plotAllTreesplotFn
#'
#' @description This function is used to creat the plots.
#'
#' @param treeList A list of trees created by treeList function.
#' @param sampleSize Sample the tree list.
#' @param sizeNode Whether to size node width by the number of observations that fall into that node.
#'
#'
#' @return A list containing vectors of the indices of observations from leaf nodes.
#'
#'

plotAllTreesPlotFn <- function(treeList, sampleSize = NULL, sizeNode = FALSE) {

  # plot a sample of trees
  if (length(treeList) > 200) {
    sampleSize <- 200
    treeList <- sample(treeList, sampleSize, replace = FALSE)
  } else if (!is.null(sampleSize)) {
    treeList <- sample(treeList, sampleSize, replace = FALSE)
  }

  # remove stumps
  treeList <- Filter(function(x) igraph::gsize(x) > 0, treeList)

  # set node colours
  nodenames <- unique(na.omit(unlist(lapply(treeList, . %>% activate(nodes) %>% pull(var)))))
  nodenames <- sort(nodenames)
  nodecolors <- setNames(scales::hue_pal(c(0, 360) + 15, 100, 64, 0, 1)(length(nodenames)), nodenames)

  allPlots <- lapply(treeList, plotFun, n = length(treeList), color = nodecolors)
  #allPlots <- lapply(treeList, plotFun, n = length(treeList), color = nodecolors, sizeNode = sizeNode)
  #allPlots <- lapply(treeList, function(x) plotFun(List = x, n = length(treeList), sizeNode = sizeNode, color = nodecolors))



  # get legend
  legend <- cowplot::get_legend(allPlots[[1]])
  # remove legends from individual plots
  allPlots <- lapply(allPlots, function(x) x + theme(legend.position = "none"))
  n <- length(allPlots)
  nRow <- floor(sqrt(n))
  allTreesPlot <- arrangeGrob(grobs=allPlots, nrow=nRow)

  cowplot::plot_grid(allTreesPlot, legend, rel_widths = c(1, .1), ncol = 2)
}


plotFun <- function(List, colors = NULL, n, sizeNode = FALSE) {

  if (sizeNode) {
  plot <- ggraph(List, "partition", weight = noObs) +
    geom_node_tile(aes(fill = var), size = 0.25) +
    geom_node_text(aes(label = ""), size = 4) +
    scale_y_reverse() +
    theme_void()
  if (!is.null(colors)) {
    plot <- plot + scale_fill_manual(values = colors, name = "Variable") +
      scale_color_manual(values = colors, na.value = "grey")
  }
} else {
  plot <- ggraph(List, "partition") +
    geom_node_tile(aes(fill = var), size = 0.25) +
    geom_node_text(aes(label = ""), size = 4) +
    scale_y_reverse() +
    theme_void()
  if (!is.null(colors)) {
    plot <- plot + scale_fill_manual(values = colors, name = "Variable") +
      scale_color_manual(values = colors, na.value = "grey")
  }
}
  return(plot)
}
