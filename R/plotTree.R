#' plotTree
#'
#' @description plots individual trees from the R-packages BART, dbarts, and
#' bartMachine
#'
#' @param treeData A data frame created by treeData function
#' @param treeNo The tree number to plot.
#' @param iter The MCMC iteration or chain to plot.
#' @param plotType What type of plot to display. either dendrogram or icicle.
#' @return A plot of an individual tree
#'
#' @import ggraph
#' @import ggplot2
#' @importFrom tidygraph activate
#' @importFrom tidygraph tbl_graph
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr one_of
#' @importFrom dplyr group_split
#' @importFrom dplyr transmute
#' @importFrom dplyr row_number
#' @importFrom dplyr pull
#' @importFrom dplyr filter
#' @importFrom dplyr n
#' @importFrom stats setNames
#' @importFrom tibble as_tibble
#' @importFrom tidygraph bfs_dist
#'
#' @export

plotTree <- function(treeData,
                     iter = 1,
                     treeNo = 1,
                     plotType = c("dendrogram", "icicle")) {

  p <- plotAnyTree(treeData, iter = iter, treeNo = treeNo)

  if (plotType == "dendrogram") {
    gp <- ggraph(p, "dendrogram") +
      geom_edge_elbow() +
      geom_node_label(aes(label = label, color = label)) +
      theme_graph() +
      theme(legend.position = "none")
  } else if (plotType == "icicle") {
    gp <- ggraph(p, "partition") +
      geom_node_tile(aes(fill = var), size = 0.25) +
      geom_node_label(aes(label = label, color = var)) +
      scale_y_reverse() +
      theme_void() +
      theme(legend.position = "none")
  }

  return(gp)
}



# -------------------------------------------------------------------------


# Main plot function:
plotAnyTree <- function(treeData, iter = 1, treeNo = 1) {
  UseMethod("plotAnyTree")
}


# -------------------------------------------------------------------------



plotAnyTree.bart <- function(treeData,
                             treeNo = 1,
                             iter = 1) {
  df <- treeData$structure

  # select tree num from iteration
  df <- df %>%
    filter(iteration == iter, treeNum == treeNo)

  # Which columns to display
  keeps <- c("var", "node", "parent", "iteration", "treeNum", "label", "value")

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


  # Turn into a table graph object
  singleTree <- tidygraph::tbl_graph(nodes = nodeList[[1]], edges = edgeList[[1]])

  singleTree <- singleTree %>%
    activate(nodes) %>%
    mutate(depth = bfs_dist(root = 1))

  return(singleTree)
}


# DBARTS ------------------------------------------------------------------


plotAnyTree.dbarts <- function(treeData,
                               treeNo = 1,
                               iter = 1) {
  noObservations <- max(treeData$structure$n)

  # Which columns to display
  keeps <- c("var", "node", "isLeaf", "iteration", "treeNum", "label", "noObs", "value")

  treeData$structure <- dplyr::select(
    treeData$structure,
    dplyr::one_of(keeps)
  )

  treeData$structure <- transform(treeData$structure, varValue = ifelse(is.na(var), -1, var))


  # select tree num from iteration
  df <- treeData$structure %>%
    filter(iteration == iter, treeNum == treeNo)

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



  # get tree max depth and add to trees df
  dS <- as.vector(treeDepth(df))
  df$maxDepth <- dS[1]


  # Get Nodes and Edges -----------------------------------------------------

  # Get all the edges
  getEdges <- function(trees) {
    maxDepth <- trees$maxDepth[1]

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


  allEdges <- getEdges(df)


  # remove unnessecery columns
  df <- df %>%
    select(-varValue, -isLeaf, -noObs)

  # Turn into a table graph object
  singleTree <- tidygraph::tbl_graph(nodes = df, edges = allEdges)

  return(singleTree)
}




# BARTMACHINE -------------------------------------------------------------

plotAnyTree.bartMach <- function(treeData,
                                 treeNo = 1,
                                 iter = 1) {
  df <- treeData$structure

  # select tree num from iteration
  df <- df %>%
    filter(iteration == iter, treeNum == treeNo)

  # create dataframe of edges
  dfOfEdges <- data.frame(
    from = df$from,
    to = df$to
  )

  dfOfEdges <- na.omit(dfOfEdges)

  # remove unnecessary columns
  df <- df %>%
    select(-isStump, -to, -from, -node, -parentNode, -treeNumID)


  # Turn into a table graph object
  singleTree <- tidygraph::tbl_graph(nodes = df, edges = dfOfEdges)

  return(singleTree)
}
