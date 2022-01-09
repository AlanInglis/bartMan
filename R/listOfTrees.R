#' listOfTrees
#'
#' @description Creates a list of tree attributes for a model
#' created by either the BART, dbarts, or bartMachine packages.
#'
#' @param treeData A list of tree attributes created by extractTreeData function.
#' @return A list of every tree and its attributes.
#'
#' @importFrom dplyr select
#' @importFrom dplyr one_of
#' @importFrom dplyr mutate
#' @importFrom dplyr group_split
#' @importFrom dplyr filter
#' @importFrom purrr map
#' @importFrom purrr map2
#' @importFrom tidygraph tbl_graph
#'
#' @importFrom dplyr transmute
#' @importFrom dplyr row_number
#' @importFrom dplyr pull
#' @importFrom dplyr n
#' @importFrom stats setNames
#'
#' @export

listOfTrees <- function(treeData){

  listTrees <- createList(treeData)
  return(listTrees)
}


# -------------------------------------------------------------------------

# Main function:
createList <- function(treeData) {
  UseMethod("createList")
}


# BART --------------------------------------------------------------------

createList.bart <- function(treeData){


  # Which columns to display
  keeps <- c("var", "node", "parent", "iteration", "treeNum", "label", "value", "depth")

  res <- dplyr::select(
    treeData$structure,
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


createList.dbarts <- function(treeData){

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
    x["varValue"] <- x["isLeaf"]  <- NULL; x
  })


  # Turn into a table graph object
  eachTree <- list()
  for (i in 1:length(treesSplit)) {
    eachTree[[i]] <- tidygraph::tbl_graph(treesSplit[[i]], as.data.frame(allEdges[, i]))
  }



  return(eachTree)
}


# bartMachine -------------------------------------------------------------


createList.bartMach <- function(treeData){

  df <- treeData$structure

  # split the dataframe into a list of dfs, one for each tree
  list_edges <- split(df, df$treeNumID)

  # remove unnecessary columns
  treesSplit <- lapply(list_edges, function(x) {
    x["isStump"] <- x["to"] <- x['from'] <-  x["node"] <- x["parentNode"] <- x["treeNumID"] <- NULL; x
  })

  # create dataframe of edges
  dfOfEdges <- lapply(list_edges, function(df_tree){
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
