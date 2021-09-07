#' dbartsTreeList
#'
#' @description Creates a list of tree attributes for a model
#' created by the dbarts package.
#'
#' @param trees A data frame created by dbartsTreeData function.
#' @return A list of every tree and its attributes.
#'
#' @importFrom dplyr select
#' @importFrom dplyr one_of
#' @importFrom dplyr group_split
#' @importFrom dplyr transmute
#' @importFrom dplyr row_number
#' @importFrom dplyr pull
#' @importFrom dplyr filter
#' @importFrom dplyr n
#' @importFrom stats setNames
#' @importFrom tidygraph tbl_graph
#'
#' @export


dbartsTreeList <- function(trees){

  noObservations <- max(trees$structure$n)

  # Which columns to display
  keeps <- c("var", "node", "isLeaf", "iteration", "treeNum", "label", "n", "value")

  trees$structure <- dplyr::select(
    trees$structure,
    dplyr::one_of(keeps)
  )

  trees$structure <- transform(trees$structure, varValue = ifelse(is.na(var), -1, var))


  # split into individual trees
  treesSplit <- trees$structure %>%
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
  treesSplit <- mapply(cbind, treesSplit, "depth" = dS, SIMPLIFY = F)


# Get Nodes and Edges -----------------------------------------------------

  # Get all the edges
  getEdges <- function(trees) {
    maxDepth <- trees$depth[1]

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

  # remove unnessecery columns
  treesSplit <- lapply(treesSplit, function(x) {
    x["depth"] <- x["varValue"] <- x["isLeaf"] <- x["n"] <- NULL; x
  })


  # Turn into a table graph object
  eachTree <- list()
  for (i in 1:length(treesSplit)) {
    eachTree[[i]] <- tidygraph::tbl_graph(treesSplit[[i]], as.data.frame(allEdges[, i]))
  }



  return(eachTree)
}
