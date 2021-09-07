#' getNodeObsIndex
#'
#' @description Gets the index of the observations contained within a node.
#'
#' @param treeData A data frame created by treeData function.
#' @param data Data frame of data excluding the response.
#' @param nRows Sequence of integers from 1 to the number of rows in the data.
#'
#' @return A list containing node attributes.
#'
#' @export


getNodeObsIndex <- function(tree, data, nRows) {

  # small set up of tree
  index <- nRows
  tree  <- as.data.frame(tree$structure)
  tree  <- transform(tree, var = ifelse(is.na(var), -1, var))

  rebuiltTreesAll <- by(
    tree[c("treeNum", "var", "value")],
    tree[c("iteration")],
    function(trees) {
      by(tree[c("var", "value")], tree[["treeNum"]],
         rebuildTree, var = data, nRows = index)
    }
  )
}
