#' proximityMatrix
#'
#' @description Creates a matrix of proximity values.
#'
#' @param treeData A data frame created by treeData function.
#' @param data Data frame of data excluding the response.
#' @param nRows Sequence of integers from 1 to the number of rows in the data.
#' @param normalize Default is TRUE. Divide the total number of pairs of observations by
#' the number of trees.
#' @param reorder Default is TRUE. Whether to sort the matrix so high values are pushed to top left.
#'
#' @return A matrix containing proximity values.
#'
#' @importFrom utils combn
#' @export

proximityMatrix <- function(tree = tree,
                            data = data,
                            nRows = nRows,
                            normalize = TRUE,
                            reorder = TRUE){
  # get no of trees
  noTrees <- tree$nTree

  # get all indicies
  allIndex <- getIndex(tree = tree,
                       data = data,
                       nRows = nRows)
  # get max value
  maxVal <- max(unlist(allIndex))

  # turn into matrix
  resMat <- Reduce(`+`, lapply(allIndex, function(x)
    crossprod(table(transform(stack(setNames(
      combn(x,
            2, simplify = FALSE), combn(x, 2, paste, collapse = "_"))),
      values = factor(values, levels = 1:maxVal))[2:1]))> 0))

  if(normalize){
    resMat <- resMat/noTrees
  }

  diag(resMat) <- 0

  if(reorder){
    resMat <- proxReorder(resMat)
  }

  return(resMat)
}
