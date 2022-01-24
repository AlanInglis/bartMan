#' proximityMatrixOLD
#'
#' @description Creates a matrix of proximity values.
#'
#' @param treeData A list of tree attributes created by treeData function.
#' @param data Data frame of data excluding the response.
#' @param nRows Number of rows to consider.
#' @param normalize Default is TRUE. Divide the total number of pairs of observations by
#' the number of trees.
#' @param reorder Default is TRUE. Whether to sort the matrix so high values are pushed to top left.
#'
#' @return A matrix containing proximity values.
#'
#' @importFrom utils combn
#' @export

proximityMatrixOLD <- function(tree,
                            data,
                            nRows,
                            normalize = TRUE,
                            reorder = TRUE){
  # get no of trees
  noTrees <- tree$nTree

  # get all indicies
  allIndex <- getIndex(tree = tree,
                       data = data,
                       nRows = seq_len(nRows))

  # get max value
  maxVal <- max(unlist(allIndex))

  # turn into matrix
  resMat <- tcrossprod(table(stack(allIndex)))

  if(normalize){
    resMat <- resMat/noTrees
  }

  diag(resMat) <- 0

  if(reorder){
    resMat <- proxReorder(resMat)
  }

  return(resMat)
}



# reorder proximity matrix ------------------------------------------------

proxReorder <- function(d) {

  prox <- as.dist(d)
  score <- apply(as.matrix(prox), 1, max)
  o <- DendSer::dser(-prox, -score, cost = DendSer::costLS)
  res <- d[o, o]
  class(res)<- class(d)
  res

}


# Function to get the index -----------------------------------------------

#' getIndex
#'
#' @description gets the indices of observations from leaf nodes.
#'
#' @param treeData A data frame created by treeData function.
#' @param data Data frame of data excluding the response.
#' @param nRows Sequence of integers from 1 to the number of rows in the data.
#'
#' @return A list containing vectors of the indices of observations from leaf nodes.
#'
#' @importFrom rrapply rrapply
#' @importFrom dplyr %>%
#' @importFrom purrr transpose
#' @importFrom purrr map2
#' @importFrom purrr map
#' @importFrom purrr keep
#' @importFrom purrr flatten
#' @importFrom stringr str_c



getIndex <- function(tree, data, nRows){

  listNodes <-  getAllNodeObsIndex(tree = tree,
                                   data = data,
                                   nRows = nRows)
  listNodes <- list(listNodes)

  out <- rrapply::rrapply(listNodes, how = 'bind')

  i1 <- grep('isLeaf', names(out))


  result <- purrr::map2(out[i1-2], out[i1], `[`) %>%
    transpose %>%
    purrr::map( ~ keep(.x, lengths(.x) > 0)) %>%
    flatten %>%
    setNames(str_c('index.', seq_along(.)))

  return(result)

}


# Get the observations ----------------------------------------------------


# This function extracts the node attributes for a single tree
# Additional info:
# tree = tree$structure
# var = data
# nRows = sequence of observations

observationNodes <- function(tree, var, nRows) {

  # small set up of tree list
  index <- nRows
  tree <- as.data.frame(tree)
  tree <- transform(tree, var = ifelse(is.na(var), -1, var))



  # Check if leaf node
  if (tree$var[1] == -1){
    return(list(value = tree$value[1],
                index = index,
                observationCount = length(index),
                isLeaf = TRUE))
  }


  # Check which observations on left
  leftSide <- var[index, tree$var[1]] <= tree$value[1]

  # get left side of tree
  left <- rebuildTree(tree[-1,], var, index[leftSide])



  if (is.null(left$n_nodes)) {
    noNodesLeft <- 1
  } else {
    noNodesLeft <- left$n_nodes
  }

  # get right side of tree
  right <- rebuildTree(tree[seq.int(2 + noNodesLeft, nrow(tree)), ], var, index[!leftSide])

  if (is.null(right$n_nodes)) {
    noNodesRight <- 1
  } else {
    noNodesRight <- right$n_nodes
  }


  # return list of attributes
  list(var = tree$var[1],
       value = tree$value[1],
       observationIndex = index,
       observationCount = nrow(var),
       left = left,
       right = right,
       n_nodes = 1 + noNodesLeft + noNodesRight)
}

# rebuild the tree
rebuildTree <- function(tree, var, nRows) {
  result <- observationNodes(tree, var, nRows)
  return(result)
}

# This function extracts all the node attributes for every tree
getAllNodeObsIndex <- function(tree, data, nRows) {

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










