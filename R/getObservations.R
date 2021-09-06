
dbartsSampleNodes <- function(tree, var, index) {

  # Check if leaf node
  if (tree$var[1] == -1){
    return(list(value = tree$value[1], index = index))
  }

  # add is leaf
  tree$isLeaf <- ifelse(tree$var < 0, TRUE, FALSE)

  # add node number
  tree$nodeNo <- c(1:nrow(tree))

  # get leaf index
  leafIndex <- which(tree$isLeaf)
  parentIndex <- which(tree$isLeaf == F)

  # Check observations on left
  leftSide <- var[index, tree$var[1]] <= tree$value[1]

  # get left side of tree
  left <- rebuildTree(tree[-1,], var, index[leftSide])

  if (is.null(left$n_nodes)) {
    noNodesLeft <- 1
  } else {
    noNodesLeft <- left$n_nodes
    left$noNodesLeft <- NULL
  }

  # get right side of tree
  right <- rebuildTree(tree[seq.int(2 + noNodesLeft, nrow(tree)), ], var, index[!leftSide])

  if (is.null(right$n_nodes)) {
    noNodesRight <- 1
  } else {
    noNodesRight <- right$n_nodes
    right$noNodesRight <- NULL
  }


  # return list of attributes
  list(var = tree$var[1],
       value = tree$value[1],
       index = index,
       left = left,
       right = right,
       n_nodes = 1 + noNodesLeft + noNodesRight,
       isLeaf = tree$isLeaf)
}

# rebuild the tree
rebuildTree <- function(tree, var, index) {
  result <- dbartsSampleNodes(tree, var, index)
  #result$n_nodes <- NULL
  return(result)
}
