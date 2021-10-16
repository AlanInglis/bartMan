
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








