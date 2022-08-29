#' vimpBart
#'
#' @description A matrix with nMCMC rows with each variable as a column.
#' Each row represents an MCMC iteration. For each variable, the total count
#' of the number of times that variable is used in a tree is given.
#'
#' @param treeData A data frame created by treeData function.
#' @param type What value to return. Either the raw count 'val', the proportion 'prop',
#' the column means of the proportions 'propMean', or the median of the proportions 'propMedian'.
#'
#' @return A matrix
#'
#' @importFrom dplyr ungroup
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr select
#'
#' @export
#'


vimpBart <- function(treeData, type = 'prop'){

  if (!(type %in% c("val", "prop", "propMean", "propMedian"))) {
    stop("type must be \"val\", \"prop\", \"propMedian\", or \"propMean\"")
  }

  df <- treeData$structure

  # get count of vars
  vCount <- df %>%
    ungroup() %>%
    select(iteration, var) %>%
    group_by(iteration) %>%
    table() %>%
    as.data.frame.matrix()

  # turn into matrix (with all values)
  nam <- treeData$varName
  mat <- matrix(0, nrow = treeData$nMCMC, ncol = length(nam))
  colnames(mat) <- nam

  namesDat <- colnames(vCount)
  matchCol <- which(nam %in% namesDat)

  for (i in matchCol) {
    mat[, nam[i]] <- vCount[nam[i]][[nam[i]]]
  }

  if(type == 'prop'){
    mat <- proportions(mat, 1)
  }else if(type == 'propMean'){
    mat <- proportions(mat, 1)
    mat <- colMeans(mat)
  }else if(type == 'propMedian'){
    mat <- proportions(mat, 1)
    mat <- apply(mat, 2, median)
  }

  return(mat)
}

