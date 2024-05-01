#' vimpBart
#'
#' @description A matrix with nMCMC rows with each variable as a column.
#' Each row represents an MCMC iteration. For each variable, the total count
#' of the number of times that variable is used in a tree is given.
#'
#' @param trees A data frame created by `extractTreeData` function.
#' @param type What value to return. Either the raw count 'val', the proportion 'prop',
#' the column means of the proportions 'propMean', or the median of the proportions 'propMedian'.
#'
#' @return A matrix of importance values
#'
#' @importFrom dplyr ungroup
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr select
#'
#' @examples
#' if(requireNamespace("dbarts", quietly = TRUE)){
#'  # Load the dbarts package to access the bart function
#'  library(dbarts)
#'  # Get Data
#'  df <- na.omit(airquality)
#'  # Create Simple dbarts Model For Regression:
#'  set.seed(1701)
#'  dbartModel <- bart(df[2:6], df[, 1], ntree = 5, keeptrees = TRUE, nskip = 10, ndpost = 10)
#'
#'  # Tree Data
#'  trees_data <- extractTreeData(model = dbartModel, data = df)
#'  vimpBart(trees_data, type = 'prop')
#'  }
#' @export
#'


vimpBart <- function(trees, type = 'prop'){

  if (!(type %in% c("val", "prop", "propMean", "propMedian"))) {
    stop("type must be \"val\", \"prop\", \"propMedian\", or \"propMean\"")
  }

  df <- trees$structure

  # get count of vars
  vCount <- df %>%
    ungroup() %>%
    select(iteration, var) %>%
    group_by(iteration) %>%
    table() %>%
    as.data.frame.matrix()

  # turn into matrix (with all values)
  nam <- trees$varName
 # mat <- matrix(0, nrow = trees$nMCMC, ncol = length(nam))
  mat <- matrix(0, nrow = length(unique(df$iteration)), ncol = length(nam))
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

