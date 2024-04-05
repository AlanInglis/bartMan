#' proximityMatrix
#'
#' @description Creates a matrix of proximity values.
#'
#' @param trees A list of tree attributes created by `extractTreeData` function.

#' @param nRows Number of rows to consider.
#' @param normalize Default is TRUE. Divide the total number of pairs of observations by
#' the number of trees.
#' @param reorder Default is TRUE. Whether to sort the matrix so high values are pushed to top left.
#' @param iter Which iteration to use, if NULL the proximity matrix is calculated over all
#' iterations.
#'
#' @return A matrix containing proximity values.
#'
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom DendSer dser
#'
#' @examples
#' \donttest{
#' df_trees <- extractTreeData(model = my_model, data = my_data)
#' bmProx <- proximityMatrix(trees = df_trees, reorder = TRUE, normalize = TRUE, iter = 1)
#' }
#' @export

proximityMatrix <- function(trees, nRows, normalize = TRUE, reorder = TRUE, iter = NULL) {


  if (!is.null(iter)) {
    trees$structure <- trees$structure %>%
      filter(iteration == iter)
    treeTotal <- max(trees$structure$treeNum)
  } else {
    treeNumber <- trees$nTree
    iterNumber <- trees$nMCMC
    treeTotal <- treeNumber*iterNumber
  }

  #data <- data[nRows, ]
  # get observations in each node

  # dfObs <- trees$structure %>%
  #   select(var, splitValue, iteration, treeNum, value, noObs) %>%
  #   group_by(iteration, treeNum) %>%
  #   mutate(obsList = evalNode(dat, var, splitValue))
  #
  # obsIndex <- lapply(dfObs$obsList, function(x) {
  #   lapply(x, row.names)
  # })
  #
  # whichObs <- lapply(obsIndex, rapply, f = c)
  # whichObs <- lapply(whichObs, as.numeric)


  dfObs <- trees$structure %>%
       select(var, splitValue, iteration, treeNum, value, noObs)
  whichObs <- trees$structure$obsNode

  # get all indicies
  idxLength <- length(whichObs)

  # rename the list elements to work with stack function
  myLetters <- function(n) {
    unlist(Reduce(paste0,
                  replicate(n %/% length(letters), letters, simplify = FALSE),
                  init = letters,
                  accumulate = TRUE
    ))[1:n]
  }

  myNames <- myLetters(idxLength)
  allIndex <- setNames(whichObs, myNames)

  #add back to dataframe to find terminal nodes
  dfObs$whichNode <- allIndex

  # filter terminal nodes
  dfTerm <- dfObs %>%
    filter(is.na(var))

  # get observations
  allIndex <- dfTerm$whichNode

  # turn into matrix
  resMat <- tcrossprod(table(utils::stack(allIndex)))
  diag(resMat) <- 0

  # normalize
  if (normalize) {
    resMat <- resMat / treeTotal
  }

  if (reorder) {
    resMat <- proxReorder(resMat)
  }


  return(resMat)
}

# reorder proximity matrix ------------------------------------------------

proxReorder <- function(d) {
  prox <- stats::as.dist(d)
  score <- apply(as.matrix(prox), 1, max)
  o <- DendSer::dser(-prox, -score, cost = DendSer::costLS)
  res <- d[o, o]
  class(res) <- class(d)
  res
}
