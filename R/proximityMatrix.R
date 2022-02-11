#' proximityMatrix
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
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom DendSer dser
#' @export

proximityMatrix <- function(treeData, data, nRows, normalize = TRUE, reorder = TRUE) {

  data <- data[nRows, ]
  # get observations in each node
  dfObs <- treeData$structure %>%
    select(var, splitValue, iteration, treeNum, value, noObs) %>%
    group_by(iteration, treeNum) %>%
    mutate(obsList = evalNode(data, var, splitValue))

  obsIndex <- lapply(dfObs$obsList, function(x) {
    lapply(x, row.names)
  })

  whichObs <- lapply(obsIndex, rapply, f = c)
  whichObs <- lapply(whichObs, as.numeric)

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
  resMat <- tcrossprod(table(stack(allIndex)))
  diag(resMat) <- 0

  # normalize
  treeNumber <- treeData$nTree
  iterNumber <- treeData$nMCMC
  treeTotal <- treeNumber*iterNumber
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
  prox <- as.dist(d)
  score <- apply(as.matrix(prox), 1, max)
  o <- DendSer::dser(-prox, -score, cost = DendSer::costLS)
  res <- d[o, o]
  class(res) <- class(d)
  res
}
