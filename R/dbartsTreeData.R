#' dbartsTreeData
#'
#' @description Creates a data frame of all tree attributes for a model
#' created by the dbarts package.
#'
#' @param model Model created from the dbarts package.
#'
#' @return A data frame of every tree and its attributes.
#'
#'
#' @importFrom dplyr rename
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @importFrom dplyr coalesce
#'
#'
#' @export

dbartsTreeData <- function(model) {

  # get all trees
  treesTotal <- model$call$ntree
  iteration  <- model$call$ndpost

  trees <- list()
  trees$structure <- model$fit$getTrees(treeNums = 1:treesTotal, sampleNums = 1:iteration)

  # add other info
  trees$nMCMC <- as.integer(iteration)
  trees$nTree <- as.integer(treesTotal)
  trees$nVar  <- as.integer(length(model$varcount))

  # Get variable names
  varNames <- colnames(model$fit$data@x)

  # set up data frame
  trees$structure$node <- c(1:(nrow(trees$structure)))
  trees$structure$value <- round(trees$structure$value, 4)
  trees$structure <- transform(trees$structure, isLeaf = ifelse(var < 0, TRUE, FALSE))
  trees$structure <- transform(trees$structure, leafValue = ifelse(isLeaf == TRUE, value, NA_integer_))
  trees$structure <- transform(trees$structure, splitValue = ifelse(isLeaf == FALSE, value, NA_integer_))
  trees$structure <- transform(trees$structure, varName = ifelse(var < 0, NA, var))
  trees$structure$varName <- varNames[trees$structure$varName]
  trees$structure <- transform(trees$structure, label = ifelse(is.na(varName), value, paste(varName, value, sep = " ≤ ")))
  trees$structure <-  trees$structure %>%
    mutate(value = coalesce(splitValue, leafValue))

  # fix node number for each tree
  # noObservations <- length(model$yhat.train)
  # trees$structure <- trees$structure %>%
  #   group_by(newSplit = cumsum(n == noObservations), .keep = FALSE) %>%
  #   mutate(node = 1:length(tree)) %>%
  #   ungroup() %>%
  #   mutate(var = varName) %>%
  #   rename(iteration = sample, treeNum = tree) %>%
  #   select(-newSplit, -.keep, - varName)

  trees$structure <- trees$structure %>%
    group_by(tree, sample) %>%
    mutate(node = row_number()) %>%
    ungroup() %>%
    mutate(var = varName) %>%
    rename(iteration = sample, treeNum = tree) %>%
    select( - varName)

  # reorder columns
  trees$structure <- trees$structure %>%
    select(
      var,
      splitValue,
      node,
      isLeaf,
      leafValue,
      iteration,
      treeNum,
      label,
      value,
      n
    )

  # add class
  class(trees) <- c("list", "dbarts")
  return(trees)
}
