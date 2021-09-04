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
#'
#'
#' @export

dbartsTreeData <- function(model) {

  # get all trees
  treesTotal <- model$call$ntree
  iteration  <- model$call$ndpost
  trees <- model$fit$getTrees(treeNums = 1:treesTotal, chainNums = 1:iteration)


  # Get variable names
  varNames <- colnames(model$fit$data@x)

  # set up data frame
  trees$node <- c(1:(nrow(trees)))
  trees$value <- round(trees$value, 4)
  trees <- transform(trees, isLeaf = ifelse(var < 0, TRUE, FALSE))
  trees <- transform(trees, leafValue = ifelse(isLeaf == TRUE, trees$value, NA_integer_))
  trees <- transform(trees, splitValue = ifelse(isLeaf == FALSE, value, NA_integer_))
  trees <- transform(trees, varName = ifelse(var < 0, NA, var))
  trees$varName <- varNames[trees$varName]
  trees <- transform(trees, label = ifelse(is.na(varName), value, paste(varName, value, sep = " ≤ ")))



  # fix node number for each tree
  noObservations <- length(model$y)
  trees <- trees %>%
    group_by(newSplit = cumsum(n == noObservations), .keep = FALSE) %>%
    mutate(node = 1:length(tree)) %>%
    ungroup() %>%
    mutate(var = varName) %>%
    rename(iteration = sample, treeNum = tree) %>%
    select(-newSplit, -.keep, - varName)

   # reorder columns
   trees <- trees %>%
     select(
     var,
     splitValue,
     node,
     isLeaf,
     leafValue,
     iteration,
     treeNum,
     label,
     n
   )

   # add class
  class(trees) <- c("tbl_df", "tbl", "data.frame", "dbarts")
  return(trees)
  }
