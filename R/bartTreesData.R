#' bartTreeData
#'
#' @description Creates a data frame of all tree attributes for a model
#' created by the bart package.
#'
#' @param model Model created from the bart package.
#'
#' @return A data frame of every tree and its attributes.
#'
#'
#' @importFrom readr read_table
#' @importFrom readr cols
#' @importFrom readr col_integer
#' @importFrom readr col_double
#' @importFrom purrr map_df
#' @importFrom dplyr tibble
#' @importFrom dplyr left_join
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @importFrom dplyr coalesce
#' @importFrom stats complete.cases
#' @importFrom dplyr one_of
#' @importFrom dplyr group_split
#' @importFrom dplyr filter
#' @importFrom purrr map
#' @importFrom purrr map2
#' @importFrom tidygraph tbl_graph
#'
#' @export



bartTreeData <- function(model) {


  # variable names:
  varNames <- names(model$varcount.mean)
  # get trees from model
  modelTrees <- model$treedraws$trees

  # extracting tree structure
  trees <- list()
  trees$structure <- suppressWarnings(
    readr::read_table(
      file = modelTrees,
      col_names = c("node", "var", "splitValue", "leafValue"),
      col_types =
        readr::cols(
          node = readr::col_integer(),
          var = readr::col_integer(),
          splitValue = readr::col_integer(),
          leaf = readr::col_double()
        ),
      skip = 1,
      na = c("")
    )
  )

  # Adding in columns
  trees$structure$var <- varNames[trees$structure$var + 1] # as vars are indexed at 0
  trees$structure$splitID <- trees$structure$splitValue + 1
  trees$structure$tier <- as.integer(floor(log2(trees$structure$node)))

  # getting split points
  splitPoints <- purrr::map_df(
    .x = model$treedraws$cutpoints,
    .f = ~ dplyr::tibble(splitValue = ., splitID = 1:length(.)),
    .id = "var"
  )

  # adding split points into tree structure
  trees$structure <- dplyr::left_join(
    dplyr::select(trees$structure, -splitValue),
    splitPoints,
    by = c("var", "splitID")
  )

  # Add in model fit info
  modelInfo <- unlist(strsplit(modelTrees, " "))[1:3]
  modelInfo <- gsub("(^\\d+)([\a-zA-Z0-9]*)", "\\1", modelInfo)

  trees$nMCMC <- as.integer(modelInfo[1])
  trees$nTree <- as.integer(modelInfo[2])
  trees$nVar  <- as.integer(modelInfo[3])

  trees$structure$uniqueTreeID <- cumsum(is.na(trees$structure$var) & is.na(trees$structure$splitValue) & is.na(trees$structure$leafValue))
  trees$structure$iteration <- ((trees$structure$uniqueTreeID - 1) %/% trees$nTree) + 1
  trees$structure$treeNum <- ((trees$structure$uniqueTreeID - 1) %% trees$nTree) + 1
  trees$structure$uniqueTreeID <- NULL

  # remove information about tree groups (i.e., rows with missing data)
  trees$structure$missingData <- complete.cases(trees$structure)
  missingIndex <- which(trees$structure$missingData == F)
  if (length(missingIndex == 0)) {
    trees$structure <- trees$structure[-missingIndex, ]
  }
  trees$structure$missingData <- NULL

  # Functions to get the left and right children nodes
  # and the parent nodes
  childLeft <- function(nodes) {
    childL <- nodes * 2
    childL[!childL %in% nodes] <- NA_integer_

    return(childL)
  }

  childRight <- function(nodes) {
    childR <- nodes * 2 + 1
    childR[!childR %in% nodes] <- NA_integer_

    return(childR)
  }

  parent <- function(nodes) {
    parents <- nodes %/% 2
    parents[parents == 0] <- NA_integer_

    return(parents)
  }

  trees$structure <- dplyr::group_by(trees$structure, iteration, treeNum)
  trees$structure <- dplyr::mutate(
    trees$structure,
    childLeft = childLeft(node),
    childRight = childRight(node)
  )

  trees$structure <- dplyr::ungroup(trees$structure)


  # Add is leaf column
  trees$structure$isLeaf <- is.na(trees$structure$childLeft) & is.na(trees$structure$childRight)
  # Remove leaf values for non-leaves
  trees$structure$leafValue <- ifelse(trees$structure$isLeaf, trees$structure$leafValue, NA_real_)
  # Remove split values for leaves
  trees$structure$splitValue <- ifelse(trees$structure$isLeaf, NA_real_, trees$structure$splitValue)
  # Remove var names for leaves
  trees$structure$var <- ifelse(trees$structure$isLeaf, NA_character_, trees$structure$var)
  # Add a label column
  trees$structure$label <- ifelse(trees$structure$isLeaf,
                                  as.character(round(trees$structure$leafValue, digits = 2)),
                                  paste(trees$structure$var, " ≤ ", round(trees$structure$splitValue, digits = 2))
  )
  # Add parent column
  trees$structure$parent <- parent(trees$structure$node)
  # Add value column
  trees$structure <- trees$structure %>%
    dplyr::mutate(value = dplyr::coalesce(splitValue, leafValue))

  # reordering the data and removing unnecessary columns
  trees$structure <- dplyr::select(
    dplyr::group_by(trees$structure, iteration, treeNum),
    var,
    splitValue,
    node,
    isLeaf,
    leafValue,
    childLeft,
    childRight,
    parent,
    iteration,
    treeNum,
    label,
    value,
    -splitID,
    -tier
  )

  # add class
  class(trees) <- c("list", "bart")
  return(trees)
}
