#' Transform tree data into a structured dataframe
#'
#' This function takes raw data and a tree structure, then processes it to form a detailed and structured dataframe.
#' The data is transformed to indicate terminal nodes, calculate leaf values, and determine split values. It then
#' assigns labels, calculates node depth, and establishes hierarchical relationships within the tree.
#' Additional metadata about the tree, such as maximum depth, parent and child node relationships, and observation
#' nodes are also included. The final dataframe is organized and enriched with necessary attributes for further analysis.
#'
#' @param data A dataframe containing the raw data used for building the tree.
#' @param trees A dataframe representing the initial tree structure, including variables and values for splits.
#' @param response Optional character of the name of the response variable in your BART model. Including the response
#'                 will remove it from the list elements `Variable names` and `nVar`.
#'
#' @return A list containing a detailed dataframe of the tree structure (`structure`) with added information such as
#' node depth, parent and child nodes, and observational data, along with meta-information about the tree like
#' variable names (`varNames`), number of MCMC iterations (`nMCMC`), number of trees (`nTree`), and number of variables (`nVar`).
#'
#'
#' @export
#'
#' @importFrom dplyr bind_rows group_by mutate ungroup select row_number
#' @importFrom stats setNames

tree_dataframe <- function(data, trees, response = NULL){
  trees <- transform(trees, terminal = ifelse(is.na(var), TRUE, FALSE))
  trees <- transform(trees, leafValue = ifelse(terminal == TRUE, value, NA_integer_))
  trees <- transform(trees, splitValue = ifelse(terminal == FALSE, value, NA_integer_))
  trees$label <- ifelse(trees$terminal,
                        as.character(round(trees$leafValue, digits = 2)),
                        paste(trees$var, " \U2264 ", round(trees$splitValue, digits = 2)))

  depthList <- lapply(split(trees, ~treeNum + iteration),
                      function(x) cbind(x, depth = node_depth(x)-1))

  trees <- dplyr::bind_rows(depthList, .id = "list_id")

  # max depth
  trees <-  trees |>
    group_by(iteration, treeNum) |>
    mutate(depthMax = max(depth)) |>
    ungroup()

  # add node number
  trees <- trees |>
    dplyr::group_by(iteration, treeNum) |>
    dplyr::mutate(node = dplyr::row_number()) |>
    dplyr::ungroup()

  # get children and parent columns
  trees <-  trees |>
    group_by(iteration, treeNum, node)
  trees <- getChildren(data = trees)
  trees <-  trees |> ungroup()

  cat("Extracting Observation Data...\n")
  # get observations
  dat <- as.data.frame(data)
  trees <- getObservations(data = dat, treeData = trees)

  # add is stump column
  trees <-  trees  |>
    mutate(isStump = is.na(childLeft) & is.na(childRight) & is.na(parent) & depth == 0)

  # reordering the data and removing unnecessary columns
  trees <- dplyr::select(
    dplyr::group_by(trees, iteration, treeNum),
    var,
    splitValue,
    terminal,
    leafValue,
    iteration,
    treeNum,
    node,
    childLeft,
    childRight,
    parent,
    depth,
    depthMax,
    isStump,
    label,
    value,
    obsNode,
    noObs)|>
    ungroup()

  # attach other info
  Variable_names <- colnames(data)
  nMCMC <- max(trees$iteration)
  nTree <- max(trees$treeNum)
  nVar <- ncol(data)

  # remove response if selected
  if(!is.null(response)){
    if(!response %in% Variable_names){stop('Response name not found in data.')}
    Variable_names <- Variable_names[Variable_names != response]
    nVar <- nVar-1
  }
  df_tree_list <- list(structure = trees, varNames = Variable_names, nMCMC = nMCMC, nTree = nTree, nVar = nVar)

  hideHelper1 <- function(df){
    class(df) <- c("hideHelper1", class(df))
    df
  }
  trees <- hideHelper1(df_tree_list)
  return(trees)
}

