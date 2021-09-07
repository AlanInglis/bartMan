#' bartTreeList
#'
#' @description Creates a list of tree attributes for a model
#' created by the bart package.
#'
#' @param trees A data frame created by bartTreesData function.
#' @return A list of every tree and its attributes.
#'
#' @importFrom dplyr select
#' @importFrom dplyr one_of
#' @importFrom dplyr mutate
#' @importFrom dplyr group_split
#' @importFrom dplyr filter
#' @importFrom purrr map
#' @importFrom purrr map2
#' @importFrom tidygraph tbl_graph
#'
#' @export



bartTreeList <- function(trees) {


  # Which columns to display
  keeps <- c("var", "node", "parent", "iteration", "treeNum", "label", "value")

  res <- dplyr::select(
    trees$structure,
    dplyr::one_of(keeps)
  )

  # Create edge and node list
  res <- dplyr::mutate(res,
                       newNode = seq_along(node),
                       newParent = newNode[match(parent, node)],
                       node = newNode,
                       parent = newParent
  )
  res <- dplyr::select(res, -newNode, -newParent)

  nodeList <- dplyr::group_split(dplyr::select(res, -parent), .keep = TRUE)
  edgeList <- purrr::map(
    dplyr::group_split(dplyr::select(
      res,
      iteration,
      treeNum,
      parent,
      node
    ), .keep = FALSE),
    ~ dplyr::filter(., !is.na(parent))
  )

  # Turn into data structure for tidy graph manipulation
  tblgList <- purrr::map2(
    .x = nodeList,
    .y = edgeList,
    .f = ~ tidygraph::tbl_graph(
      nodes = .x,
      edges = .y,
      directed = TRUE
    )
  )

  return(tblgList)
}
