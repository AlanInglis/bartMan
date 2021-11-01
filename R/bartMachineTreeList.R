#' bartMachineTreeList
#'
#' @description Creates a list of tree attributes for a model
#' created by the bartMachine package.
#'
#' @param trees A data frame created by bartMachineTreeData function.
#' @return A list of every tree and its attributes.
#'
#' @importFrom tidygraph tbl_graph
#'
#' @export

bartMachineTreeList <- function(trees){

  df <- trees

  # split the dataframe into a list of dfs, one for each tree
  list_edges <- split(df, df$treeNum)

  # rename some columns
  x <- c("varID", "nodeID", "splitValue", "isLeaf", "leafValue",  "isStump", "var", "treeNum",
         "depth", "to",  "parentNode", "from", "iteration", "value", "label")
  list_edges <- lapply(list_edges, setNames, x)

  # reorder columns
  list_edges <- lapply(list_edges, function(x){
    x[c("var",
        "varID",
        "nodeID",
        "splitValue",
        "isLeaf",
        "leafValue",
        "isStump",
        "treeNum",
        "iteration",
        "depth",
        "parentNode",
        "value",
        "label",
        "from",
        "to")]})

  # remove unnessecery columns
  treesSplit <- lapply(list_edges, function(x) {
    x["isStump"] <- x["isLeaf"] <- x["parentNode"] <- x["value"] <- NULL; x
  })


  # create dataframe of edges
  dfOfEdges <- lapply(list_edges, function(df_tree){
    res <- data.frame(
      from = df_tree$from,
      to = df_tree$to
    )
    # delete NAs from result
    res <- res[!is.na(res$from),]
    return(res)
  })


  # Turn into a table graph object
  eachTree <- list()
  for (i in 1:length(treesSplit)) {
    eachTree[[i]] <- tidygraph::tbl_graph(treesSplit[[i]], dfOfEdges[[i]])
  }

  return(eachTree)

}
