#' treeData
#'
#' @description Creates a data frame of tree attributes for a model
#' created from either the bart or dbarts packages.
#'
#' @param model A dmodel created from either the bart or dbarts packages.
#' @return A tibble of every tree and its attributes.
#'
#'
#' @export

treeData <- function(model){

  if(class(model) == "bart"){
    dataTrees <- dbartsTreeData(model)
  }else if(class(model) == "wbart"){
    dataTrees <- bartTreeData(model)
  }
  return(dataTrees)

}


