#' treeData
#'
#' @description Creates a data frame of tree attributes for a model
#' created from either the bart, dbarts, or bartMachine packages.
#'
#' @param model A model created from either the bart, dbarts, or bartMachine packages.
#' @return A tibble of every tree and its attributes.
#'
#'
#' @export

treeData <- function(model){

  if(class(model) == "bart"){
    dataTrees <- dbartsTreeData(model)
  }else if(class(model) == "wbart" || class(model) == "lbart" || class(model) == "pbart"){
    dataTrees <- bartTreeData(model)
  }else if(class(model) == "bartMachine"){
    dataTrees <- bartMachineTreeData(model)
  }
  return(dataTrees)

}


