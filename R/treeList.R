#' treeList
#'
#' @description Creates a list of tree attributes for a model
#' created from either the bart or dbarts packages.
#'
#' @param trees A data frame created by treeData function.
#' @return A list of every tree and its attributes.
#'
#'
#' @export

treeList <- function(trees){
  if(any(class(trees) == "bart")){
    listOfTrees <- bartTreeList(trees)
  }else if(any(class(trees) == "dbarts")){
    listOfTrees<- dbartsTreeList(trees)
  }else if(any(class(trees) == "bartMachine")){
    listOfTrees<- bartMachineTreeList(trees)
  }
  return(listOfTrees)
}




