#' getIndex
#'
#' @description gets the indices of observations from leaf nodes.
#'
#' @param treeData A data frame created by treeData function.
#' @param data Data frame of data excluding the response.
#' @param nRows Sequence of integers from 1 to the number of rows in the data.
#'
#' @return A list containing vectors of the indices of observations from leaf nodes.
#'
#' @importFrom rrapply rrapply
#' @importFrom dplyr %>%
#' @importFrom purrr transpose
#' @importFrom purrr map2
#' @importFrom purrr map
#' @importFrom purrr keep
#' @importFrom purrr flatten
#' @importFrom stringr str_c



getIndex <- function(tree, data, nRows){

 listNodes <-  getAllNodeObsIndex(tree = tree,
                                  data = data,
                                  nRows = nRows)
 listNodes <- list(listNodes)

  out <- rrapply::rrapply(listNodes, how = 'bind')

  i1 <- grep('isLeaf', names(out))


  result <- purrr::map2(out[i1-2], out[i1], `[`) %>%
    transpose %>%
    purrr::map( ~ keep(.x, lengths(.x) > 0)) %>%
    flatten %>%
    setNames(str_c('index.', seq_along(.)))

  return(result)

}

