#' treeMap
#'
#' @description Creates a treemap displaying the frequency of different tree structures.
#'
#' @param treeList A list of trees created using the treeList function.
#'
#' @return A treemap plot.
#'
#'
#' @import ggplot2
#' @importFrom purrr map_chr
#' @importFrom tidygraph activate
#' @importFrom tidyr separate
#' @importFrom dplyr as_tibble
#' @importFrom dplyr arrange
#' @importFrom treemapify geom_treemap
#' @importFrom treemapify geom_treemap_text
#'
#'
#' @export


treeMap <- function(treeList){

  # get edge frequency
  edgeFreq <- treeList %>%
    map_chr(~as_tibble(activate(.x, edges)) %>%
              map_chr(str_c, collapse = " ") %>%
              toString())%>%
    table() %>%
    as_tibble() %>%
    setNames(c("data", "frequency")) %>%
    separate(data, c("from", "to"), ", ") %>%
    arrange(-frequency)

  gp <- ggplot(edgeFreq,
         aes(fill = frequency,
             area = frequency,
             label = to)) +
    treemapify::geom_treemap() +
    treemapify::geom_treemap_text(colour = "white",
                      place = "centre")

  return(gp)

}
