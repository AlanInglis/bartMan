#' treeDepth
#'
#' @description A plot of tree depth over iterations.
#'
#' @param trees A list of tree attributes created using the extractTreeData function.
#'
#' @return A plot of average tree depths over iteration
#'
#' @import ggplot2
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr summarize
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom dplyr as_tibble
#'
#' @examples
#' if(requireNamespace("dbarts", quietly = TRUE)){
#'  # Load the dbarts package to access the bart function
#'  library(dbarts)
#'  # Get Data
#'  df <- na.omit(airquality)
#'  # Create Simple dbarts Model For Regression:
#'  set.seed(1701)
#'  dbartModel <- bart(df[2:6], df[, 1], ntree = 5, keeptrees = TRUE, nskip = 10, ndpost = 10)
#'
#'  # Tree Data
#'  trees_data <- extractTreeData(model = dbartModel, data = df)
#'  treeDepth(trees = trees_data)
#' }
#'
#' @export

treeDepth <- function(trees) {
  maxIter <- max(trees$structure$iteration)

  newTrees <- trees$structure %>%
    select(iteration, treeNum, depthMax) %>%
    as_tibble() %>%
    group_by(iteration, treeNum) %>%
    summarize(maxDepth = max(depthMax)) %>%
    ungroup() %>%
    group_by(iteration) %>%
    summarize(avgDepth = mean(maxDepth))

  ylimMax <- max(newTrees$avgDepth)

  p <-  ggplot(newTrees, aes(iteration, avgDepth)) +
    geom_point(alpha = 0.5, colour = 'blue') +
    #geom_line(alpha = 0.5, colour = 'blue') +
    geom_smooth(formula = y ~ x, method = "loess", colour = "black", se = F) +
    theme_bw() +
    xlab("Iteration") +
    ylab("Mean Tree Depth")



  # p<- ggplot(newTrees, aes(iteration, avgDepth)) +
  #   geom_point(alpha = 0.5, colour = 'blue') +
  #   geom_line(alpha = 0.5, colour = 'blue') +
  #   theme_bw() +
  #   xlab("Iteration") +
  #   ylab("Average Tree Depth")


  # p <- ggplot(newTrees[1:burnIn, ], aes(iteration, avgDepth)) +
  #   geom_vline(xintercept = burnIn, linetype = 5, alpha = 0.5) +
  #   geom_point(alpha = 0.5) +
  #   geom_line(alpha = 0.5) +
  #   #geom_smooth(formula = y ~ x, method = "loess", colour = "red", se = F) +
  #   geom_point(data = newTrees[burnIn:maxIter, ], colour = "blue", alpha = 0.5) +
  #   geom_line(data = newTrees[burnIn:maxIter, ], colour = "blue", alpha = 0.5) +
  #   #geom_smooth(formula = y ~ x, data = newTrees[burnIn:maxIter, ], color = "black", method = "loess", se = F) +
  #   theme_bw() +
  #   xlab("Iteration") +
  #   ylab("Average Tree Depth")


  return(p)
}
