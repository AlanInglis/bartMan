#' treeDepth
#'
#' @description A plot of tree depth over iterations.
#'
#' @param treeData A list of tree attributes created using the extractTreeData function.
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
#' @export

treeDepth <- function(treeData) {
  maxIter <- max(treeData$structure$iteration)

  newTrees <- treeData$structure %>%
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
