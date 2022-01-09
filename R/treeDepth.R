#' treeDepth
#'
#' @description A plot of tree depth over iterations.
#'
#' @param treeData A list of tree attributes created using the extractTreeData function.
#' @param burnIn Selected burn in value.
#'
#' @return A plot of average tree depths over iteration
#'
#' @import ggplot2
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr summarize
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom tibble as_tibble
#' @export

treeDepth <- function(treeData, burnIn) {
  maxIter <- max(treeData$structure$iteration)

  newTrees <- treeData$structure %>%
    select(iteration, treeNum, depth) %>%
    as_tibble() %>%
    group_by(iteration, treeNum) %>%
    summarize(maxDepth = max(depth)) %>%
    ungroup() %>%
    group_by(iteration) %>%
    summarize(avgDepth = mean(maxDepth))

  ylimMax <- max(newTrees$avgDepth)

  p <- ggplot(newTrees[1:burnIn, ], aes(iteration, avgDepth)) +
    geom_vline(xintercept = burnIn, linetype = 5, alpha = 0.5) +
    geom_point(alpha = 0.5) +
    geom_line(alpha = 0.5) +
    geom_smooth(formula = y ~ x, method = "lm", colour = "red", se = F) +
    geom_point(data = newTrees[burnIn:maxIter, ], colour = "blue", alpha = 0.5) +
    geom_line(data = newTrees[burnIn:maxIter, ], colour = "blue", alpha = 0.5) +
    geom_smooth(formula = y ~ x, data = newTrees[burnIn:maxIter, ], color = "black", method = "lm", se = F) +
   # ylim(0, ylimMax) +
    theme_bw() +
    xlab("Iteration") +
    ylab("Average Tree Depth")


  return(p)
}