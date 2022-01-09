#' treeNodes
#'
#' @description A plot of number of nodes over iterations.
#'
#' @param treeData A list of tree attributes created using the extractTreeData function.
#' @param burnIn Selected burn in value.
#'
#' @return A plot of tree depth over iteration
#'
#' @import ggplot2
#' @importFrom dplyr %>%
#' @importFrom dplyr summarize
#' @importFrom dplyr group_by
#' @export


treeNodes <- function(treeData, burnIn) {
  df <- treeData$structure
  maxIter <- max(df$iteration)

  df <- df %>%
    group_by(iteration, treeNum) %>%
    summarize(count = n()) %>%
    group_by(iteration) %>%
    summarize(new = mean(count))

  p <- ggplot(df[1:burnIn, ], aes(iteration, new)) +
    geom_vline(xintercept = burnIn, linetype = 5, alpha = 0.5) +
    geom_point(alpha = 0.5) +
    geom_line(alpha = 0.5) +
    geom_smooth(formula = y ~ x, method = "lm", colour = "red", se = F) +
    geom_point(data = df[burnIn:maxIter, ], colour = "blue", alpha = 0.5) +
    geom_line(data = df[burnIn:maxIter, ], colour = "blue", alpha = 0.5) +
    geom_smooth(formula = y ~ x, data = df[burnIn:maxIter, ], color = "black", method = "lm", se = F) +
    theme_bw() +
    xlab("Iteration") +
    ylab("Average Tree Nodes")

  return(p)
}