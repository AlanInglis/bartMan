#' treeNodes
#'
#' @description A plot of number of nodes over iterations.
#'
#' @param trees A list of tree attributes created using the extractTreeData function.
#'
#' @return A plot of tree number of nodes over iterations.
#'
#' @import ggplot2
#' @importFrom dplyr %>%
#' @importFrom dplyr summarize
#' @importFrom dplyr group_by
#'
#' @examples
#' \dontrun{
#' df_trees <- extractTreeData(model = my_model, data = my_data)
#' treeNodes(trees = df_trees)
#' }
#'
#' @export


treeNodes <- function(trees) {
  df <- trees$structure
  maxIter <- max(df$iteration)

  df <- df %>%
    group_by(iteration, treeNum) %>%
    summarize(count = n()) %>%
    group_by(iteration) %>%
    summarize(new = mean(count))

  p <- ggplot(df, aes(iteration, new)) +
    geom_point(alpha = 0.5, colour = 'blue') +
    # geom_line(alpha = 0.5, colour = 'blue') +
    geom_smooth(formula = y ~ x, method = "loess", colour = "black", se = F) +
    theme_bw() +
    xlab("Iteration") +
    ylab("Mean Tree Nodes")



  #
  # p <- ggplot(df, aes(iteration, new)) +
  #   geom_point(alpha = 0.5, colour = 'blue') +
  #   geom_line(alpha = 0.5, colour = 'blue') +
  #   theme_bw() +
  #   xlab("Iteration") +
  #   ylab("Average Tree Nodes")

  # p <- ggplot(df[1:burnIn, ], aes(iteration, new)) +
  #   geom_vline(xintercept = burnIn, linetype = 5, alpha = 0.5) +
  #   geom_point(alpha = 0.5) +
  #   geom_line(alpha = 0.5) +
  #   #geom_smooth(formula = y ~ x, method = "loess", colour = "red", se = F) +
  #   geom_point(data = df[burnIn:maxIter, ], colour = "blue", alpha = 0.5) +
  #   geom_line(data = df[burnIn:maxIter, ], colour = "blue", alpha = 0.5) +
  #   #geom_smooth(formula = y ~ x, data = df[burnIn:maxIter, ], color = "black", method = "loess", se = F) +
  #   theme_bw() +
  #   xlab("Iteration") +
  #   ylab("Average Tree Nodes")

  return(p)
}
