#' acceptRate
#'
#' @description Plots the acceptence rate of trees from a BART model.
#'
#' @param treeData A data frame created by treeData function.
#' @param burnIn Numerical value of the burn-in.
#'  Displays a division on the plot to separate prior and post burn-in iterations.
#'
#' @return A plot of acceptence rate.
#'
#' @importFrom dplyr filter
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr mutate
#' @importFrom dplyr lag
#' @importFrom tidyr replace_na
#' @import ggplot2
#'
#' @export

acceptRate <- function(treeData, burnIn = 0) {
  df <- treeData$structure

  maxIter <- max(df$iteration)

  acceptence <- df %>%
    filter(!is.na(var)) %>%
    group_by(iteration, treeNum) %>%
    summarise(values = paste0(sort(unique(label)), collapse = ",")) %>%
    group_by(treeNum) %>%
    mutate(changed = values != lag(values)) %>%
    replace_na(list(changed = TRUE)) %>%
    group_by(iteration) %>%
    summarise(percent_change = mean(changed))

  p <- ggplot(acceptence[1:burnIn, ], aes(iteration, percent_change)) +
    geom_vline(xintercept = burnIn, linetype = 5, alpha = 0.5) +
    geom_point(alpha = 0.5) +
    geom_smooth(formula = y ~ x, method = "lm", colour = "red", se = F) +
    geom_point(data = acceptence[burnIn:maxIter, ], colour = "blue", alpha = 0.5) +
    geom_smooth(formula = y ~ x, data = acceptence[burnIn:maxIter, ], method = "lm", colour = "black", se = F) +
    theme_bw() +
    xlab("Iteration") +
    ylab("% Acceptence Rate of Trees")

  return(p)
}