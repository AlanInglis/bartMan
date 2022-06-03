#' acceptRate
#'
#' @description Plots the acceptance rate of trees from a BART model.
#'
#' @param treeData A data frame created by treeData function.
#'  Displays a division on the plot to separate prior and post burn-in iterations.
#'
#' @return A plot of acceptance rate.
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

acceptRate <- function(treeData) {
  df <- treeData$structure

  maxIter <- max(df$iteration)

  acceptance <- df %>%
    filter(!is.na(var)) %>%
    group_by(iteration, treeNum) %>%
    summarise(values = paste0(sort(unique(label)), collapse = ",")) %>%
    group_by(treeNum) %>%
    mutate(changed = values != lag(values)) %>%
    replace_na(list(changed = TRUE)) %>%
    group_by(iteration) %>%
    summarise(percent_change = mean(changed))

  p <- ggplot(acceptance, aes(x = iteration, y = percent_change)) +
    geom_point(alpha = 0.5, colour = 'blue') +
    geom_smooth(formula = y ~ x, data = acceptance,
                method = "lm", colour = "black", se = F) +
    theme_bw() +
    xlab("Iteration") +
    ylab("% Acceptence Rate of Trees")


  return(p)
}
