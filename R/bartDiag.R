#' bartDiag
#'
#' @description Displays a selection of diagnostic plots for a BART model.
#'
#' @param model a model created from the BART package
#' @param response The name of the response for the fit.
#' @param burnIn Trace plot will only show iterations above selected burn in value.
#'
#' @return A selection of diagnostic plots
#'

#' @import ggplot2
#' @importFrom patchwork area
#' @importFrom patchwork plot_layout
#' @importFrom tidybayes residual_draws
#' @importFrom tidytreatment variance_draws
#' @importFrom tidybayes point_interval
#' @importFrom tidybayes geom_pointinterval
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr summarise
#' @importFrom dplyr tibble
#' @importFrom dplyr first
#' @importFrom tibble as_tibble
#' @export


bartDiag <- function(model,
                    response,
                    burnIn = 0
                    ){

  qq        <- bartQQ(model, response)
  trace     <- bartTrace(model, burnIn = burnIn)
  residual  <- bartResiduals(model, response = response)
  histogram <- bartHist(model)
  fitVSact  <- bartFitted(model, response)
  vImp      <- bartVimp(model)

  design <- c(
    area(1,1,3,3),
    area(1,5,3,7),
    area(5,1,7,3),
    area(5,5,7,7),
    area(9,1,11,3),
    area(9,5,11,7)
  )

  diagPlot <- qq + trace + residual + histogram +
    fitVSact + vImp + plot_layout(design = design)

  return(diagPlot)


}


# QQ plot -----------------------------------------------------------------

bartQQ <- function(model, response){

  res <- tidybayes::residual_draws(model, response = y, include_newdata = FALSE)

  p<- res %>% summarise(.residual = mean(.residual)) %>%
    ggplot(aes(sample = .residual)) +
    geom_qq(color = "blue", alpha = 0.2) +
    geom_qq_line() +
    xlab('Theoretical') +
    ylab("Sample") +
    theme_bw() +
    ggtitle("Q-Q plot")

  return(p)

}


# Trace plot --------------------------------------------------------------


bartTrace <- function(model, burnIn = 0){

  # get values
  varDraws <- tidytreatment::variance_draws(model, value = "siqsq")

  # filter values
  p <- varDraws %>%
    filter(.draw > burnIn) %>%
    ggplot(aes(x = .draw, y = siqsq)) +
    geom_line(color = "blue") +
    theme_bw() +
    xlab('Iteration') +
    ylab("SiqSq") +
    ggtitle("Trace plot of model variance")

  return(p)
}

# Residuals vs Fitted --------------------------------------------------------------


bartResiduals <- function(model,
                          response){

  res <- tidybayes::residual_draws(model, response = response, include_newdata = FALSE)

  p <- res %>%
    tidybayes::point_interval(.fitted, .residual, .width = c(0.95) ) %>%
    select(-.fitted.lower, -.fitted.upper) %>%
    ggplot() +
    tidybayes::geom_pointinterval(aes(x = .fitted, y = .residual, ymin = .residual.lower,  ymax = .residual.upper),
                                  alpha = 0.2,
                                  color = "blue") +
    theme_bw() +
    xlab('Fitted') +
    ylab("Residual") +
    ggtitle("Fitted vs Residuals")

  return(p)

}


# Histogram Residuals -----------------------------------------------------

bartHist <- function(model){

  res <- tidybayes::residual_draws(model, response = y, include_newdata = FALSE)

  p <- ggplot(data = res, aes(.residual)) +
  geom_histogram(bins = 50, color="blue", fill="white") +
    theme_bw()+
    xlab('Residual') +
    ylab("Frequency") +
    ggtitle("Histogram")

  return(p)

}


# Fitted Vs Actual --------------------------------------------------------


bartFitted <- function(model, response){

  res <- tidybayes::residual_draws(model, response = response, include_newdata = FALSE)


  # p <- res %>%
  # # tidybayes::point_interval(.fitted, y, .width = c(0.95) ) %>%
  #  # select(-y.lower, -y.upper) %>%
  #   summarise(.fitted = mean(.fitted), y = dplyr::first(y)) %>%
  #   ggplot(aes(x = y, y = .fitted)) +
  #   geom_point(color = "blue", alpha = 0.2) +
  #   geom_smooth(method = "lm", color = "black", formula = y ~ x) +
  #   xlab('Actual') +
  #   ylab("Fitted") +
  #   theme_bw() +
  #   ggtitle("Actual vs Fitted")


 p <-  res %>%
    tidybayes::point_interval(.fitted, y, .width = c(0.95) ) %>%
    select(-y.lower, -y.upper) %>%
    ggplot() +
    tidybayes::geom_pointinterval(aes(x = y, y = .fitted, ymin = .fitted.lower,  ymax = .fitted.upper),
                                  alpha = 0.2,
                                  color = "blue") +
    geom_smooth(aes(x = y, y = .fitted), method = "lm", color = "black", formula = y ~ x) +
    xlab('Actual') +
    ylab("Fitted") +
    theme_bw() +
    ggtitle("Actual vs Fitted")


  return(p)
}


# Variable Importance -----------------------------------------------------

bartVimp <- function(model)
{

  # get variable importance
  vImp <- model$varcount.mean
  vImp <- dplyr::tibble(Variable = names(vImp), Importance = vImp)


 p <- vImp %>%
    arrange(Importance) %>%
    mutate(Variable = factor(Variable, unique(Variable))) %>%
    ggplot() + aes(x=Variable, y=Importance) +
    geom_segment(aes(x=Variable, xend=Variable, y=0, yend=Importance), color="blue") +
    geom_point(color="blue") +
    theme_light() +
    coord_flip()+
    ggtitle(label = "BART Variable Importance") +
    theme_bw() +
    xlab('Variable') +
    ylab("Importance") +
    theme(axis.title.y = element_text(angle = 90, vjust = 0.5),
          legend.key.size = unit(0.5, "cm"))


  return(p)
}













