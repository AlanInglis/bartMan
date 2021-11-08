#' bartDiag
#'
#' @description Displays a selection of diagnostic plots for a BART model.
#'
#' @param model a model created from the BART package
#' @param response The name of the response for the fit.
#' @param burnIn Trace plot will only show iterations above selected burn in value.
#' @param pal A vector of colours to show importance, for use with scale_fill_gradientn.
#' @param impLims Specifies the fit range for the color map for importance.
#' @param reorder If TRUE, then the variable importance plot is reordered in terms of importance.
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
#' @importFrom tibble as_tibble
#' @export


bartDiag <- function(model,
                    response,
                    burnIn = 0,
                    pal = rev(colorspace::sequential_hcl(palette = "Blues 3", n = 100)),
                    impLims = NULL,
                    reorder = FALSE
                    ){

  qq <- bartQQ(model, response)
  trace <- bartTrace(model, burnIn = burnIn)
  residual <- bartResiduals(model, response = response, pal = pal)
  histogram <- bartHist(model)
  fitVSact <- bartFitted(model)
  vImp <- bartVimp(model,
                   impLims = impLims,
                   pal = pal,
                   reorder = reorder)

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
                          response,
                          pal = rev(colorspace::sequential_hcl(palette = "Reds 3", n = 100))){

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
    ggtitle("Residuals vs Fitted")

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


bartFitted <- function(model){

  res <- tidybayes::residual_draws(model, response = y, include_newdata = FALSE)

  p <- res %>%
    summarise(.fitted = mean(.fitted), y = first(y)) %>%
    ggplot(aes(x = y, y = .fitted)) +
    geom_point(color = "blue", alpha = 0.2) +
    geom_smooth(method = "lm", color = "black") +
    xlab('Actual') +
    ylab("Fitted") +
    theme_bw() + ggtitle("Observations vs Fitted")

  return(p)
}


# Variable Importance -----------------------------------------------------

bartVimp <- function(model,
                     impLims = NULL,
                     pal = rev(colorspace::sequential_hcl(palette = "Reds 3", n = 100)),
                     reorder = FALSE)
{

  # get variable importance
  vImp <- model$varcount.mean
  vImp <- dplyr::tibble(Variable = names(vImp), Importance = vImp)

  # set limits
  if (is.null(impLims)) {
    impLims <- range(vImp$Importance)
    limitsImp <- range(labeling::rpretty(impLims[1], impLims[2]))
  } else {
    limitsImp <- impLims
  }

  # reorder plots
  if(reorder){
    vImp$Variable <- forcats::fct_rev(forcats::fct_reorder(vImp$Variable, vImp$Importance))
  }

  # plot
  p <- ggplot(vImp, aes(x = Variable, y = Importance)) +
    geom_col(aes(fill = Importance)) +
    scale_fill_gradientn(
      colors = pal, limits = limitsImp, name = "Vimp",
      guide = guide_colorbar(
        frame.colour = "black",
        ticks.colour = "black"
      ), oob = scales::squish
    ) +
    ggtitle(label = "BART Variable Importance") +
    theme_bw() +
    xlab('Variable') +
    ylab("Inclusion") +
    theme(axis.title.y = element_text(angle = 90, vjust = 0.5),
          legend.key.size = unit(0.5, "cm"))

  return(p)
}













