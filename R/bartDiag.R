#' bartDiag
#'
#' @description Displays a selection of diagnostic plots for a BART model.
#'
#' @param model a model created from either the BART, dbarts, or bartMachine package.
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
                     burnIn = 0) {
  qq <- bartQQ(model, response)
  trace <- bartTrace(model, burnIn = burnIn)
  residual <- bartResiduals(model, response = response)
  histogram <- bartHist(model, response)
  fitVSact <- bartFitted(model, response)
  vImp <- bartVimp(model)

  design <- c(
    area(1, 1, 3, 3),
    area(1, 5, 3, 7),
    area(5, 1, 7, 3),
    area(5, 5, 7, 7),
    area(9, 1, 11, 3),
    area(9, 5, 11, 7)
  )

  diagPlot <- qq + trace + residual + histogram +
    fitVSact + vImp + plot_layout(design = design)

  return(diagPlot)
}


# QQ plot -----------------------------------------------------------------

bartQQ <- function(model, response) {
  if (class(model) == "wbart" || class(model) == "bartMachine") {
    res <- tidybayes::residual_draws(model, response = response, include_newdata = FALSE)
    res <- res %>% summarise(.residual = mean(.residual))
  } else {
    res <- data.frame(.residual = residuals(model))
  }

  p <- res %>%
    ggplot(aes(sample = .residual)) +
    geom_qq(color = "blue", alpha = 0.2) +
    geom_qq_line() +
    xlab("Theoretical") +
    ylab("Sample") +
    theme_bw() +
    ggtitle("Q-Q plot")

  return(p)
}


# Trace plot --------------------------------------------------------------


bartTrace <- function(model, burnIn = 0) {
  if (class(model) == "wbart" || class(model) == "bartMachine") {
    # get values
    varDraws <- tidytreatment::variance_draws(model, value = "siqsq")
    varDraws$sigma <- sqrt(varDraws$siqsq)
  } else {
    varDraws <- data.frame(
      sigma = model$sigma,
      .draw = c(1:model$call$ndpost)
    )
  }


  # filter values
  p <- varDraws %>%
    filter(.draw > burnIn) %>%
    ggplot(aes(x = .draw, y = sigma)) +
    geom_line(color = "blue") +
    theme_bw() +
    xlab("Iteration") +
    ylab("Sigma") +
    ggtitle("Trace plot of model variance")

  return(p)
}

# Residuals vs Fitted --------------------------------------------------------------


bartResiduals <- function(model,
                          response) {
  if (class(model) == "wbart" || class(model) == "bartMachine") {
    res <- tidybayes::residual_draws(model, response = response, include_newdata = FALSE)
  } else {
    res <- data.frame(
      .residual = residuals(model),
      .fitted = fitted(model)
    )
  }


  if (class(model) == "wbart" || class(model) == "bartMachine") {
    p <- res %>%
      tidybayes::point_interval(.fitted, .residual, .width = c(0.95)) %>%
      select(-.fitted.lower, -.fitted.upper) %>%
      ggplot() +
      tidybayes::geom_pointinterval(aes(x = .fitted, y = .residual, ymin = .residual.lower, ymax = .residual.upper),
        alpha = 0.2,
        color = "blue"
      ) +
      theme_bw() +
      xlab("Fitted") +
      ylab("Residual") +
      ggtitle("Fitted vs Residuals")
  } else {
    p <- ggplot(res, aes(.fitted, .residual)) +
      geom_point(color = "blue", alpha = 0.2) +
      theme_bw() +
      xlab("Fitted") +
      ylab("Residual") +
      ggtitle("Fitted vs Residuals")
  }

  return(p)
}


# Histogram Residuals -----------------------------------------------------

bartHist <- function(model, response) {

  if (class(model) == "wbart" || class(model) == "bartMachine") {
    res <- tidybayes::residual_draws(model, response = response, include_newdata = FALSE)
  } else {
    res <- data.frame(.residual = residuals(model))
  }

  p <- ggplot(data = res, aes(.residual)) +
    geom_histogram(bins = 50, color = "blue", fill = "white") +
    theme_bw() +
    xlab("Residual") +
    ylab("Frequency") +
    ggtitle("Histogram")

  return(p)
}


# Fitted Vs Actual --------------------------------------------------------


bartFitted <- function(model, response) {


  if (class(model) == "wbart" || class(model) == "bartMachine") {
    res <- tidybayes::residual_draws(model, response = response, include_newdata = FALSE)
  } else {
    plquants = c(.05,.95)
    cols = c('blue', 'black')

    qLow <- apply(dB$yhat.train, length(dim(dB$yhat.train)), quantile, probs = plquants[1])
    qMid <- apply(dB$yhat.train, length(dim(dB$yhat.train)), quantile, probs = 0.5)
    qUpp <- apply(dB$yhat.train, length(dim(dB$yhat.train)), quantile, probs=plquants[2])

    res <- data.frame(y = y,
                      qMid = qMid,
                      qLow = qLow,
                      qUpp = qUpp)
  }


  if (class(model) == "wbart" || class(model) == "bartMachine") {
  p <- res %>%
    tidybayes::point_interval(.fitted, y, .width = c(0.95)) %>%
    select(-y.lower, -y.upper) %>%
    ggplot() +
    tidybayes::geom_pointinterval(aes(x = y, y = .fitted, ymin = .fitted.lower, ymax = .fitted.upper),
      alpha = 0.2,
      color = "blue"
    ) +
    geom_smooth(aes(x = y, y = .fitted), method = "lm", color = "black", formula = y ~ x) +
    xlab("Actual") +
    ylab("Fitted") +
    theme_bw() +
    ggtitle("Actual vs Fitted")

  }else{
    p <-  ggplot(data = res, aes(y, qMid)) +
      geom_point(color = "blue", alpha = 0.2) +
      theme_bw() +
      xlab("Actual") +
      ylab("Fitted") +
      ggtitle("Actial vs Fitted") +
      geom_smooth(aes(x = y, y = qMid), method = "lm", color = "black", formula = y ~ x) +
      tidybayes::geom_pointinterval(aes(x = y, y = qMid, ymin = qLow, ymax = qUpp),
                                    alpha = 0.2,
                                    color = "blue"
      )
  }


  return(p)
}


# Variable Importance -----------------------------------------------------

bartVimp <- function(model) {

  if (class(model) == "wbart") {
    # get variable importance
    vImp <- model$varcount.mean
    vImp <- dplyr::tibble(Variable = names(vImp), Importance = vImp)
  } else if(class(model) == "bartMachine"){
    bmVimp <- bartMachine::investigate_var_importance(model, num_replicates_for_avg = 5)
    vimpVals <- bmVimp$avg_var_props
    vImp <- data.frame(Variable = names(vimpVals), Importance = vimpVals, row.names = NULL)
  }else{
    varCount <- NULL
    for(i in 1:length(model$varcount[,1])){
      varCount[[i]] <-  prop.table(model$varcount[i,])
    }
    vImp <- varCount %>%
      bind_rows() %>%
      colMeans()
      vImp <- tibble(Variable = names(vImp), Importance = vImp)
  }

    p <- vImp %>%
      arrange(Importance) %>%
      mutate(Variable = factor(Variable, unique(Variable))) %>%
      ggplot() +
      aes(x = Variable, y = Importance) +
      geom_segment(aes(x = Variable, xend = Variable, y = 0, yend = Importance), color = "blue") +
      geom_point(color = "blue") +
      theme_light() +
      coord_flip() +
      ggtitle(label = "BART Variable Importance") +
      theme_bw() +
      xlab("Variable") +
      ylab("Importance") +
      theme(
        axis.title.y = element_text(angle = 90, vjust = 0.5),
        legend.key.size = unit(0.5, "cm")
      )

  return(p)
}
