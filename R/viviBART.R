#' viviBART
#'
#' @description Creates a matrix displaying variable importance on the diagonal
#'  and variable interaction on the off-diagonal with the uncertainty included.
#'
#' @param fit A supervised machine learning model, which understands condvis2::CVpredict
#' @param data Data frame used for fit.
#' @param response The name of the response for the fit.
#' @param noReplications How many time to repeat the calculation.
#' @param gridSize The size of the grid for evaluating the predictions.
#' @param nmax Maximum number of data rows to consider. Default is 500. Use all rows if NULL.
#' @param reorder If TRUE (default) uses DendSer to reorder the matrix of interactions and variable importances.
#' @param class Category for classification, a factor level, or a number indicating which factor level.
#' @param normalized Should Friedman's H-statistic be normalized or not. Default is FALSE.
#' @return A list of matrices. One matrix is of interaction values with importance on the diagonal. The second
#' is a matrix of their uncertainties.
#'
#' @importFrom vivid vivi
#' @importFrom vivid vividReorder
#' @importFrom dplyr bind_rows
#' @importFrom stats xtabs
#' @import ggplot2
#'
#' @examples
#'
#' aq <- na.omit(airquality)
#' f <- lm(Ozone ~ ., data = aq)
#' m <- vivi(fit = f, data = aq, response = "Ozone") # as expected all interactions are zero
#' viviHeatmap(m)
#'
#'\donttest{
#' library(ranger)
#' rf <- ranger(Species ~ ., data = iris, importance = "impurity", probability = TRUE)
#' vivi(fit = rf, data = iris, response = "Species")
#' }
#' @export

# -------------------------------------------------------------------------

viviBART <- function(model,
                     response,
                     data,
                     noReplications = 2,
                     gridSize = 50,
                     nmax = 500,
                     reorder = TRUE,
                     class = 1,
                     normalized = FALSE){

  if(noReplications < 2){
    stop("noReplications must be > 1")
  }

  classif <- is.factor(data[[response]])
  if(classif){
    message("Importance works for numeric and numeric binary response only; setting importance to 1.")
  }

  # if(class(model) == "pbart"){
  #  pFun <- function(fit, data, prob=TRUE) as.numeric(condvis2::CVpredict(fit, data[,-response]))
  # } else if(SOMETHING){
  #   pFun <- function(fit, data, prob=TRUE) apply(predict(fit, data[,-response]), 2, mean)
  # } else {
  #   pFun <- NULL
  # }


  # need to fix... need to make response actual response
  if(class(model) == "pbart"){
    pFun <- function(fit, data, prob=TRUE) as.numeric(condvis2::CVpredict(fit, data[,-response]))
  } else {
    pFun <- NULL
  }

  mat <- NULL
  vimp = NULL

  # calc uncertainty for the variable importance:
  for(i in 1:noReplications){
    message("Replication ", i, "...")
    mat[[i]] <- vivid::vivi(fit = model,
                            data = data,
                            response = response,
                            reorder = F,
                            gridSize = gridSize,
                            nmax = nmax,
                            normalized = normalized,
                            class = class,
                            predictFun = pFun)
    vimp[[i]] <- diag(mat[[i]])
  }
  vimps <- bind_rows(vimp)
  apSD <- apply(vimps, 2, sd)
  uncVimp  <- 1.96 * apSD/sqrt(noReplications)

  avgVimp <- apply(vimps, 2, mean) # average importance

  # get uncertainty for interactions
  getIntValues <- function(matrix){
    matrix[lower.tri(matrix)] <- 0
    diag(matrix) <- 0
    dfInt <- as.data.frame(matrix)
    dfInt <- dfInt[dfInt$Value != 0,] %>%
      mutate(name = paste(Variable_1, Variable_2)) %>%
      select(name, Value)
    vintAll <- setNames(dfInt$Value, dfInt$name)
    return(vintAll)
  }

  vInt <- lapply(mat, getIntValues) %>% bind_rows()
  apSDvint <- apply(vInt, 2, sd)
  uncVint <- 1.96 * apSDvint/sqrt(noReplications)

  avgVint <- apply(vInt, 2, mean) # average interaction


  # Turn back into matrix ---------------------------------------------------

  # 1st uncertainty matrix
  uncMat <- read.table(text = names(uncVint))
  names(uncMat) <- c("Var1", "Var2")
  uncMat <- xtabs(uncVint ~ ., cbind(rbind(uncMat, setNames(rev(uncMat), names(uncMat))), uncVint = rep(uncVint, 2)))
  diag(uncMat) <- uncVimp

  # 2nd actual matrix
  actMat <- read.table(text = names(avgVint))
  names(actMat) <- c("Var1", "Var2")
  actMat <- xtabs(avgVint ~ ., cbind(rbind(actMat, setNames(rev(actMat), names(actMat))), avgVint = rep(avgVint, 2)))
  diag(actMat) <- avgVimp

  # reorder actual values matrix
  if(reorder){
  actualVals <- vivid::vividReorder(actMat)
  }else{
    actualVals <- actMat
  }

  # reorder uncertainty matrix to match
  actValsColOrder <- colnames(actualVals)
  uncMat <- uncMat[actValsColOrder, actValsColOrder]

  class(actualVals) <- c('vivid', 'matrix', 'array')
  class(uncMat) <- c('vivid', 'matrix', 'array')

  myList <- list(actualMatrix = actualVals,
                 uncertaintyMatrix = uncMat)

  return(myList)

}


# palette number has to be 2^x
#' viviBartPlot
#'
#' @description Plots a Heatmap showing variable importance on the diagonal
#' and variable interaction on the off-diagonal with uncertainty included.
#'
#' @param mat Matrices, such as that returned by viviBART, of values to be plotted.
#' @param intPal A vector of colours to show interactions, for use with scale_fill_gradientn.
#' @param impPal A vector of colours to show importance, for use with scale_fill_gradientn.
#' @param intLims Specifies the fit range for the color map for interaction strength.
#' @param impLims Specifies the fit range for the color map for importance.
#' @param border Logical. If TRUE then draw a black border around the diagonal elements.
#' @param angle The angle to rotate the x-axis labels. Defaults to zero.
#'
#' @import ggplot2
#' @importFrom ggnewscale new_scale_fill
#' @importFrom stats as.dist
#'
#' @return A heatmap plot showing variable importance on the diagonal
#' and variable interaction on the off-diagonal.
#' @export

viviBartPlot <- function(matrix,
                         intPal = rev(colorspace::sequential_hcl(palette = "Purples 3", n =  2^4/2)),
                         impPal = rev(colorspace::sequential_hcl(palette = "Greens 3", n =  2^4/2)),
                         intLims = NULL,
                         impLims = NULL,
                         uncIntLims = NULL,
                         uncImpLims = NULL,
                         border = FALSE
                          ){

  # get values
  actualMatrix <- matrix$actualMatrix
  uncertMatrix <- matrix$uncertaintyMatrix


  # Limits and Breaks ------------------------------------------------------------------

  # set the limits for actual importance
  if (is.null(impLims)) {
    impLims <- range(diag(actualMatrix))
    limitsImp <- range(labeling::rpretty(impLims[1], impLims[2]))
  } else {
    limitsImp <- impLims
  }

  # set the limits for actual interactions
  if (is.null(intLims)) {
    intLims <- range(as.dist(actualMatrix))
    limitsInt <- range(labeling::rpretty(intLims[1], intLims[2]))
  } else {
    limitsInt <- intLims
  }


  # set the limits for uncert importance
  if (is.null(uncImpLims)) {
    uncImpLims <- range(diag(uncertMatrix))
    limitsImpUnc <- range(labeling::rpretty(uncImpLims[1], uncImpLims[2]))
  } else {
    limitsImpUnc <- uncImpLims
  }

  # set the limits for uncert interactions
  if (is.null(uncIntLims)) {
    uncIntLims <- range(as.dist(uncertMatrix))
    limitsIntUnc <- range(labeling::rpretty(uncIntLims[1], uncIntLims[2]))
  } else {
    limitsIntUnc <- uncIntLims
  }

  # making sure the breaks are inside the limits
  vintBreaks <- list(c(limitsInt), c(limitsIntUnc))
  vintLims <- vintBreaks
  vintLims[[1]][1] <- vintLims[[1]][1] - 0.001
  vintLims[[1]][2] <- vintLims[[1]][2] + 0.001
  vintLims[[2]][1] <- vintLims[[2]][1] - 0.001
  vintLims[[2]][2] <- vintLims[[2]][2] + 0.001

  vintBreaks <- lapply(vintBreaks, function(x){
    quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1))
    }
  )
  vintBreaks <- lapply(vintBreaks, function(x){
    unname(x)
  })
  vintBreaksLabel <- lapply(vintBreaks, function(x){
    round(x, 2)
  })


  vimpsBreaks <- list(c(limitsImp), c(limitsImpUnc))
  vimpLims <- vimpsBreaks
  vimpLims[[1]][1] <- vimpLims[[1]][1] - 0.001
  vimpLims[[1]][2] <- vimpLims[[1]][2] + 0.001
  vimpLims[[2]][1] <- vimpLims[[2]][1] - 0.001
  vimpLims[[2]][2] <- vimpLims[[2]][2] + 0.001

  vimpsBreaks <- lapply(vimpsBreaks, function(x){
    quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1))
    })
  vimpsBreaks <- lapply(vimpsBreaks, function(x){
    unname(x)
  })
  vimpsBreaksLabel <- lapply(vimpsBreaks, function(x){
    round(x, 2)
  })



  # Create dataframe  -------------------------------------------------------

  # turn into dataframe for plotting
  meltedMat <- vivid:::as.data.frame.vivid(actualMatrix)
  meltedUnc <- vivid:::as.data.frame.vivid(uncertMatrix)

  # add uncertainty to actual dataframe
  meltedMat$Uncert <- meltedUnc$Value

  # get actual int vals
  dfInt <- meltedMat[which(meltedMat$Measure == "Vint"), ]
  # get actual imp vals
  dfImp <- meltedMat[which(meltedMat$Measure == "Vimp"), ]

  # get names
  nam <- colnames(actualMatrix)
  # order factors
  dfInt$Variable_1 <- factor(dfInt$Variable_1, levels = nam)
  dfInt$Variable_2 <- factor(dfInt$Variable_2, levels = nam)



  # create plot for Vint ----------------------------------------------------

  pInt <- ggplot(dfInt) +
    geom_raster(aes(x = Variable_1, y = Variable_2, fill = zip(Value, Uncert))) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(limits = rev(levels(dfInt$Variable_2))) +
    coord_equal() +
    bivariate_scale(
      name = c("Vint", "Uncertainty"),
      aesthetics = "fill",
      limits = vintLims,
      breaks = vintBreaks,
      labels = vintBreaksLabel,
      oob = scales::squish,
      palette = pal_vsup(
        values = intPal,
        unc_levels = 4,
        max_desat = 0.6,
        pow_desat = 0.2,
        max_light = 0.6,
        pow_light = 1
      ),
      guide = "colorfan"
    ) +
    theme_bw() +
    theme(
      legend.title.align = 0.5,
      legend.key.size = grid::unit(0.8, "cm"),
      plot.margin = ggplot2::margin(5.5, 20, 5.5, 5.5)
    )

# create plot for Vimp ----------------------------------------------------

  suppressMessages(
  newPlt <- pInt +
    new_scale_fill() +
    geom_raster(data = dfImp, aes(x = Variable_1, y = Variable_2, fill = zip(Value, Uncert))) +
    coord_equal() +
    bivariate_scale(
      name = c("Vimp", "Uncertainty"),
      aesthetics = "fill",
      limits = vimpLims,
      breaks = vimpsBreaks,
      labels = vimpsBreaksLabel,
      oob = scales::squish,
      palette = pal_vsup(
        values = impPal,
        unc_levels = 4,
        max_desat = 0.6,
        pow_desat = 0.2,
        max_light = 0.6,
        pow_light = 1
      ),
      guide = "colorfan"
    ) +
    xlab("") +
    ylab("") +
    theme_light() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    theme(axis.text = element_text(size = 11)) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0)) +
    theme(aspect.ratio = 1)
  )
  return(newPlt)
}


