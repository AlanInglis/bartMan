#' viviBartPlot
#'
#' @description Plots a Heatmap showing variable importance on the diagonal
#' and variable interaction on the off-diagonal with uncertainty included.
#'
#' @param matrix Matrices, such as that returned by viviBartMatrix, of values to be plotted.
#' @param intPal A vector of colours to show interactions, for use with scale_fill_gradientn. Palette number has to be 2^x/2
#' @param impPal A vector of colours to show importance, for use with scale_fill_gradientn. Palette number has to be 2^x/2
#' @param intLims Specifies the fit range for the color map for interaction strength.
#' @param impLims Specifies the fit range for the color map for importance.
#' @param uncIntLims Specifies the fit range for the color map for interaction strength uncertainties.
#' @param uncImpLims Specifies the fit range for the color map for importance uncertainties.
#' @param angle The angle to rotate the x-axis labels. Defaults to zero.
#' @param border Logical. If TRUE then draw a black border around the diagonal elements.
#'
#' @import vivid
#' @importFrom ggnewscale new_scale_fill
#'
#' @return TBD
#' @export
#'
#'

viviBartPlot <- function(matrix,
                         intPal = rev(colorspace::sequential_hcl(palette = "Purples 3", n =  2^4/2)),
                         impPal = rev(colorspace::sequential_hcl(palette = "Greens 3", n =  2^4/2)),
                         intLims = NULL,
                         impLims = NULL,
                         uncIntLims = NULL,
                         uncImpLims = NULL,
                         angle = 0,
                         border = FALSE){

  p <- viviPlot(matrix = matrix,
                intPal = intPal,
                impPal = impPal,
                intLims = intLims,
                impLims = impLims,
                uncIntLims = uncIntLims,
                uncImpLims = uncImpLims,
                angle = angle,
                border = border)
  return(p)
}


# -------------------------------------------------------------------------

# Main function:
viviPlot <- function(matrix,
                     intPal = rev(colorspace::sequential_hcl(palette = "Purples 3", n =  2^4/2)),
                     impPal = rev(colorspace::sequential_hcl(palette = "Greens 3", n =  2^4/2)),
                     intLims = NULL,
                     impLims = NULL,
                     uncIntLims = NULL,
                     uncImpLims = NULL,
                     angle = 0,
                     border = FALSE) {
  UseMethod("viviPlot", matrix)
}



# -------------------------------------------------------------------------
# Standard plot -----------------------------------------------------------
# -------------------------------------------------------------------------

viviPlot.standardMat <-function(matrix,
                                intPal = rev(colorspace::sequential_hcl(palette = "Purples 3", n =  2^4/2)),
                                impPal = rev(colorspace::sequential_hcl(palette = "Greens 3", n =  2^4/2)),
                                intLims = NULL,
                                impLims = NULL,
                                uncIntLims = NULL,
                                uncImpLims = NULL,
                                angle = 0,
                                border = FALSE){

  p <- vivid::viviHeatmap(mat = matrix,
                          intPal = intPal,
                          impPal = impPal,
                          intLims = intLims,
                          impLims = impLims,
                          angle = angle,
                          border = border)
  return(p)
}


# -------------------------------------------------------------------------
# VSUP plot ---------------------------------------------------------------
# -------------------------------------------------------------------------

viviPlot.vsup <- function(matrix,
                          intPal = rev(colorspace::sequential_hcl(palette = "Purples 3", n =  2^4/2)),
                          impPal = rev(colorspace::sequential_hcl(palette = "Greens 3", n =  2^4/2)),
                          intLims = NULL,
                          impLims = NULL,
                          uncIntLims = NULL,
                          uncImpLims = NULL,
                          angle = 0,
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
  # vintBreaksLabel <- lapply(vintBreaks, function(x){
  #   round(x, 2)
  # })


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
  # vimpsBreaksLabel <- lapply(vimpsBreaks, function(x){
  #   round(x, 2)
  # })



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
      labels = vimpsBreaks,
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
        labels = vimpsBreaks,
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
      theme(axis.text.x = element_text(angle = angle, hjust = 0)) +
      theme(aspect.ratio = 1)
  )
  return(newPlt)
}



# -------------------------------------------------------------------------
# Quantile plot -----------------------------------------------------------
# -------------------------------------------------------------------------

viviPlot.quantiles <- function(matrixList,
                               intPal = rev(colorspace::sequential_hcl(palette = "Purples 3", n = 100)),
                               impPal = rev(colorspace::sequential_hcl(palette = "Greens 3", n =  100)),
                               intLims = NULL,
                               impLims = NULL,
                               uncIntLims = NULL,
                               uncImpLims = NULL,
                               angle = 0,
                               border = FALSE
){

  # get each matirx
  quant.05 <- matrixList$lowerQuantile
  quant.50 <- matrixList$median
  quant.95 <- matrixList$upperQuantile

  # Limits and Breaks ------------------------------------------------------------------

  limitsFun <- function(matrix){

    # set the limits for actual importance
    if (is.null(impLims)) {
      impLims <- range(diag(matrix))
      limitsImp <- range(labeling::rpretty(impLims[1], impLims[2]))
    } else {
      limitsImp <- impLims
    }

    # set the limits for actual interactions
    if (is.null(intLims)) {
      intLims <- range(as.dist(matrix))
      limitsInt <- range(labeling::rpretty(intLims[1], intLims[2]))
    } else {
      limitsInt <- intLims
    }

    limsList <- list(limitsImp = limitsImp, limitsInt = limitsInt)
    return(limsList)
  }


  quant.05Lim <- limitsFun(quant.05)
  quant.50Lim <- limitsFun(quant.50)
  quant.95Lim <-limitsFun(quant.95)


  # get max and min limits
  allLims <- data.frame(impLims = c(quant.05Lim$limitsImp, quant.50Lim$limitsImp, quant.95Lim$limitsImp),
                        intLims = c(quant.05Lim$limitsInt, quant.50Lim$limitsInt, quant.95Lim$limitsInt))

  lims <- list(limitsImp = c(min(allLims$impLims), max(allLims$impLims)),
               limitsInt = c(min(allLims$intLims), max(allLims$intLims))
  )



  # Create dataframe  -------------------------------------------------------

  createDataFrame <- function(matrix){

    meltedMat <- vivid:::as.data.frame.vivid(matrix)
    # get int vals
    dfInt <- meltedMat[which(meltedMat$Measure == "Vint"), ]
    # get imp vals
    dfImp <- meltedMat[which(meltedMat$Measure == "Vimp"), ]
    # get names
    nam <- colnames(matrix)
    # order factors
    dfInt$Variable_1 <- factor(dfInt$Variable_1, levels = nam)
    dfInt$Variable_2 <- factor(dfInt$Variable_2, levels = nam)

    dfList <- list(dfInt = dfInt, dfImp = dfImp)
    return(dfList)
  }

  df.05 <- createDataFrame(quant.05)
  df.50 <- createDataFrame(quant.50)
  df.95 <- createDataFrame(quant.95)


  # Create plots -----------------------------------------------------------


  plotfun <- function(dat, lims){

    dfInt <- dat$dfInt
    dfImp <- dat$dfImp
    limitsImp <- lims$limitsImp
    limitsInt <- lims$limitsInt


    p <- ggplot(dfInt, aes(.data[["Variable_1"]], .data[["Variable_2"]])) +
      geom_tile(aes(fill = .data[["Value"]])) +
      scale_x_discrete(position = "top") +
      scale_y_discrete(limits = rev(levels(dfInt$Variable_2))) +
      scale_fill_gradientn(
        colors = intPal, limits = limitsInt, name = "Vint",
        guide = guide_colorbar(
          order = 1,
          frame.colour = "black",
          ticks.colour = "black"
        ), oob = scales::squish
      ) +
      new_scale_fill() +
      geom_tile(data = dfImp,
                aes(fill = .data[["Value"]])
      ) +
      scale_fill_gradientn(
        colors = impPal, limits = limitsImp, name = "Vimp",
        guide = guide_colorbar(
          order = 2,
          frame.colour = "black",
          ticks.colour = "black"
        ), oob = scales::squish
      ) +
      xlab("") +
      ylab("") +
      theme_light() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      theme(axis.text = element_text(size = 11)) +
      theme(axis.text.x = element_text(angle = angle, hjust = 0)) +
      theme(aspect.ratio = 1)

    if(border){
      p$layers[[2]]$aes_params$colour = 'black'
      p$layers[[2]]$aes_params$size = 0.2
    }

    return(p)
  }

  p1 <- plotfun(df.05, lims = lims)
  p2 <- plotfun(df.50, lims = lims)
  p3 <- plotfun(df.95, lims = lims) + theme(legend.position = "bottom")
  #theme(legend.key.size = unit(0.5, "cm"))

  legendFinal <- cowplot::get_legend(p3)

  p1 <- p1 + ggtitle("5% quantile") + theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))
  p2 <- p2 + ggtitle("Median") + theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))
  p3 <- p3 + ggtitle("95% quantile") + theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))


  design <- c(
    area(1,1),
    area(1,2),
    area(1,3),
    area(2,2)
  )

  allPlots <- cowplot::plot_grid(p1,p2,p3,
                                 NULL, legendFinal, ncol = 3, nrow = 2,
                                 rel_heights = c(1.5,0.5))

  return(allPlots)
}

