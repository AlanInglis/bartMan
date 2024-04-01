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
#' @param unc_levels The number of uncertainty levels
#' @param max_desat The maximum desaturation level.
#' @param pow_desat The power of desaturation level.
#' @param max_light The maximum light level.
#' @param pow_light The power of light level.
#' @param label legend label for the uncertainty measure.
#'
#' @importFrom ggnewscale new_scale_fill
#' @importFrom stats as.dist
#'
#' @return Either a heatmap, VSUP, or quantile heatmap plot.
#'
#' @examples
#' \dontrun{
#' df_trees <- extractTreeData(model = my_model, data = my_data)
#' vsupMat <- viviBartMatrix(df_trees, type = 'vsup', metric = 'propMean', metricError = "CV")
#' viviBartPlot(vsupMat, label = 'CV')
#' }
#'
#'
#' @export


viviBartPlot <- function(matrix,
                         intPal = NULL,
                         impPal = NULL,
                         intLims = NULL,
                         impLims = NULL,
                         uncIntLims = NULL,
                         uncImpLims = NULL,
                         unc_levels = 4,
                         max_desat = 0.6,
                         pow_desat = 0.2,
                         max_light = 0.6,
                         pow_light = 1,
                         angle = 0,
                         border = FALSE,
                         label = NULL){

  if(is.null(intPal)){
    intPal <- scales::colour_ramp(
      colors = c(blue = '#FFFFCC', red = '#800026')
    )((0:7)/7)
  }
  if(is.null(impPal)){
    impPal <-  c("#E0F3DB", "#CCEBC5", "#A8DDB5", "#7BCCC4",
                 "#4EB3D3", "#2B8CBE", "#0868AC", "#084081")

  }

  p <- viviPlot(matrix = matrix,
                intPal = intPal,
                impPal = impPal,
                intLims = intLims,
                impLims = impLims,
                uncIntLims = uncIntLims,
                uncImpLims = uncImpLims,
                unc_levels = unc_levels,
                max_desat = max_desat,
                pow_desat = pow_desat,
                max_light = max_light,
                pow_light = pow_light,
                angle = angle,
                border = border,
                label = label)
  return(p)
}


# -------------------------------------------------------------------------

# Main function:
viviPlot <- function(matrix,
                     intPal = NULL,
                     impPal = NULL,
                     intLims = NULL,
                     impLims = NULL,
                     uncIntLims = NULL,
                     uncImpLims = NULL,
                     unc_levels = 4,
                     max_desat = 0.6,
                     pow_desat = 0.2,
                     max_light = 0.6,
                     pow_light = 1,
                     angle = 0,
                     border = FALSE,
                     label = NULL) {
  UseMethod("viviPlot", matrix)
}



# -------------------------------------------------------------------------
# Standard plot -----------------------------------------------------------
# -------------------------------------------------------------------------

viviPlot.standardMat <-function(matrix,
                                intPal = NULL,
                                impPal = NULL,
                                intLims = NULL,
                                impLims = NULL,
                                uncIntLims = NULL,
                                uncImpLims = NULL,
                                angle = 0,
                                border = FALSE,
                                unc_levels = 4,
                                max_desat = 0.6,
                                pow_desat = 0.2,
                                max_light = 0.6,
                                pow_light = 1,
                                label = NULL){

  p <- bartHeatmap(mat = matrix,
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
                          unc_levels = 4,
                          max_desat = 0.6,
                          pow_desat = 0.2,
                          max_light = 0.6,
                          pow_light = 1,
                          angle = 0,
                          border = FALSE,
                          label = NULL
){

  if(is.null(label)){
    label <- 'Uncertainty'
  }

  # get values
  actualMatrix <- matrix$actualMatrix
  uncertMatrix <- matrix$uncertaintyMatrix


  # Limits and Breaks ------------------------------------------------------------------

  # set the limits for actual importance
  if (is.null(impLims)) {
    impLims <- range(diag(actualMatrix))
    limitsImp <- range(pretty(c(impLims[1], impLims[2])))#range(labeling::rpretty(impLims[1], impLims[2]))
  } else {
    limitsImp <- impLims
  }

  # set the limits for actual interactions
  if (is.null(intLims)) {
    intLims <- range(stats::as.dist(actualMatrix))
    limitsInt <-  range(pretty(c(intLims[1], intLims[2])))#range(labeling::rpretty(intLims[1], intLims[2]))
  } else {
    limitsInt <- intLims
  }


  # set the limits for uncert importance
  if (is.null(uncImpLims)) {
    uncImpLims <- range(diag(uncertMatrix))
    limitsImpUnc <-  range(pretty(c(uncImpLims[1], uncImpLims[2])))#range(labeling::rpretty(uncImpLims[1], uncImpLims[2]))
  } else {
    limitsImpUnc <- uncImpLims
  }

  # set the limits for uncert interactions
  if (is.null(uncIntLims)) {
    uncIntLims <- range(stats::as.dist(uncertMatrix))
    limitsIntUnc <-  range(pretty(c(uncIntLims[1], uncIntLims[2])))#range(labeling::rpretty(uncIntLims[1], uncIntLims[2]))
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

  # vintBreaks <- lapply(vintBreaks, function(x){
  #   quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1))
  # }
  # )
  # vintBreaks <- lapply(vintBreaks, function(x){
  #   unname(x)
  # })

  vintBreaks[[1]] <- seq(vintBreaks[[1]][1], vintBreaks[[1]][2], length.out = 5)
  vintBreaks[[2]] <- seq(vintBreaks[[2]][1], vintBreaks[[2]][2], length.out = 5)


  vintBreaksLabel <- vintBreaks
  vintBreaksLabel[[1]] <- round(vintBreaksLabel[[1]], 3)
  vintBreaksLabel[[2]] <- round(vintBreaksLabel[[2]], 5)
  vintBreaksLabel[[2]] <-  rev(vintBreaksLabel[[2]])

  # vintBreaksLabel <- lapply(vintBreaks, function(x){
  #   round(x, 4)
  # })


  vimpsBreaks <- list(c(limitsImp), c(limitsImpUnc))
  vimpLims <- vimpsBreaks
  vimpLims[[1]][1] <- vimpLims[[1]][1] - 0.001
  vimpLims[[1]][2] <- vimpLims[[1]][2] + 0.001
  vimpLims[[2]][1] <- vimpLims[[2]][1] - 0.001
  vimpLims[[2]][2] <- vimpLims[[2]][2] + 0.001

  # vimpsBreaks <- lapply(vimpsBreaks, function(x){
  #   quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1))
  # })
  # vimpsBreaks <- lapply(vimpsBreaks, function(x){
  #   unname(x)
  # })

  vimpsBreaks[[1]] <- seq(vimpsBreaks[[1]][1], vimpsBreaks[[1]][2], length.out = 5)
  vimpsBreaks[[2]] <- seq(vimpsBreaks[[2]][1], vimpsBreaks[[2]][2], length.out = 5)

  vimpBreaksLabel <- vimpsBreaks
  vimpBreaksLabel[[1]] <- round(vimpBreaksLabel[[1]], 3)
  vimpBreaksLabel[[2]] <- round(vimpBreaksLabel[[2]], 5)
  vimpBreaksLabel[[2]] <-rev( vimpBreaksLabel[[2]])
  # vimpBreaksLabel <- lapply(vimpsBreaks, function(x){
  #   round(x, 4)
  # })



  # Create dataframe  -------------------------------------------------------

  # turn into dataframe for plotting
  meltedMat <- as.data.frame.bartMan(actualMatrix)
  meltedUnc <- as.data.frame.bartMan(uncertMatrix)

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

  # label name
  # if(is.null(labelName)){
  #   labelName <- ()
  # }

  # create plot for Vint ----------------------------------------------------

  pInt <- ggplot(dfInt) +
    geom_raster(aes(x = Variable_1, y = Variable_2, fill = zip(Value, Uncert))) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(limits = rev(levels(dfInt$Variable_2))) +
    coord_equal() +
    bivariate_scale(
      name = c("Vint", label),
      aesthetics = "fill",
      limits = vintLims,
      breaks = vintBreaks,
      labels = vintBreaksLabel,
      oob = scales::squish,
      palette = pal_vsup(
        values = intPal,
        unc_levels = unc_levels,
        max_desat = max_desat,
        pow_desat = pow_desat,
        max_light = max_light,
        pow_light = pow_light
      ),
      guide = "colorfan"
    ) +
    theme_bw()

  # create plot for Vimp ----------------------------------------------------

  suppressMessages(
    newPlt <- pInt +
      new_scale_fill() +
      geom_raster(data = dfImp, aes(x = Variable_1, y = Variable_2, fill = zip(Value, Uncert))) +
      coord_equal() +
      bivariate_scale(
        name = c("Vimp", label),
        aesthetics = "fill",
        limits = vimpLims,
        breaks = vimpsBreaks,
        labels = vimpBreaksLabel,
        oob = scales::squish,
        palette = pal_vsup(
          values = impPal,
          unc_levels = unc_levels,
          max_desat = max_desat,
          pow_desat = pow_desat,
          max_light = max_light,
          pow_light = pow_light
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
      theme(axis.text.x = element_text(angle = angle)) +
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
                               border = FALSE,
                               unc_levels = 4,
                               max_desat = 0.6,
                               pow_desat = 0.2,
                               max_light = 0.6,
                               pow_light = 1,
                               label = NULL
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
      limitsImp <-  range(pretty(c(impLims[1], impLims[2])))#range(labeling::rpretty(impLims[1], impLims[2]))
    } else {
      limitsImp <- impLims
    }

    # set the limits for actual interactions
    if (is.null(intLims)) {
      intLims <- range(stats::as.dist(matrix))
      limitsInt <-  range(pretty(c(intLims[1], intLims[2])))#range(labeling::rpretty(intLims[1], intLims[2]))
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

    meltedMat <- as.data.frame.bartMan(matrix)
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

  if(angle > 10){
    hj <- 0
  }else{
    hj <- 0.5
  }


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
      theme(axis.text.x = element_text(angle = angle, hjust = hj)) +
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


as.data.frame.bartMan <- function(x, row.names = NULL, optional = FALSE, ...) {

  matrix <- x
  df <- cbind(expand.grid(dimnames(matrix), stringsAsFactors = FALSE), value = as.vector(matrix) )

  # get the row and column index
  Row <- as.vector(row(matrix))
  Col <- as.vector(col(matrix))

  # Create measure column
  df$Measure <- with(df, ifelse(Var1 == Var2, "Vimp", "Vint"))

  # cbind them together
  viviDataFrame <- cbind(df, Row, Col)

  # set names
  names(viviDataFrame)[1] <- "Variable_1"
  names(viviDataFrame)[2] <- "Variable_2"
  names(viviDataFrame)[3] <- "Value"

  return(viviDataFrame)
}


bartHeatmap <- function(mat,
                        intPal = rev(colorspace::sequential_hcl(palette = "Purples 3", n = 100)),
                        impPal = rev(colorspace::sequential_hcl(palette = "Greens 3", n = 100)),
                        intLims = NULL,
                        impLims = NULL,
                        border = FALSE,
                        angle = 0) {




  # Small set-up ------------------------------------------------------------

  # get label names
  labelNames <- colnames(mat)

  # Limits ------------------------------------------------------------------

  # set the limits for importance
  if (is.null(impLims)) {
    impLims <- range(diag(mat))
    limitsImp <- range(pretty(c(impLims[1], impLims[2])))#range(labeling::rpretty(impLims[1], impLims[2]))
  } else {
    limitsImp <- impLims
  }

  # set the limits for interactions
  if (is.null(intLims)) {
    intLims <- range(stats::as.dist(mat))
    limitsInt <- range(pretty(c(intLims[1], intLims[2])))#range(labeling::rpretty(intLims[1], intLims[2]))
  } else {
    limitsInt <- intLims
  }



  # Set up plot -------------------------------------------------------

  df <- as.data.frame.bartMan(mat)
  # get int vals
  dfInt <- df[which(df$Measure == "Vint"), ]
  # get imp vals
  dfImp <- df[which(df$Measure == "Vimp"), ]


  # Create Plot ------------------------------------------------------------

  # order factors
  dfInt$Variable_1 <- factor(dfInt$Variable_1, levels = labelNames)
  dfInt$Variable_2 <- factor(dfInt$Variable_2, levels = labelNames)

  if(angle > 10){
    hj <- 0
  }else{
    hj <- 0.5
  }


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
    theme(axis.text.x = element_text(angle = angle, hjust = hj)) +
    theme(aspect.ratio = 1)

  if(border){
    p$layers[[2]]$aes_params$colour = 'black'
    p$layers[[2]]$aes_params$size = 0.2
  }


  return(p)
}




