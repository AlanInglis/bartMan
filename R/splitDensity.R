#' splitDensity
#'
#' @description Density plots of the split value for each variable.
#'
#' @param trees A list of trees created using the trees function.
#' @param data Data frame containing variables from the model.
#' @param bandWidth Bandwidth used for density calculation. If not provided, is estimated from the data.
#' @param panelScale If TRUE, the default, relative scaling is calculated separately for each panel.
#'  If FALSE, relative scaling is calculated globally.
#'  @param scaleFactor A scaling factor to scale the height of the ridgelines relative to the spacing between them.
#'   A value of 1 indicates that the maximum point of any ridgeline touches the baseline right above,
#'   assuming even spacing between baselines.
#' @param display Choose how to display the plot. Either histogram, facet wrap, ridges
#' or display both the split value and density of the predictor by using dataSplit.
#' @param scaleFactor A numerical value to scale the plot.
#'
#' @return A faceted group of density plots
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @import ggplot2
#' @examples
#' \donttest{
#' df_trees <- extractTreeData(model = my_model, data = my_data)
#' splitDensity(trees = df_trees, data = my_data, display = 'dataSplit')
#' }
#'
#' @export


splitDensity <- function(trees,
                         data,
                         bandWidth = NULL,
                         panelScale = NULL,
                         scaleFactor = NULL,
                         display = "histogram") {

  if (!(display %in% c("histogram", "ridge", "density", 'dataSplit'))) {
    stop("display must be \"histogram\", \"ridge\", \"density\", or \"dataSplit\"")
  }

  # get just the variable and split value
  tt <- trees$structure %>%
    ungroup() %>%
    select(var, splitValue) %>%
    stats::na.omit()

  # create plotting order to match order of data
  nam <- trees$varName
  tt$var <- factor(tt$var, levels = nam)

  varNames <- unique(tt$var)

  # create plot

  if (display == "density") {
    dPlot <- tt %>%
      ggplot(aes(x = splitValue)) +
      geom_density(aes(colour = var, fill = var)) +
      facet_wrap(~var) +
      theme_bw() +
      ylab("Density") +
      xlab("Split value") +
      theme(legend.position = "none")
  } else if (display == "ridge") {
    dPlot <- tt %>%
      ggplot(aes(x = splitValue, y = var, fill = stat(x))) +
      ggridges::geom_density_ridges(aes(fill = var, alpha = 0.1)) +
      ylab("Variable") +
      xlab("Split value") +
      theme_bw() +
      theme(legend.position = "none")
  } else if(display == "histogram") {
    dPlot <- tt %>%
      ggplot(aes(x = splitValue)) +
      geom_histogram(aes(colour = var, fill = var), bins = 30) +
      facet_wrap(~var) +
      theme_bw() +
      ylab("Density") +
      xlab("Split value") +
      theme(legend.position = "none")
  }else if(display == 'dataSplit'){
    if (!requireNamespace("ggridges", quietly = TRUE)) {
      stop("Package \"ggridges\" needed for this function to work. Please install it.",
           call. = FALSE)
    }

    dataIdx <- which((names(data) %in% varNames))
    dat <- data[, dataIdx]

    meltDat <- utils::stack(dat)
    colnames(meltDat) <- c('value', 'variable')
    names(tt) <- c('variable', 'value')

    dataList <- list(meltDat, tt)
    names(dataList)  <- c('data', 'split_value')
    #dfList <- plyr::ldply(dataList)

    dfList <- rbind(dataList$data, dataList$split_value)
    dfList$.id <- c(rep('data', length(dataList$data$value)),
                    rep('split  \nvalue', length(dataList$split_value$value)))

    dfList <- dfList |> select(.id, value, variable)

    if(is.null(bandWidth)){
      bandWidth = 0.2
    }else{
      bandWidth = bandWidth
    }

    if(is.null(scaleFactor)){
      scaleFactor = 0.5
    }else{
      scaleFactor = scaleFactor
    }

    dPlot <- dfList %>%
      ggplot(aes(x = value, y = .id, fill = stat(x))) +
      ggridges::geom_density_ridges(bandwidth = bandWidth,
                                    scale = scaleFactor,
                                    aes(fill = .id,
                                        alpha = 0.1),
                                    panel_scaling = panelScale) +
      facet_wrap(~variable)+
      ylab("") +
      xlab("Value") +
      theme_bw() +
      theme(legend.position = "none")
  }

  suppressMessages(print(dPlot))
  #return(dPlot)
}
