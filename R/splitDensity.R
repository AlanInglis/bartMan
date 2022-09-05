#' splitDensity
#'
#' @description Density plots of the split value for each variable.
#'
#' @param treeData A list of trees created using the treeData function.
#' @param data Data frame containing variables from the model.
#' @param display Choose how to display the plot. Either histogram, facet wrap, ridges
#' or display both the split value and density of the predictor by using dataSplit.
#'
#' @return A faceted group of density plots
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @import ggplot2
#'
#' @export


splitDensity <- function(treeData, data, display = "histogram") {

  if (!(display %in% c("histogram", "ridge", "density", 'dataSplit'))) {
    stop("display must be \"histogram\", \"ridge\", \"density\", or \"dataSplit\"")
  }

  # get just the variable and split value
  tt <- treeData$structure %>%
    ungroup() %>%
    select(var, splitValue) %>%
    stats::na.omit()

  # create plotting order to match order of data
  nam <- treeData$varName
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

    dPlot <- dfList %>%
      ggplot(aes(x = value, y = .id, fill = stat(x))) +
      ggridges::geom_density_ridges(aes(fill = .id, alpha = 0.1), panel_scaling = F) +
      facet_wrap(~variable)+
      ylab("") +
      xlab("Value") +
      theme_bw() +
      theme(legend.position = "none")
  }

  suppressMessages(print(dPlot))
  #return(dPlot)
}
