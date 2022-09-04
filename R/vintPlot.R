#' vintPlot
#'
#' @description Plot the pair-wise variable interactions inclusion porportions
#' for a BART model with the 25% and 75% quantile.
#'
#' @param treeData A data frame created by treeData function.
#' @param plotType Which type of plot to return. Either a barplot 'barplot' with the
#' quantiles shown as a line, a point plot with the quantiles shown as a gradient 'pointGrad', or a
#' letter-value plot 'lvp'.
#' @param combineFact If a variable is a factor in a data frame, when building the BART model it is replaced with dummies.
#' Note that q dummies are created if q>2 and one dummy is created if q=2, where q is the number of levels of the factor.
#' If combineFact = TRUE, then the importance is calculated for the entire factor by aggregating the dummy variables’
#' inclusion proportions.
#' @param top Display only the top X metrics (does not apply to the letter-value plot).
#'
#' @return A plot of variable importance.
#'
#' @import ggplot2
#' @importFrom dplyr tibble
#' @importFrom dplyr %>%
#' @importFrom dplyr arrange
#' @importFrom dplyr mutate
#' @importFrom purrr map
#'
#'
#' @export
#'


vintPlot <- function(treeData,
                      combineFact = FALSE,
                      plotType = 'barplot',
                      top = NULL){



  df <- treeData$structure
  nam <- treeData$varName

  # cycle through trees and create list of Vints
  mkTree <- function(x, pos = 1L) {
    var <- x[pos]
    if (is.na(var)) {
      list(NA_character_, NULL, NULL, 1L)
    } else {
      node <- vector("list", 4L)
      node[[1L]] <- var
      node[[2L]] <- l <- Recall(x, pos + 1L)
      node[[3L]] <- r <- Recall(x, pos + 1L + l[[4L]])
      node[[4L]] <- 1L + l[[4L]] + r[[4L]]
      node
    }
  }


  tabTree <- function(tree, sep = ":") {
    x <- rep.int(NA_character_, tree[[4L]])
    pos <- 1L
    recurse <- function(subtree) {
      var1 <- subtree[[1L]]
      if (!is.na(var1)) {
        for (i in 2:3) {
          var2 <- subtree[[c(i, 1L)]]
          if (!is.na(var2)) {
            x[pos] <<- paste0(var1, sep, var2)
            pos <<- pos + 1L
            Recall(subtree[[i]])
          }
        }
      }
    }
    recurse(tree)
    x <- x[!is.na(x)]
    if (length(x)) {
      x <- factor(x)
      setNames(tabulate(x), levels(x))
    } else {
      integer(0L)
    }
  }

  f <- function(x) tabTree(mkTree(x))
  L <- tapply(df[["var"]], df[c("treeNum", "iteration")], f, simplify = FALSE)

  g <- function(l) {
    x <- unlist(unname(l))
    tapply(x, names(x), sum)
  }

  listVint <- apply(L, 2L, g)
  listVint <- listVint[lengths(listVint)>0] # remove empty list element


  # turn into df
  dfVint <- as.matrix(dplyr::bind_rows(listVint))
  dfVint[is.na(dfVint)] <- 0


  # create a matrix of all possible combinations
  nam <- treeData$varName
  namDF <- expand.grid(nam, nam)

  newName <- NULL
  for(i in 1:length(namDF$Var1)){
    newName[i] <- paste0(namDF$Var2[i], ":", namDF$Var1[i])
  }

  allCombMat <- matrix(NA, nrow = treeData$nMCMC, ncol = length(newName))
  colnames(allCombMat) <- newName

  # join actual values into matirx of all combinations
  oIdx <- match(colnames(dfVint), colnames(allCombMat))

  if(nrow(dfVint) != nrow(allCombMat)){
    missingRows <-  nrow(allCombMat) - nrow(dfVint)
    dfVint <- rbind(dfVint, matrix(data = 0, ncol=ncol(dfVint), nrow=missingRows))
  }

  allCombMat[ ,oIdx] <- dfVint
  allCombMat[is.na(allCombMat)] <- 0
  dfVint <- allCombMat

  # reorder names to make symmetrical
  vintNames <- utils::stack(as.data.frame(dfVint))
  colnames(vintNames) = c('value', 'Var2')
  vintNames <- vintNames[,2:1]

  dfName <- data.frame(nam = unique(vintNames$Var2))

  dfName$nam <- as.character(dfName$nam)
  newNames <- dfName %>%
    mutate(nam = map(
      strsplit(nam, ":", fixed = T),
      ~ sort(.x) %>% trimws(.) %>% paste0(., collapse = ':')
    ))


  colnames(dfVint) <- newNames$nam

  # add symmetrical columns together
  dfVint <- t(apply(dfVint, 1, \(x) stats::ave(x, names(x), FUN = sum)))

  # get proportions
  propMatVint <- proportions(dfVint, 1)
  propMatVint[is.nan(propMatVint)] <- 0

  if(combineFact){
    propMatVint <- combineFactorsInt(propMatVint, treeData$data)
    dfVint <- combineFactorsInt(dfVint, treeData$data)
  }
  propMatVintMean <- colMeans(propMatVint)

  # turn into df
  dfProps <- utils::stack(propMatVintMean)
  colnames(dfProps) = c('props', 'var')
  dfProps$var <- as.character(dfProps$var)



  # add counts
  countMean <- colMeans(dfVint)

  # turn into df
  dfCountMean <- utils::stack(countMean)
  colnames(dfCountMean) = c('count', 'var')
  dfCountMean$var <- as.character(dfCountMean$var)

  # put together
  dfPropCount <- data.frame(
    var = dfCountMean$var,
    count = dfCountMean$count,
    props = dfProps$props

  )

  # Get uncertainty metrics -------------------------------------------------

  vintSD <- apply(propMatVint, 2, sd)

  vintSD <- utils::stack(vintSD)
  colnames(vintSD) = c('SD', 'var')
  vintSD$var <- as.character(vintSD$var)

  vintSE <- vintSD$SD/sqrt(treeData$nMCMC)
  names(vintSE) <- vintSD$var

  vintSE <- utils::stack(vintSE)
  colnames(vintSE) = c('SE', 'var')
  vintSE$var <- as.character(vintSE$var)


  # get quantiles of proportions
  vint25 <- apply(propMatVint, 2, function(x) quantile(x, c(.25)))
  vint50 <- apply(propMatVint, 2, function(x) quantile(x, c(.50)))
  vint75 <- apply(propMatVint, 2, function(x) quantile(x, c(.75)))

  vint25 <- utils::stack(vint25)
  colnames(vint25) <- c('value', 'var')
  vint25$var <- as.character(vint25$var)

  vint50 <- utils::stack(vint50)
  colnames(vint50) <- c('value', 'var')
  vint50$var <- as.character(vint50$var)

  vint75 <- utils::stack(vint75)
  colnames(vint75) <- c('value', 'var')
  vint75$var <- as.character(vint75$var)

  # put together in df
  errorDF <- data.frame(
    var = vintSD$var,
    SD = vintSD$SD,
    SE = vintSE$SE,
    q25 = vint25$value,
    q50 = vint50$value,
    q75 = vint75$value
  )

  dfPropCount$SD <-  errorDF$SD
  dfPropCount$SE <-  errorDF$SE
  dfPropCount$Q25 <- errorDF$q25
  dfPropCount$Q50 <- errorDF$q50
  dfPropCount$Q75 <- errorDF$q75

  # the matrix of values
  dfFinal <- dfPropCount %>%
    group_by(var) %>%
    mutate(count = count,
           propMean = mean(props),
           SD = mean(SD),
           SE = mean(SE),
           Q25 = mean(Q25),
           Q50 = mean(Q50),
           Q75 = mean(Q75)) %>%
    select(var, count, propMean, SD, SE, Q25, Q50, Q75, -props) %>%
    distinct() %>%
    ungroup()
  #dfFinal$var <- unlist(dfFinal$var)

  # add coeff of variance
  dfFinal$CV <- dfFinal$SD / dfFinal$propMean
  dfFinal$CV[is.nan(dfFinal$CV)] <- 0

  # reorder
  dfFinal <-  dfFinal |> select(var, count, propMean, SD, CV, SE, Q25, Q50, Q75)


  if(!is.null(top)){
    dfFinal <- dfFinal |> arrange(-propMean) |> filter(row_number() %in% 1:top)
  }



  if (plotType == "barplot") {
    p <-  dfFinal %>%
      arrange(propMean) %>%
      mutate(Variable = factor(var, unique(var))) %>%
      ggplot() +
      aes(x = Variable, y = propMean) +
      geom_bar(aes(x = Variable, y = propMean), stat = "identity", fill = "steelblue", col = "black") +
      geom_segment(aes(x = Variable, xend = Variable, y = Q25, yend = Q75), color = "black") +
      theme_light() +
      coord_flip() +
      theme_bw() +
      xlab("Variable") +
      ylab("Interaction") +
      theme(
        axis.title.y = element_text(angle = 90, vjust = 0.5),
        legend.key.size = unit(0.5, "cm")
      )
  } else if (plotType == "pointGrad") {

    if (!requireNamespace("ggforce", quietly = TRUE)) {
      stop("Package \"ggforce\" needed for this function to work. Please install it.",
           call. = FALSE)
    }

    propFinal <- transform(dfFinal, valToPlot = ifelse(Q25 == 0 & Q75 == 0, 0, propMean))

    p <- dfFinal %>%
      arrange(propMean) %>%
      mutate(Variable = factor(var, unique(var))) %>%
      ggplot(aes(x = Variable, y = propMean)) +
      ggforce::geom_link(aes(
        x = Variable, xend = Variable, yend = Q75,
        colour = "gray50", alpha = rev(stat(index))
      ),
      size = 5, n = 1000
      ) +
      ggforce::geom_link(aes(
        x = Variable, xend = Variable, yend = Q25,
        colour = 'gray50', alpha = rev(stat(index))
      ),
      size = 5, n = 1000
      ) +
      geom_point(aes(x = Variable, y = propMean), shape = 18, size = 2, color = "black") +
      coord_flip() +
      theme_bw() +
      labs(x = "Variable", y = "Importance") +
      theme(legend.position = "none")

  } else if (plotType == "lvp") {

    if (!requireNamespace("lvplot", quietly = TRUE)) {
      stop("Package \"lvplot\" needed for this function to work. Please install it.",
           call. = FALSE)
    }

    dfVint <- as.data.frame(dfVint)
    # suppressMessages(
    #   lvpVint <- reshape::melt(dfVint)
    # )

    uniNames <- unique(colnames(dfVint))
    dfVint <- dfVint[,uniNames]

    lvpVint <- utils::stack(as.data.frame(dfVint))
    colnames(lvpVint) <- c('value', 'variable')
    lvpVint$rowID <- as.numeric(lvpVint$variable)

    lvpVint <- lvpVint %>%
      arrange(-value)

    ###
    pal <- c("#08306B", "#084D96", "#1B69AF", "#3787C0", "#58A1CE",
             "#81BADA", "#ABCFE5", "#CBDEF0","#E0ECF7", "#F7FBFF")


    p <- ggplot(lvpVint, aes(stats::reorder(variable, value), value)) +
      lvplot::geom_lv(aes(fill = ..LV..),
                      conf = 0.5,
                      outlier.colour = "blue",
                      outlier.shape = 5,
                      varwidth = TRUE,
                      col = 'black'
      ) +
      scale_fill_manual(values = pal) +
      labs(x = "", y = "Importance") +
      theme_bw() +
      theme(legend.position = "none") + coord_flip()




  }

return(p)
}






