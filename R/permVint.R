#' permVint
#'
#' @description A variable interaction evaluation which creates a null model by
#' permuting the response, rebuilding the model, and calculating the inclusion proportion (IP)
#' of adjacent splits on the null model.
#' The final result displayed is the original model's IP minus the null IP.
#'
#' @param model Model created from either the BART, dbarts or bartMachine packages.
#' @param data A data frame containing variables in the model.
#' @param treeData A data frame created by extractTreeData function.
#' @param response The name of the response for the fit.
#' @param numTreesPerm The number of trees to be used in the null model.
#' As suggested by Chipman (2009), a small number of trees is recommended (~20) to force important
#' variables to used in the model. If NULL, then the number of trees from the true model is used.
#'
#' @return A variable interaction plot.
#'
#'
#' @importFrom BART wbart
#' @importFrom dbarts bart
#' @importFrom bartMachine bartMachine
#' @importFrom bartMachine get_var_counts_over_chain
#' @importFrom dplyr tibble
#' @importFrom dplyr mutate
#' @importFrom dplyr arrange
#' @import ggplot2
#'
#' @export
#'


permVint <- function(model, data, treeData, response, numTreesPerm = NULL, plotType = 'barplot') {

  # get og model vints
  actualVint <- vints(treeData = treeData)

  # get null permutation vints
  permVint <- permBartVint(
  model = model,
  data = data,
  response = response,
  numTreesPerm = numTreesPerm
  )

  # final vints
  finalVintDF <- actualVint$dfVint
  finalVintMat <- actualVint$propMat - permVint$propMat

  # turn into df
  vintDF <- permVintDF(propData = finalVintMat,
                       dataF = finalVintDF,
                       treeData = treeData)

  # vintPlot <- permPlotFn(data = vimp,
  #                        plotType = plotType)
  vintDF$meanNew <- pmax(vintDF$propMean, 0)
  vintDF$low = pmax(vintDF$meanNew - 2 * vintDF$SE, 0)
  vintDF$high = pmax(vintDF$meanNew + 2 * vintDF$SE, 0)


  p <- vintDF %>%
    arrange(meanNew) %>%
    mutate(Variable = factor(var, unique(var))) %>%
    ggplot(aes(x = Variable, y = meanNew)) +
    ggforce::geom_link(aes(
      x = Variable, xend = Variable, yend = low,
      col = Variable, alpha = rev(stat(index))
    ),
    size = 5, n = 1000
    ) +
    ggforce::geom_link(aes(
      x = Variable, xend = Variable, yend = high,
      col = Variable, alpha = rev(stat(index))
    ),
    size = 5, n = 1000
    ) +
    geom_point(aes(x = Variable, y = meanNew), shape = 18, size = 2, color = "black") +
    coord_flip() +
    theme_bw() +
    labs(x = "Variable", y = "Importance") +
    theme(legend.position = "none")


  return(p)
}

# -------------------------------------------------------------------------

# Main function:
permBartVint <- function(model, data, response, numTreesPerm = NULL) {
  UseMethod("perBart")
}




# -------------------------------------------------------------------------

# Vint function

vints <- function(treeData){

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
  emptyRows <- unique(which(lengths(listVint)==0))
  listVint <- listVint[lengths(listVint)>0] # remove empty list element


  # turn into df
  dfVint <- as.matrix(bind_rows(listVint))
  if(length(emptyRows > 0)){
    dfVint <- rbind(dfVint, matrix(data=NA, ncol=ncol(dfVint), nrow=length(emptyRows)))

  }
  dfVint[is.na(dfVint)] <- 0


  # create a matirx of all possible combinations
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

  allCombMat[ ,oIdx] <- dfVint
  allCombMat[is.na(allCombMat)] <- 0
  dfVint <- allCombMat

  # get proportions
  propMatVint <- proportions(dfVint, 1)

  myList <- list(dfVint = dfVint, propMat = propMatVint)
  return(myList)

}



# -------------------------------------------------------------------------

# Dataframe Function

permVintDF <- function(propData, dataF, treeData){

  dfVint <- dataF
  propMatVint <- propData
  # turn into df
  dfProps <- reshape2::melt(propMatVint)
  colnames(dfProps) <- c("iteration", "var", "props")

  # add counts
  countM <- reshape2::melt(dfVint)
  colnames(countM) <- c("iteration", "var", "count")
  dfProps$count <- countM$count

  # make symmetrical interactions (ie x1:x2 == x2:x1)
  # dfFinal <- dfProps %>%
  #   mutate(
  #     var = str_extract_all(var, "\\d+"),
  #     var = map_chr(var, ~ str_glue("x{sort(.x)[[1]]}:x{sort(.x)[[2]]}"))
  #   )
  dfFinal <- dfProps %>%
    mutate(var = map(
      stringr::str_split(var, pattern = ":"),
      ~ sort(.x) %>% trimws(.) %>% paste0(., collapse = ':')
    ))

  dfFinal[is.na(dfFinal)] <- 0

  dfFinal <- dfFinal %>%
    group_by(var) %>%
    mutate(count = sum(count),
           propMean = mean(props),
           SD = sd(props),
           Q25 = quantile(props, c(.25)),
           Q50 = quantile(props, c(.50)),
           Q75 = quantile(props, c(.75)),
           SE =  sd(props)/sqrt(treeData$nMCMC)) %>%
    select(-iteration, - props) %>%
    distinct() %>%
    ungroup()


  # add adjustment
  vimpsAdj <- bartMan::vimpBart(treeData, type = 'propMean')

  splitN <- t(simplify2array(strsplit(as.character(dfFinal[["var"]]), ":")))
  values <- t(apply(splitN, 1, function(x){c(vimpsAdj[x[1]], vimpsAdj[x[2]])}))
  suppressMessages(
    res <- cbind(dfFinal, values)
  )
  names(res) <- c('var', 'count', 'propMean', 'SD',  'lowerQ', 'median',
                  'upperQ', 'SE', 'propVimp1', 'propVimp2')

  trans <- apply(res[,c('propMean', 'propVimp1', 'propVimp2')], 1, function(x){
    (x[1]-x[2]*x[3])/sqrt(x[2]*x[3])
  })

  dfFinal$adjusted <- trans

  dfFinal$adjusted[dfFinal$adjusted <=  0] <- 0
  names(dfFinal) <- c('var', 'count', 'propMean', 'SD',  'lowerQ', 'median',
                      'upperQ', 'SE', 'adjusted')

  dfFinal$var <- unlist(dfFinal$var)

  return(dfFinal)
}






# BART --------------------------------------------------------------------


perBart.wbart <- function(model, data,  response, numTreesPerm = NULL) {

  # get model info
  modelTrees <- model$treedraws$trees
  modelInfo <- unlist(strsplit(modelTrees, " "))[1:3]
  modelInfo <- gsub("(^\\d+)([\a-zA-Z0-9]*)", "\\1", modelInfo)

  nMCMC <- as.integer(modelInfo[1])
  nTree <- as.integer(modelInfo[2])
  nVar <- as.integer(modelInfo[3])
  burnIn <- length(model$sigma) - nMCMC

  # null model info
  responseIdx <- which((names(data) %in% response))
  if(is.null(numTreesPerm)){
    numTreesPerm <- nTree
  }

  # null model function
  permuteBARTFn <- function(data) {
    yPerm <- sample(data[, responseIdx], replace = FALSE)
    x <- data[, -responseIdx]

    bmodelPerm <- wbart(
      x.train = x,
      y.train = yPerm,
      nskip = burnIn,
      ndpost = nMCMC,
      nkeeptreedraws = nMCMC,
      ntree = numTreesPerm
    )

    permDF <- extractTreeData(bmodelPerm, data)
    permVints <- vints(permDF)
    return(permVints)
  }

  perMats <- permuteBARTFn(data)

  return(perMats)
}





