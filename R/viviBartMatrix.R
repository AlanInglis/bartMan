#' viviBartMatrix
#'
#' @description Returns a matrix or list of matrices. If type = 'standard' a
#' matrix filled with vivi values is returned. If type = 'vsup' two matrices are returned.
#' One with the actual values and another matrix of uncertainty values.
#' If type = 'quantiles, three matrices are returned. One for the 25%, 50%, and 75% quantiles.
#'
#'  @param treeData A data frame created by extractTreeData function.
#'  @param type Which type of matrix to return. Either 'standard', 'vsup', 'quantiles'
#'  @param metric Which metric to use to fill the actual values matrix. Either 'propMean' or 'count'.
#'  @param metricError Which metric to use to fill the uncertainty matrix. Either 'SD' or 'SEofCI'.
#'  @param reorder LOGICAL. If TRUE then the matrix is reordered so high values are pushed to the top left.
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr rename
#' @importFrom tibble rownames_to_column
#' @importFrom purrr map_chr
#' @importFrom vivid vividReorder
#'
#' @return A heatmap plot showing variable importance on the diagonal
#' and variable interaction on the off-diagonal.
#' @export


viviBartMatrix <- function(treeData, type = "standard", metric = "propMean", metricError = "SD", reorder = FALSE){

  if (!(metric %in% c("propMean", "count", "adjusted"))) {
    stop("metric must be \"propMean\", \"count\", or \"adjusted\"")
  }

  if (!(metricError %in% c("SD", "SEofCI"))) {
    stop("metricError must be \"SD\" or \"SEofCI\"")
  }

  if (!(type %in% c("standard", "vsup", "quantiles"))) {
    stop("type must be \"standard\", \"vsup\", or \"quantiles\"")
  }

  viviDf <- viviBartInternal(treeData)

  if(type == 'standard'){
   viviMat <-  viviBartStd(treeData = treeData,
                                 data = viviDf,
                                 metric = metric,
                                 reorder = reorder
                                 )
  }else if(type == 'vsup'){
    viviMat <- viviBartVSUP(treeData = treeData,
                                 data = viviDf,
                                 metricError = metricError,
                                 metric = metric)
  }else if(type == 'quantiles'){
    viviMat <- viviBartQuantile(treeData = treeData,
                                   data = viviDf,
                                   reorder = reorder)
  }

  return(viviMat)

}




# -------------------------------------------------------------------------
# VIVI dataframe ----------------------------------------------------------
# -------------------------------------------------------------------------

viviBartInternal <- function(treeData){

  # Vimps -------------------------------------------------------------------

  # get vimps
  vimps <- bartMan::vimpBart(treeData, type = 'prop')
  vimpsVal <- bartMan::vimpBart(treeData, type = 'val')
  #propVimp <- proportions(vimps, 1)
  vImp <- colMeans(vimps)
  vimpsVal <- colMeans(vimpsVal)

  # get SE
  vimpSD <- apply(vimps, 2, sd)
  upperVimp  <- vImp + 1.96 * vimpSD/sqrt(treeData$nMCMC)
  lowerVimp  <- vImp - 1.96 * vimpSD/sqrt(treeData$nMCMC)
  SEvimp <- (upperVimp - lowerVimp) / 3.92 # SE of 95% CI

  # get quantiles of proportions
  vimp25 <- apply(vimps, 2, function(x) quantile(x, c(.25)))
  vimp50 <- apply(vimps, 2, function(x) quantile(x, c(.50)))
  vimp75 <- apply(vimps, 2, function(x) quantile(x, c(.75)))

  vimpData <- cbind(vimpsVal, vImp, vimpSD, lowerVimp, upperVimp, SEvimp, vimp25, vimp50, vimp75)
  vimpData <- vimpData %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "variable") %>%
    rename(count = vimpsVal, propMean = vImp, SD = vimpSD,  lowerCI = lowerVimp, upperCI = upperVimp, SEofCI = SEvimp,
           lowerQ = vimp25, median = vimp50, upperQ = vimp75)


  # Vints -------------------------------------------------------------------

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
  dfVint <- as.matrix(bind_rows(listVint))
  dfVint[is.na(dfVint)] <- 0


  # get SE
  # vintSD <- apply(dfVint, 2, sd)
  # uncVint  <- 1.96 * vintSD/sqrt(treeData$nMCMC)
  #
  # get proportions
  propMatVint <- proportions(dfVint, 1)
  propMatVintMean <- colMeans(propMatVint)
  propMatVintMedian <- apply(propMatVint, 2, median)

  # get SE
  vintSD <- apply(dfVint, 2, sd)
  upperVint  <- propMatVintMean + 1.96 * vintSD/sqrt(treeData$nMCMC)
  lowerVint  <- propMatVintMean - 1.96 * vintSD/sqrt(treeData$nMCMC)
  SEvint <- (upperVint - lowerVint) / 3.92 # SE of 95% CI

  # get quantiles of proportions
  vint25 <- apply(propMatVint, 2, function(x) quantile(x, c(.25)))
  vint50 <- apply(propMatVint, 2, function(x) quantile(x, c(.50)))
  vint75 <- apply(propMatVint, 2, function(x) quantile(x, c(.75)))

  # turn into df
  propMM <- reshape::melt(propMatVintMean)
  propMM <- tibble::rownames_to_column(propMM, "var")

  countVint <- colSums(dfVint)
  countM <- reshape::melt(countVint)
  countM <- tibble::rownames_to_column(countM, "var")

  propMM$count <- countM$value

  # make symmetrical interactions (ie x1:x2 == x2:x1)

  ##NOTE MAYBE DELETE NOW BECAUSE TAKING CARE OF LATER
  # propMM <- propMM %>%
  #   mutate(
  #     var = str_extract_all(var, "\\d+"),
  #     var = map_chr(var, ~ str_glue("x{sort(.x)[[1]]}:x{sort(.x)[[2]]}"))
  #   )

  # add in quantile, sd, and se columns
  propMM$lowerQ <- vint25
  propMM$median <- vint50
  propMM$upperQ <- vint75
  propMM$SEofCI <- SEvint
  propMM$SD <- vintSD

  propFinal <- propMM %>%
    group_by(var) %>%
    mutate(propMean = mean(value))  %>%
    mutate(count = sum(count)) %>%
    select(var, count, propMean, SD, lowerQ, median, upperQ, SEofCI) %>%
    #distinct() %>%
    ungroup %>%
    arrange(-count)

  # add adjustment
  vimpsAdj <- bartMan::vimpBart(treeData, type = 'propMean')

  splitN <- t(simplify2array(strsplit(as.character(propFinal[["var"]]), ":")))
  values <- t(apply(splitN, 1, function(x){c(vimpsAdj[x[1]], vimpsAdj[x[2]])}))
  suppressMessages(
    res <- cbind(propFinal, values)
  )
  names(res) <- c('var', 'count', 'meanProp', 'SD',  'lowerQ', 'median',
                  'upperQ', 'SEofCI', 'propVimp1', 'propVimp2')

  trans <- apply(res[,c('meanProp', 'propVimp1', 'propVimp2')], 1, function(x){
    (x[1]-x[2]*x[3])/sqrt(x[2]*x[3])
  })

  propFinal$adjusted <- trans


  myList <- list(Vimp = vimpData, Vint = propFinal)


  return(myList)

}





# -------------------------------------------------------------------------
# Standard Matrix ---------------------------------------------------------
# -------------------------------------------------------------------------

viviBartStd <- function(treeData, data, reorder = TRUE, metric = "propMean"){

  propFinal <- data$Vint
  vars2  <- t(simplify2array(strsplit(as.character(propFinal[["var"]]), ":")))
  ovars <- treeData$varName
  mat <- matrix(0, length(ovars), length(ovars)) # create matrix
  rownames(mat) <- colnames(mat) <- ovars # set names
  mat[vars2] <- propFinal[[metric]] # set values
  mat <- mat[lower.tri(mat, diag = T)[nrow(mat):1], ] + t(mat)[lower.tri(mat, diag = T)[nrow(mat):1], ]

  if(metric == 'count'){
    mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)] # make symmetrical
  }

  if(metric == 'adjusted'){
    vimps <- data$Vimp$propMean
  }else{
    vimps <- data$Vimp[,metric] # add vimps to matirx
  }
  diag(mat) <- vimps

  # reorder actual values matrix
  if(reorder){
    mat <- vivid::vividReorder(mat)
  }

  class(mat) <- c('vivid', 'matrix', 'array', 'standardMat')

  return(mat)

}


# -------------------------------------------------------------------------
# VSUP --------------------------------------------------------------------
# -------------------------------------------------------------------------



viviBartVSUP <- function(treeData, data, reorder = TRUE, metricError = 'SD', metric = "propMean"){

  # get matrix of uncertainty values
  propFinal <- data$Vint
  vars2  <- t(simplify2array(strsplit(as.character(propFinal[["var"]]), ":")))
  ovars <- treeData$varName
  mat <- matrix(0, length(ovars), length(ovars)) # create matrix
  rownames(mat) <- colnames(mat) <- ovars # set names
  mat[vars2] <- propFinal[[metricError]] # set values
  mat <- (mat[lower.tri(mat, diag = T)[nrow(mat):1], ] +
            t(mat)[lower.tri(mat, diag = T)[nrow(mat):1], ]) / 2  # get average of symmetric values



  vimps <- data$Vimp # add vimps to matirx
  diag(mat) <- vimps[,metricError]
  uncertaintyMatrix <- mat

  # get actual values matrix
  actualMatrix <- viviBartStd(treeData, data, reorder = FALSE, metric = metric)

  # reorder actual values matrix
  if(reorder){
    actualMatrix <- vivid::vividReorder(actualMatrix)
  }

  # reorder uncertainty matrix to match
  actValsColOrder <- colnames(actualMatrix)
  uncertaintyMatrix <- uncertaintyMatrix[actValsColOrder, actValsColOrder]

  class(actualMatrix) <- c('vivid', 'matrix', 'array', 'vsup')
  class(uncertaintyMatrix) <- c('vivid', 'matrix', 'array', 'vsup')

  myList <- list(actualMatrix = actualMatrix,
                 uncertaintyMatrix = uncertaintyMatrix)

  class(myList) <- c('list', 'vsup')

  return(myList)
}


# -------------------------------------------------------------------------
# Quantile matrix --------------------------------------------------------
# -------------------------------------------------------------------------


viviBartQuantile <-  function(treeData, data, reorder = TRUE){

  propFinal <- data$Vint

  matFun <- function(treeData, dataVint, data, metric){
    vars2  <- t(simplify2array(strsplit(as.character(dataVint[["var"]]), ":")))
    ovars <- treeData$varName
    mat <- matrix(0, length(ovars), length(ovars)) # create matrix
    rownames(mat) <- colnames(mat) <- ovars # set names
    mat[vars2] <- dataVint[[metric]] # set values
    mat <- (mat[lower.tri(mat, diag = T)[nrow(mat):1], ] +
              t(mat)[lower.tri(mat, diag = T)[nrow(mat):1], ]) / 2  # get average of symmetric values

    # add vimps
    if(metric == 'adjusted'){
      vimps <- data$Vimp$propMean
    }else{
      vimps <- data$Vimp[,metric] # add vimps to matirx
    }
    diag(mat) <- vimps
    return(mat)
  }

  lowerQuantile  <- matFun(treeData, propFinal, data,  metric = 'lowerQ')
  median <- matFun(treeData, propFinal, data, metric = 'median')
  upperQuantile  <- matFun(treeData, propFinal, data, metric = 'upperQ')

  # reorder median matrix
  if(reorder){
    median <- vivid::vividReorder(median)

    # reorder lower and higher matrix to match
    reorderMat <- function(targetMat, orderedMat){
      medianValsColOrder <- colnames(targetMat)
      orderedMat <- orderedMat[medianValsColOrder, medianValsColOrder]
      return(orderedMat)
    }

    lowerQuantile <- reorderMat(median, lowerQuantile)
    upperQuantile <- reorderMat(median, upperQuantile)
  }



  class(lowerQuantile) <- c('vivid', 'matrix', 'array', 'quant')
  class(median) <- c('vivid', 'matrix', 'array', 'quant')
  class(upperQuantile) <- c('vivid', 'matrix', 'array', 'quant')



  myList <- list(lowerQuantile = lowerQuantile,
                 median = lowerQuantile,
                 upperQuantile= upperQuantile)

  class(myList) <- c('list', 'quantiles')

  return(myList)

}










