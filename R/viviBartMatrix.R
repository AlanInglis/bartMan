#' viviBartMatrix
#'
#'@description Returns a matrix or list of matrices. If type = 'standard' a
#' matrix filled with vivi values is returned. If type = 'vsup' two matrices are returned.
#' One with the actual values and another matrix of uncertainty values.
#' If type = 'quantiles', three matrices are returned. One for the 25%, 50%, and 75% quantiles.
#'
#'@param trees A data frame created by `extractTreeData` function.
#'@param type Which type of matrix to return. Either 'standard', 'vsup', 'quantiles'
#'@param metric Which metric to use to fill the actual values matrix. Either 'propMean' or 'count'.
#'@param metricError Which metric to use to fill the uncertainty matrix. Either 'SD', 'CV' or 'SE'.
#'@param reorder LOGICAL. If TRUE then the matrix is reordered so high values are pushed to the top left.
#'
#'
#'
#'@importFrom dplyr %>%
#'@importFrom dplyr group_by
#'@importFrom dplyr ungroup
#'@importFrom dplyr mutate
#'@importFrom dplyr select
#'@importFrom dplyr rename
#'@importFrom dplyr distinct
#'
#'@return A heatmap plot showing variable importance on the diagonal
#' and variable interaction on the off-diagonal.
#'
#' @examples
#' if(requireNamespace("dbarts", quietly = TRUE)){
#'  # Load the dbarts package to access the bart function
#'  library(dbarts)
#'  # Get Data
#'  df <- na.omit(airquality)
#'  # Create Simple dbarts Model For Regression:
#'  set.seed(1701)
#'  dbartModel <- bart(df[2:6], df[, 1], ntree = 5, keeptrees = TRUE, nskip = 10, ndpost = 10)
#'
#'  # Tree Data
#'  trees_data <- extractTreeData(model = dbartModel, data = df)
#'
#'  # VSUP Matrix
#'  vsupMat <- viviBartMatrix(trees = trees_data,
#'                            type = 'vsup',
#'                            metric = 'propMean',
#'                             metricError = 'CV')
#'  }
#'
#'@export


viviBartMatrix <- function(trees,
                           type = "standard",
                           metric = "propMean",
                           metricError = "CV",
                           reorder = FALSE
){

  if (!(metric %in% c("propMean", "count", "adjusted"))) {
    stop("metric must be \"propMean\", \"count\", or \"adjusted\"")
  }

  if (!(metricError %in% c("SD", 'SE', "CV"))) {
    stop("metricError must be \"SD\", \'SE\', or \"CV\"")
  }

  if (!(type %in% c("standard", "vsup", "quantiles"))) {
    stop("type must be \"standard\", \"vsup\", or \"quantiles\"")
  }

  viviDf <- viviBartInternal(trees)

  if(type == 'standard'){
    viviMat <-  viviBartStd(trees = trees,
                            data = viviDf,
                            metric = metric,
                            reorder = reorder
    )
  }else if(type == 'vsup'){
    viviMat <- viviBartVSUP(trees = trees,
                            data = viviDf,
                            metricError = metricError,
                            metric = metric,
                            reorder = reorder)
  }else if(type == 'quantiles'){
    viviMat <- viviBartQuantile(trees = trees,
                                data = viviDf,
                                reorder = reorder)
  }

  return(viviMat)

}




# -------------------------------------------------------------------------
# VIVI dataframe ----------------------------------------------------------
# -------------------------------------------------------------------------
viviBartInternal <- function(trees) {
  trees <- trees_data

  # Vimps -------------------------------------------------------------------

  # get vimps
  vimps <- bartMan::vimpBart(trees, type = 'prop')
  vimpsVal <- bartMan::vimpBart(trees, type = 'val')

  vImp <- colMeans(vimps)
  vimpsVal <- colSums(vimpsVal)

  n_iterations <- nrow(vimps)

  if (n_iterations > 1) {
    # get uncertainty measures
    vimpSD <- apply(vimps, 2, sd)
    upperVimp <- vImp + 1.96 * vimpSD / sqrt(n_iterations)
    lowerVimp <- vImp - 1.96 * vimpSD / sqrt(n_iterations)
    SEvimp <- sapply(as.data.frame(vimps), function(x) sd(x) / sqrt(length(x)))
    CVvimp <- vimpSD / vImp
    CVvimp[is.nan(CVvimp)] <- 0

    # get quantiles of proportions
    vimp25 <- apply(vimps, 2, function(x) quantile(x, c(.25)))
    vimp50 <- apply(vimps, 2, function(x) quantile(x, c(.50)))
    vimp75 <- apply(vimps, 2, function(x) quantile(x, c(.75)))
  } else {
    # handle single iteration case
    vimpSD <- rep(sd(vimps), length(vImp))
    upperVimp <- vImp
    lowerVimp <- vImp
    SEvimp <- rep(sd(vimps) / sqrt(length(vimps)), length(vImp))
    CVvimp <- vimpSD / vImp
    vimp25 <- vImp
    vimp50 <- vImp
    vimp75 <- vImp
  }

  # put together in dataframe
  vimpData <- cbind(vimpsVal, vImp, vimpSD, CVvimp, SEvimp,
                    lowerVimp, upperVimp, vimp25, vimp50, vimp75)

  vDf <- vimpData %>%
    as.data.frame()

  vDf$variable <- row.names(vimpData)
  vimpData <- vDf %>%
    rename(count = vimpsVal, propMean = vImp, SD = vimpSD,
           CV = CVvimp, SE = SEvimp,
           lowerCI = lowerVimp, upperCI = upperVimp,
           lowerQ = vimp25, median = vimp50, upperQ = vimp75) %>%
    select(variable, count, propMean, SD, CV, SE, lowerCI, upperCI, lowerQ, median, upperQ)

  # VInts -------------------------------------------------------------------

  df <- trees$structure
  nam <- trees$varName

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

  if (n_iterations > 1) {
    listVint <- apply(L, 2L, g)
    listVint <- listVint[lengths(listVint) > 0] # remove empty list element
  } else {
    # Handle the single iteration case
    listVint <- list()
    for (i in 1:length(L)) {
      single_iter_data <- L[[i]]
      if (length(single_iter_data) > 0) {
        single_iter_data <- single_iter_data[!is.na(single_iter_data)]
        if (length(single_iter_data) > 0) {
          listVint[[as.character(i)]] <- tapply(single_iter_data, names(single_iter_data), sum)
        }
      }
    }
  }

  nam <- trees$varName
  namDF <- expand.grid(nam, nam)

  newName <- NULL
  for (i in 1:length(namDF$Var1)) {
    newName[i] <- paste0(namDF$Var2[i], ":", namDF$Var1[i])
  }

  # turn into df
  if (length(listVint) > 0) {
    dfVint <- as.matrix(bind_rows(listVint))
    dfVint[is.na(dfVint)] <- 0
  } else {
    dfVint <- matrix(0, nrow = 1, ncol = length(newName))
    colnames(dfVint) <- newName
  }

  # create a matrix of all possible combinations
  allCombMat <- matrix(NA, nrow = max(n_iterations, nrow(dfVint)), ncol = length(newName))
  colnames(allCombMat) <- newName

  # join actual values into matrix of all combinations
  oIdx <- match(colnames(dfVint), colnames(allCombMat))

  if (nrow(dfVint) < nrow(allCombMat)) {
    missingRows <- nrow(allCombMat) - nrow(dfVint)
    dfVint <- rbind(dfVint, matrix(data = 0, ncol = ncol(dfVint), nrow = missingRows))
  }

  allCombMat[, oIdx] <- dfVint
  allCombMat[is.na(allCombMat)] <- 0
  dfVint <- allCombMat

  # reorder names to make symmetrical
  vintNames <- utils::stack(as.data.frame(dfVint))
  colnames(vintNames) = c('value', 'Var2')
  vintNames <- vintNames[, 2:1]

  dfName <- data.frame(nam = unique(vintNames$Var2))
  dfName$nam <- as.character(dfName$nam)
  newNames <- dfName %>%
    mutate(nam = map(
      strsplit(nam, ":", fixed = TRUE),
      ~ sort(.x) %>% trimws(.) %>% paste0(., collapse = ':')
    ))

  colnames(dfVint) <- newNames$nam

  # add symmetrical columns together
  dfVint <- t(apply(dfVint, 1, function(x) stats::ave(x, names(x), FUN = sum)))

  # get proportions
  propMatVint <- proportions(dfVint, 1)
  propMatVint[is.nan(propMatVint)] <- 0

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

  if (n_iterations > 1) {
    vintSD <- apply(propMatVint, 2, sd)
    vintSE <- vintSD / sqrt(n_iterations)

    vint25 <- apply(propMatVint, 2, function(x) quantile(x, c(.25)))
    vint50 <- apply(propMatVint, 2, function(x) quantile(x, c(.50)))
    vint75 <- apply(propMatVint, 2, function(x) quantile(x, c(.75)))
  } else {
    vintSD <- rep(sd(propMatVintMean), length(propMatVintMean))
    vintSE <- rep(sd(propMatVintMean) / sqrt(length(propMatVintMean)), length(propMatVintMean))
    vint25 <- propMatVintMean
    vint50 <- propMatVintMean
    vint75 <- propMatVintMean
  }

  # Ensure the vectors have names
  names(vintSD) <- names(propMatVintMean)
  names(vintSE) <- names(propMatVintMean)
  names(vint25) <- names(propMatVintMean)
  names(vint50) <- names(propMatVintMean)
  names(vint75) <- names(propMatVintMean)

  vintSD <- utils::stack(vintSD)
  colnames(vintSD) = c('SD', 'var')
  vintSD$var <- as.character(vintSD$var)

  vintSE <- utils::stack(vintSE)
  colnames(vintSE) = c('SE', 'var')
  vintSE$var <- as.character(vintSE$var)

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

  dfPropCount$SD <- errorDF$SD
  dfPropCount$SE <- errorDF$SE
  dfPropCount$Q25 <- errorDF$q25
  dfPropCount$Q50 <- errorDF$q50
  dfPropCount$Q75 <- errorDF$q75

  # Dont need to average error metrics here as it's being averaged when creating
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

  # add coeff of variance
  dfFinal$CV <- dfFinal$SD / dfFinal$propMean
  dfFinal$CV[is.infinite(dfFinal$CV)] <- 0
  dfFinal$CV[is.nan(dfFinal$CV)] <- 0

  # reorder
  dfFinal <- dfFinal %>% select(var, count, propMean, SD, CV, SE, Q25, Q50, Q75)

  # add adjustment
  vimpsAdj <- bartMan::vimpBart(trees, type = 'propMean')

  splitN <- t(simplify2array(strsplit(as.character(dfFinal[["var"]]), ":")))
  values <- t(apply(splitN, 1, function(x) { c(vimpsAdj[x[1]], vimpsAdj[x[2]]) }))
  suppressMessages(
    res <- cbind(dfFinal, values)
  )
  names(res) <- c('var', 'count', 'propMean', 'SD', 'CV',
                  'SE', 'lowerQ', 'median', 'upperQ',
                  'propVimp1', 'propVimp2')

  trans <- apply(res[, c('propMean', 'propVimp1', 'propVimp2')], 1, function(x) {
    (x[1] - x[2] * x[3]) / sqrt(x[2] * x[3])
  })

  dfFinal$adjusted <- trans
  dfFinal$adjusted[dfFinal$adjusted <= 0] <- 0
  dfFinal$adjusted[is.na(dfFinal$adjusted)] <- 0

  names(dfFinal) <- c('var', 'count', 'propMean', 'SD', "CV",
                      'SE', 'lowerQ', 'median',
                      'upperQ', 'adjusted')

  myList <- list(Vimp = vimpData, Vint = dfFinal)

  return(myList)
}





# -------------------------------------------------------------------------
# Standard Matrix ---------------------------------------------------------
# -------------------------------------------------------------------------

viviBartStd <- function(trees,
                        data,
                        reorder = FALSE,
                        metric = "propMean"
){

  propFinal <- data$Vint
  vars2  <- t(simplify2array(strsplit(as.character(propFinal[["var"]]), ":")))
  ovars <- trees$varName

  mat <- matrix(0, length(ovars), length(ovars)) # create matrix
  rownames(mat) <- colnames(mat) <- ovars # set names
  mat[vars2] <- propFinal[[metric]] # set values
  mat <- mat[lower.tri(mat, diag = T)[nrow(mat):1], ] + t(mat)[lower.tri(mat, diag = T)[nrow(mat):1], ]


  if(metric == 'adjusted'){
    vimps <- data$Vimp$propMean
  }else{
    vimps <- data$Vimp[,metric] # add vimps to matirx
  }
  diag(mat) <- vimps

  # reorder actual values matrix
  if(reorder){
    mat <- bartReorder(mat)
  }

  mat[is.nan(mat)] <- 0

  class(mat) <- c('matrix', 'array', 'standardMat')

  return(mat)

}


# -------------------------------------------------------------------------
# VSUP --------------------------------------------------------------------
# -------------------------------------------------------------------------



viviBartVSUP <- function(trees,
                         data,
                         reorder = FALSE,
                         metricError = 'CV',
                         metric = "propMean"){

  # get matrix of uncertainty values
  propFinal <- data$Vint
  vars2  <- t(simplify2array(strsplit(as.character(propFinal[["var"]]), ":")))

  ovars <- trees$varName


  mat <- matrix(0, length(ovars), length(ovars)) # create matrix
  rownames(mat) <- colnames(mat) <- ovars # set names
  mat[vars2] <- propFinal[[metricError]] # set values
  mat <- (mat[lower.tri(mat, diag = T)[nrow(mat):1], ] +
            t(mat)[lower.tri(mat, diag = T)[nrow(mat):1], ])



  vimps <- data$Vimp # add vimps to matirx
  diag(mat) <- vimps[,metricError]
  uncertaintyMatrix <- mat

  # get actual values matrix
  actualMatrix <- viviBartStd(trees = trees, data = data,
                              reorder = reorder, metric = metric)

  # reorder actual values matrix
  if(reorder){
    actualMatrix <- bartReorder(actualMatrix)
  }

  # reorder uncertainty matrix to match
  actValsColOrder <- colnames(actualMatrix)
  uncertaintyMatrix <- uncertaintyMatrix[actValsColOrder, actValsColOrder]


  # class(actualMatrix) <- c('vivid', 'matrix', 'array', 'vsup')
  # class(uncertaintyMatrix) <- c('vivid', 'matrix', 'array', 'vsup')

  myList <- list(actualMatrix = actualMatrix,
                 uncertaintyMatrix = uncertaintyMatrix)

  class(myList) <- c('list', 'vsup')

  return(myList)
}



# -------------------------------------------------------------------------
# Quantile matrix --------------------------------------------------------
# -------------------------------------------------------------------------


viviBartQuantile <-  function(trees, data, reorder = FALSE){

  propFinal <- data$Vint

  matFun <- function(trees, dataVint, data, metric){
    vars2  <- t(simplify2array(strsplit(as.character(dataVint[["var"]]), ":")))
    ovars <- trees$varName
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

  lowerQuantile  <- matFun(trees, propFinal, data,  metric = 'lowerQ')
  median <- matFun(trees, propFinal, data, metric = 'median')
  upperQuantile  <- matFun(trees, propFinal, data, metric = 'upperQ')

  # reorder median matrix
  if(reorder){
    median <- bartReorder(median)

    # reorder lower and higher matrix to match
    reorderMat <- function(targetMat, orderedMat){
      medianValsColOrder <- colnames(targetMat)
      orderedMat <- orderedMat[medianValsColOrder, medianValsColOrder]
      return(orderedMat)
    }

    lowerQuantile <- reorderMat(median, lowerQuantile)
    upperQuantile <- reorderMat(median, upperQuantile)
  }



   class(lowerQuantile) <- c('matrix', 'array', 'quant')
   class(median) <- c('matrix', 'array', 'quant')
   class(upperQuantile) <- c('matrix', 'array', 'quant')



  myList <- list(lowerQuantile = lowerQuantile,
                 median = lowerQuantile,
                 upperQuantile= upperQuantile)

  class(myList) <- c('list', 'quantiles')

  return(myList)

}


bartReorder <- function(d) {
  vImp <- diag(d)
  rvImp <- range(vImp)
  if (rvImp[2] != rvImp[1]) {
    vImp <- (vImp - rvImp[1]) / (rvImp[2] - rvImp[1])
  }
  vInt <- stats::as.dist(d)
  rvInt <- range(vInt)
  if (rvInt[2] != rvInt[1]) {
    vInt <- (vInt - rvInt[1]) / (rvInt[2] - rvInt[1])
  }
  score <- apply(as.matrix(vInt), 1, max) + vImp
  o <- DendSer::dser(-vInt, -score, cost = DendSer::costLS)
  res <- d[o, o]
  class(res)<- class(d)
  res
}








