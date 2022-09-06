#' viviBartMatrix
#'
#'@description Returns a matrix or list of matrices. If type = 'standard' a
#' matrix filled with vivi values is returned. If type = 'vsup' two matrices are returned.
#' One with the actual values and another matrix of uncertainty values.
#' If type = 'quantiles', three matrices are returned. One for the 25%, 50%, and 75% quantiles.
#'
#'@param treeData A data frame created by extractTreeData function.
#'@param type Which type of matrix to return. Either 'standard', 'vsup', 'quantiles'
#'@param metric Which metric to use to fill the actual values matrix. Either 'propMean' or 'count'.
#'@param metricError Which metric to use to fill the uncertainty matrix. Either 'SD', 'CV' or 'SE'.
#'@param reorder LOGICAL. If TRUE then the matrix is reordered so high values are pushed to the top left.
#'@param combineFact If a variable is a factor in a data frame, when building the BART model it is replaced with dummies.
#' Note that q dummies are created if q>2 and one dummy is created if q=2, where q is the number of levels of the factor.
#' If combineFact = TRUE, then both the importance and interactions are calculated for the entire factor by aggregating the dummy variables’
#' inclusion proportions.
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
#'@export


viviBartMatrix <- function(treeData,
                           type = "standard",
                           metric = "propMean",
                           metricError = "CV",
                           reorder = FALSE,
                           combineFact = FALSE){

  if (!(metric %in% c("propMean", "count", "adjusted"))) {
    stop("metric must be \"propMean\", \"count\", or \"adjusted\"")
  }

  if (!(metricError %in% c("SD", 'SE', "CV"))) {
    stop("metricError must be \"SD\", \'SE\', or \"CV\"")
  }

  if (!(type %in% c("standard", "vsup", "quantiles"))) {
    stop("type must be \"standard\", \"vsup\", or \"quantiles\"")
  }

  viviDf <- viviBartInternal(treeData, combineFact = combineFact)

  if(type == 'standard'){
   viviMat <-  viviBartStd(treeData = treeData,
                                 data = viviDf,
                                 metric = metric,
                                 reorder = reorder,
                                 combineFact = combineFact
                                 )
  }else if(type == 'vsup'){
    viviMat <- viviBartVSUP(treeData = treeData,
                                 data = viviDf,
                                 metricError = metricError,
                                 metric = metric,
                                 combineFact = combineFact,
                                 reorder = reorder)
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

viviBartInternal <- function(treeData, combineFact = FALSE){

  # Vimps -------------------------------------------------------------------

  # get vimps
  vimps <- bartMan::vimpBart(treeData, type = 'prop')
  vimpsVal <- bartMan::vimpBart(treeData, type = 'val')

  if(combineFact){
    vimps <- combineFactors(treeData = treeData, dataCombine = vimps)
    vimpsVal <- combineFactors(treeData = treeData, dataCombine =  vimpsVal)
  }

  vImp <- colMeans(vimps)
  vimpsVal <- colSums(vimpsVal)

  # get uncertainty measures
  vimpSD <- apply(vimps, 2, sd)
  upperVimp  <- vImp + 1.96 * vimpSD/sqrt(treeData$nMCMC)
  lowerVimp  <- vImp - 1.96 * vimpSD/sqrt(treeData$nMCMC)
  SEvimp <- sapply(as.data.frame(vimps), function(x) sd(x)/sqrt(length(x)))
  CVvimp <- vimpSD / vImp
  CVvimp[is.nan(CVvimp)] <- 0

  # get quantiles of proportions
  vimp25 <- apply(vimps, 2, function(x) quantile(x, c(.25)))
  vimp50 <- apply(vimps, 2, function(x) quantile(x, c(.50)))
  vimp75 <- apply(vimps, 2, function(x) quantile(x, c(.75)))

  # put together in dataframe
  vimpData <- cbind(vimpsVal, vImp, vimpSD, CVvimp, SEvimp,
                    lowerVimp, upperVimp, vimp25, vimp50, vimp75)

  # vimpData <- vimpData %>%
  #   as.data.frame() %>%
  #   tibble::rownames_to_column(var = "variable") %>%
  #   rename(count = vimpsVal, propMean = vImp, SD = vimpSD,
  #          CV = CVvimp,  SE = SEvimp,
  #          lowerCI = lowerVimp, upperCI = upperVimp,
  #          lowerQ = vimp25, median = vimp50, upperQ = vimp75)


  vDf <- vimpData %>%
    as.data.frame()

  vDf$variable <- row.names(vimpData)
  vimpData <- vDf |>
    rename(count = vimpsVal, propMean = vImp, SD = vimpSD,
           CV = CVvimp,  SE = SEvimp,
           lowerCI = lowerVimp, upperCI = upperVimp,
           lowerQ = vimp25, median = vimp50, upperQ = vimp75) |>
    select(variable, count, propMean, SD, CV, SE, lowerCI, upperCI, lowerQ, median, upperQ)



# VInts -------------------------------------------------------------------


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
  #vintNames <- reshape2::melt(dfVint)
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
  # dfProps <- reshape2::melt(propMatVintMean)
  # dfProps$var <- names(propMatVintMean)
  # colnames(dfProps) <- c( "props", "var")
  dfProps <- utils::stack(propMatVintMean)
  colnames(dfProps) = c('props', 'var')
  dfProps$var <- as.character(dfProps$var)



  # add counts
  countMean <- colMeans(dfVint)

  # turn into df
  # dfCountMean <- reshape2::melt(countMean)
  # dfCountMean$var <- names(propMatVintMean)
  # colnames(dfCountMean) <- c( "count", "var")
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
  # vintSD <- vintSD  |>
  #   reshape2::melt()  |>
  #   mutate(var = names(propMatVintMean))
  # colnames(vintSD) <- c( "SD", "var")
  vintSD <- utils::stack(vintSD)
  colnames(vintSD) = c('SD', 'var')
  vintSD$var <- as.character(vintSD$var)

  vintSE <- vintSD$SD/sqrt(treeData$nMCMC)
  names(vintSE) <- vintSD$var
  # vintSE <- vintSE  |>
  #   reshape2::melt()  |>
  #   mutate(var = names(propMatVintMean))
  # colnames(vintSE) <- c( "SE", "var")
  vintSE <- utils::stack(vintSE)
  colnames(vintSE) = c('SE', 'var')
  vintSE$var <- as.character(vintSE$var)


  # get quantiles of proportions
  vint25 <- apply(propMatVint, 2, function(x) quantile(x, c(.25)))
  vint50 <- apply(propMatVint, 2, function(x) quantile(x, c(.50)))
  vint75 <- apply(propMatVint, 2, function(x) quantile(x, c(.75)))

  # vint25 <- vint25  |> reshape2::melt()  |> mutate(var = names(propMatVintMean))
  # vint50 <- vint50  |> reshape2::melt()  |> mutate(var = names(propMatVintMean))
  # vint75 <- vint75  |> reshape2::melt()  |> mutate(var = names(propMatVintMean))

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
  #dfFinal$var <- unlist(dfFinal$var)

  # add coeff of variance
  dfFinal$CV <- dfFinal$SD / dfFinal$propMean
  dfFinal$CV[is.nan(dfFinal$CV)] <- 0

  # reorder
  dfFinal <-  dfFinal |> select(var, count, propMean, SD, CV, SE, Q25, Q50, Q75)

  # add adjustment
  vimpsAdj <- bartMan::vimpBart(treeData, type = 'propMean')

  splitN <- t(simplify2array(strsplit(as.character(dfFinal[["var"]]), ":")))
  values <- t(apply(splitN, 1, function(x){c(vimpsAdj[x[1]], vimpsAdj[x[2]])}))
  suppressMessages(
    res <- cbind(dfFinal, values)
  )
  names(res) <- c('var', 'count', 'propMean', 'SD', 'CV',
                  'SE', 'lowerQ', 'median', 'upperQ',
                  'propVimp1', 'propVimp2')

  trans <- apply(res[,c('propMean', 'propVimp1', 'propVimp2')], 1, function(x){
    (x[1]-x[2]*x[3])/sqrt(x[2]*x[3])
  })

  dfFinal$adjusted <- trans
  dfFinal$adjusted[dfFinal$adjusted <=  0] <- 0
  dfFinal$adjusted[is.na( dfFinal$adjusted)] <- 0

  names(dfFinal) <- c('var', 'count', 'propMean', 'SD', "CV",
                      'SE',  'lowerQ', 'median',
                      'upperQ', 'adjusted')

  myList <- list(Vimp = vimpData, Vint = dfFinal)

  return(myList)

}





# -------------------------------------------------------------------------
# Standard Matrix ---------------------------------------------------------
# -------------------------------------------------------------------------

viviBartStd <- function(treeData,
                        data,
                        reorder = FALSE,
                        metric = "propMean",
                        combineFact = FALSE){

  propFinal <- data$Vint
  vars2  <- t(simplify2array(strsplit(as.character(propFinal[["var"]]), ":")))
  if(combineFact){
    vimps <- bartMan::vimpBart(treeData, type = 'prop')
    vimps <- combineFactors(treeData = treeData, dataCombine = vimps)
    ovars <- names(vimps)
  }else{
    ovars <- treeData$varName
  }
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



viviBartVSUP <- function(treeData,
                         data,
                         reorder = FALSE,
                         metricError = 'CV',
                         metric = "propMean",
                         combineFact = FALSE){

  # get matrix of uncertainty values
  propFinal <- data$Vint
  vars2  <- t(simplify2array(strsplit(as.character(propFinal[["var"]]), ":")))

  if(combineFact){
    vimps <- bartMan::vimpBart(treeData, type = 'prop')
    vimps <- combineFactors(treeData = treeData, dataCombine = vimps)
    ovars <- names(vimps)
  }else{
    ovars <- treeData$varName
  }

  mat <- matrix(0, length(ovars), length(ovars)) # create matrix
  rownames(mat) <- colnames(mat) <- ovars # set names
  mat[vars2] <- propFinal[[metricError]] # set values
  mat <- (mat[lower.tri(mat, diag = T)[nrow(mat):1], ] +
            t(mat)[lower.tri(mat, diag = T)[nrow(mat):1], ])



  vimps <- data$Vimp # add vimps to matirx
  diag(mat) <- vimps[,metricError]
  uncertaintyMatrix <- mat

  # get actual values matrix
  actualMatrix <- viviBartStd(treeData = treeData, data = data,
                              reorder = reorder, metric = metric,
                              combineFact = combineFact)

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


viviBartQuantile <-  function(treeData, data, reorder = FALSE){

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








