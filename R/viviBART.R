#' viviBart
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
#' @export


viviBart <- function(model,
                     response,
                     data,
                     noReplications = 2,
                     type = NULL,
                     gridSize = 50,
                     nmax = 500,
                     reorder = TRUE,
                     class = 1,
                     normalized = FALSE){

  if(type == 'vsup'){

    viviMat <- viviBartV(model = model,
                         respons = response,
                         data = data,
                         noReplications = noReplications,
                         gridSize = gridSize,
                         nmax = nmax,
                         reorder = reorder,
                         class = class,
                         normalized = normalized)

      class(viviMat) <- c("list", "vsup")
      return(viviMat)

  }else if(type == 'quant'){

    viviMat <- viviBartQ(model = model,
                         respons = response,
                         data = data,
                         noReplications = noReplications,
                         gridSize = gridSize,
                         nmax = nmax,
                         reorder = reorder,
                         class = class,
                         normalized = normalized)

    class(viviMat) <- c("list", "quant")
    return(viviMat)
  }
}



# VSUP function for viviBart ----------------------------------------------


  viviBartV <- function(model,
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

    # check if classification or not
    responseIdx <- which(colnames(data) == response)
    if(class(model) == "wbart" || class(model) == "pbart"){
      pFun <- function(fit, data, prob=TRUE) as.numeric(condvis2::CVpredict(fit, data[,-responseIdx]))
    } else if(class(model) == 'bart'){
      if(is.null(model$sigma)){
        pFun <- function(fit, data, prob=TRUE) apply(predict(fit, data[,-responseIdx]), 2, mean)
      } else {
        pFun <- NULL
      }
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
    vimps <- dplyr::bind_rows(vimp)
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

    vInt <- lapply(mat, getIntValues) %>% dplyr::bind_rows()
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

    class(actualVals) <- c('vivid', 'matrix', 'array', 'vsup')
    class(uncMat) <- c('vivid', 'matrix', 'array', 'vsup')

    myList <- list(actualMatrix = actualVals,
                   uncertaintyMatrix = uncMat)

    return(myList)

  }


# Quantile function for viviBart ------------------------------------------

  viviBartQ <- function(model,
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

    # check if classification or not
    responseIdx <- which(colnames(data) == response)
    if(class(model) == "wbart" || class(model) == "pbart"){
      pFun <- function(fit, data, prob=TRUE) as.numeric(condvis2::CVpredict(fit, data[,-responseIdx]))
    } else if(class(model) == 'bart'){
      if(is.null(model$sigma)){
        pFun <- function(fit, data, prob=TRUE) apply(predict(fit, data[,-responseIdx]), 2, mean)
      } else {
        pFun <- NULL
      }
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
    vimps <- dplyr::bind_rows(vimp)
    vimpsQuant <- apply(vimps, 2, function(x) quantile(x, c(0.05, 0.5, 0.95)))

    #avgVimp <- apply(vimps, 2, mean) # average importance

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

    vInt <- lapply(mat, getIntValues) %>% dplyr::bind_rows()
    vintsQuant <- apply(vInt, 2, function(x) quantile(x, c(0.05, 0.5, 0.95)))

    #avgVint <- apply(vInt, 2, mean) # average interaction


    # spereate vImp and vInt into each quantile
    vimpQuant.05 <- vimpsQuant[1,]
    vimpQuant.50 <- vimpsQuant[2,]
    vimpQuant.95 <- vimpsQuant[3,]

    vintQuant.05 <- vintsQuant[1,]
    vintQuant.50 <- vintsQuant[2,]
    vintQuant.95 <- vintsQuant[3,]


    # Turn back into matrix ---------------------------------------------------


    matrixTrans <- function(vintData, vimpData){
      mat <- read.table(text = names(vintData))
      names(mat) <- c("Var1", "Var2")
      mat <- xtabs(vintData ~ ., cbind(rbind(mat, setNames(rev(mat), names(mat))), vintData = rep(vintData, 2)))
      diag(mat) <- vimpData

      return(mat)
    }

    matQuant.05 <- matrixTrans(vintQuant.05, vimpQuant.05)
    matQuant.50 <- matrixTrans(vintQuant.50, vimpQuant.50)
    matQuant.95 <- matrixTrans(vintQuant.95, vimpQuant.95)


    # reorder median matrix
    if(reorder){
      matQuant.50 <- vivid::vividReorder(matQuant.50)
      # match other matrices with new order
      ValsColOrder <- colnames(matQuant.50)
      matQuant.05 <- matQuant.05[ValsColOrder, ValsColOrder]
      matQuant.95 <- matQuant.95[ValsColOrder, ValsColOrder]
    }


    class(matQuant.05) <- c('vivid', 'matrix', 'array', 'quant')
    class(matQuant.50) <- c('vivid', 'matrix', 'array', 'quant')
    class(matQuant.95) <- c('vivid', 'matrix', 'array', 'quant')

    myList <- list(quant.05 = matQuant.05,
                   quant.50 = matQuant.50,
                   quant.95 = matQuant.95)

    return(myList)

  }





