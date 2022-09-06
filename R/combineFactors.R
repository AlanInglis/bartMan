#' combineFactors
#'
#' @description If a variable is a factor in a data frame, when building the BART model it is replaced with dummies.
#' Note that q dummies are created if q>2 and one dummy is created if q=2, where q is the number of levels of the factor.
#' This function combines the factor levels so that the inclusion proportions
#' are aggregated so the importance can be assessed for the the entire factor.
#'
#' @param dataCombine A data frame with dummy factors to be combined into single factors.
#' @param treeData A data frame created by extractTreeData function.
#'
#' @export

combineFactors <- function(dataCombine, treeData){

  #df  <- treeData$data
  df <- treeData
  factorColNam <- names(which(!(sapply(df[colnames(df)], is.numeric))))
  factorCols <- which((colnames(df) %in% factorColNam))
  dfnew <- list()

  # for BART
  if(any(class(treeData) == 'bart')){

    for (i in 1:length(factorCols)) {
      facLevels <- unique(df[ ,factorCols[i]])
      dfnew[[i]] <- paste0(factorColNam[i],  as.numeric(facLevels))
    }


    factList <- lapply(dfnew, function(x){
      apply(dataCombine[ ,colnames(dataCombine) %in% x], 1, sum)
    })

    names(factList) <- factorColNam

    df2 <- cbind(dataCombine, as.data.frame(factList))
    df2 <- df2[, (colnames(df2) %in% colnames(df))]

  }else if(any(class(treeData) == 'dbarts')){


    dfnew <- list()
    for (i in 1:length(factorCols)) {
      facLevels <- unique(df[,factorCols[i]])
      dfnew[[i]] <- paste0(factorColNam[i], ".", facLevels)
    }

    # find if any variables have only 2 factors
    only2Factors <- which(do.call(rbind, lapply(dfnew, length)) == 2)
    singleFact <- dfnew[only2Factors]

    # remove from list
    for(i in  seq_along(only2Factors)){
      dfnew[[i]] <- NULL#dfnew[[i]][-1]
    }

    factList <- lapply(dfnew, function(x){
      apply(dataCombine[ ,colnames(dataCombine) %in% x], 1, sum)
    })

    names(factList) <- factorColNam[-only2Factors]

    # rename single factor columns
    sft <- list()
    for(i in seq_along(singleFact)){
      sft[[i]] <- dataCombine[,singleFact[[i]][length(singleFact[[i]])]]
    }

    singleFactIdx <- which(!(factorColNam %in% names(factList)))
    singleFactNames <- factorColNam[singleFactIdx]
    names(sft) <- singleFactNames

    df2 <- cbind(dataCombine, as.data.frame(factList), as.data.frame(sft))
    df2 <- df2[, (colnames(df2) %in% colnames(df))]
    df2 <- df2[ , order(names(df2))]

  }else if(any(class(treeData) == 'bartMach')){

    for (i in 1:length(factorCols)) {
      facLevels <- unique(df[,factorCols[i]])
      dfnew[[i]] <- paste0(factorColNam[i], "_", facLevels)
    }


    factList <- lapply(dfnew, function(x){
      apply(dataCombine[ ,colnames(dataCombine) %in% x], 1, sum)
    })

    names(factList) <- factorColNam

    df2 <- cbind(dataCombine, as.data.frame(factList))
    df2 <- df2[, (colnames(df2) %in% colnames(df))]


  }

  return(df2)

}

#' combineFactorsDiag
#'
#' @description If a variable is a factor in a data frame, when building the BART model it is replaced with dummies.
#' Note that q dummies are created if q>2 and one dummy is created if q=2, where q is the number of levels of the factor.
#' These functions combine the factor levels so that the inclusion proportions
#' are aggregated so the importance can be assessed for the the entire factor.
#'
#' @param data A data frame
#' @param dataCombine A data frame with dummy factors to be combined into single factors.
#' @param model a BART model
#'
#'
#' @export
combineFactorsDiag <- function(data, dataCombine, model){

  df  <- data
  factorColNam <- names(which(!(sapply(df[colnames(df)], is.numeric))))
  factorCols <- which((colnames(df) %in% factorColNam))
  dfnew <- list()

  # for BART
  if(class(model) == 'wbart' || class(model) == 'pbart'){

    for (i in 1:length(factorCols)) {
      facLevels <- unique(df[ ,factorCols[i]])
      dfnew[[i]] <- paste0(factorColNam[i],  as.numeric(facLevels))
    }

    factList <- lapply(dfnew, function(x){
      apply(dataCombine[ ,colnames(dataCombine) %in% x], 1, sum)
    })

    names(factList) <- factorColNam

    df2 <- cbind(dataCombine, as.data.frame(factList))
    df2 <- df2[, (colnames(df2) %in% colnames(df))]


  }else if(class(model) == 'bart'){


    for (i in 1:length(factorCols)) {
      facLevels <- unique(df[,factorCols[i]])
      dfnew[[i]] <- paste0(factorColNam[i], ".", facLevels)
    }

    # find if any variables have only 2 factors
    only2Factors <- which(do.call(rbind, lapply(dfnew, length)) == 2)
    singleFact <- dfnew[only2Factors]

    # remove from list
    for(i in  seq_along(only2Factors)){
      dfnew[[i]] <- NULL#dfnew[[i]][-1]
    }

    factList <- lapply(dfnew, function(x){
      apply(dataCombine[ ,colnames(dataCombine) %in% x], 1, sum)
    })

    names(factList) <- factorColNam[-only2Factors]

    # rename single factor columns
    sft <- list()
    for(i in seq_along(singleFact)){
      sft[[i]] <- dataCombine[,singleFact[[i]][length(singleFact[[i]])]]
    }

    singleFactIdx <- which(!(factorColNam %in% names(factList)))
    singleFactNames <- factorColNam[singleFactIdx]
    names(sft) <- singleFactNames


    df2 <- cbind(dataCombine, as.data.frame(factList), as.data.frame(sft))
    df2 <- df2[, (colnames(df2) %in% colnames(df))]
    df2 <- df2[ , order(names(df2))]

  }else if(class(model) == 'bartMachine'){

    for (i in 1:length(factorCols)) {
      facLevels <- unique(df[,factorCols[i]])
      dfnew[[i]] <- paste0(factorColNam[i], "_", facLevels)
    }

    factList <- lapply(dfnew, function(x){
      apply(dataCombine[ ,colnames(dataCombine) %in% x], 1, sum)
    })

    names(factList) <- factorColNam

    df2 <- cbind(dataCombine, as.data.frame(factList))
    df2 <- df2[, (colnames(df2) %in% colnames(df))]

  }


  return(df2)

}



#' combineFactorsInt
#'
#' @description If a variable is a factor in a data frame, when building the BART model it is replaced with dummies.
#' Note that q dummies are created if q>2 and one dummy is created if q=2, where q is the number of levels of the factor.
#' These functions combine the factor levels so that the inclusion proportions
#' are aggregated so the importance can be assessed for the the entire factor.
#'
#' @param dataInt A data frame of interactions
#' @param data Original data frame used to build model.
#'
#' @importFrom dplyr %>%
#' @importFrom stats na.omit
#'
#' @export

combineFactorsInt <- function(dataInt, data){

  str_ext <- function(string, pattern) {
    regmatches(string, gregexpr(pattern, string))
  }


  dataInt <- dataInt %>%
    dplyr::as_tibble(rownames = "id") %>%
    tidyr::pivot_longer(-id) %>%
    dplyr::mutate(name = purrr::map_chr(str_ext(name, paste(c(colnames(data), ":"), collapse = "|")),
                                        paste0, collapse = "")) %>%
    dplyr::group_by(id, name) %>%
    dplyr::summarise(value = sum(value)) |>
    dplyr::na_if("") %>%
    na.omit %>%
    tidyr::pivot_wider()


  dataInt <- dataInt[,-1]
  dataInt <- as.matrix(dataInt)


  return(dataInt)

}

