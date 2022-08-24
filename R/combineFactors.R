#' combineFactors
#'
#' @description If a variable is a factor in a data frame, when building the BART model it is replaced with dummies.
#' Note that q dummies are created if q>2 and one dummy is created if q=2, where q is the number of levels of the factor.
#' These functions combine the factor levels so that the inclusion proportions
#' are aggregated so the importance can be assessed for the the entire factor.
#'
#' @param treeData A data frame created by extractTreeData function.
#' @param data2 A second data frame that has the factor split into individual levels.
#' Used to compare the split data frame with the original data frame.
#'
#'
#' @export

combineFactors <- function(treeData, data2){

  df  <- treeData$data
  dff <- names(which(!(sapply(df[colnames(df)], is.numeric))))

  whichCols <- which(!(colnames(df) %in% colnames(data2)))
  factorColNam <- names(which(!(sapply(df[colnames(df)], is.numeric))))
  factorCols <- which((colnames(df) %in% dff))
  dfnew <- list()

  # for BART
  if(any(class(treeData) == 'bart')){

    for (i in 1:length(factorCols)) {
      facLevels <- unique(df[ ,factorCols[i]])
      dfnew[[i]] <- paste0(factorColNam[i],  as.numeric(facLevels))
    }


  }else if(any(class(treeData) == 'dbarts')){

    for (i in 1:length(factorCols)) {
      facLevels <- unique(df[,factorCols[i]])
      dfnew[[i]] <- paste0(factorColNam[i], ".", facLevels)
    }

    only2Factors <- which(do.call(rbind, lapply(dfnew, length)) == 2)

    for(i in only2Factors){
      dfnew[[i]] <- NULL#dfnew[[i]][-1]
    }

  }else if(any(class(treeData) == 'bartMach')){

    for (i in 1:length(factorCols)) {
      facLevels <- unique(df[,factorCols[i]])
      dfnew[[i]] <- paste0(factorColNam[i], "_", facLevels)
    }

  }


  factList <- lapply(dfnew, function(x){
    apply(data2[ ,colnames(data2) %in% x], 1, sum)
  })

  names(factList) <- factorColNam

  df2 <- cbind(data2, as.data.frame(factList))
  df2 <- df2[, (colnames(df2) %in% colnames(df))]

  return(df2)

}

#' @export
combineFactorsDiag <- function(data, df2, model){

df  <- data
dff <- names(which(!(sapply(df[colnames(df)], is.numeric))))

whichCols <- which(!(colnames(df) %in% colnames(df2)))
factorColNam <- names(which(!(sapply(df[colnames(df)], is.numeric))))
factorCols <- which((colnames(df) %in% dff))
dfnew <- list()

# for BART
if(class(model) == 'wbart' || class(model) == 'pbart'){

  for (i in 1:length(factorCols)) {
    facLevels <- unique(df[ ,factorCols[i]])
    dfnew[[i]] <- paste0(factorColNam[i],  as.numeric(facLevels))
  }


}else if(class(model) == 'bart'){

  for (i in 1:length(factorCols)) {
    facLevels <- unique(df[,factorCols[i]])
    dfnew[[i]] <- paste0(factorColNam[i], ".", facLevels)
  }

  only2Factors <- which(do.call(rbind, lapply(dfnew, length)) == 2)

  for(i in only2Factors){
    dfnew[[i]] <- NULL#dfnew[[i]][-1]
  }

}else if(class(model) == 'bartMachine'){

  for (i in 1:length(factorCols)) {
    facLevels <- unique(df[,factorCols[i]])
    dfnew[[i]] <- paste0(factorColNam[i], "_", facLevels)
  }

}


factList <- lapply(dfnew, function(x){
  apply(df2[ ,colnames(df2) %in% x], 1, sum)
})

names(factList) <- factorColNam

df2 <- cbind(df2, as.data.frame(factList))
df2 <- df2[, (colnames(df2) %in% colnames(df))]

return(df2)

}

#' @export

# dfCompare <- vimpsTest
combineFactorsInt <- function(propData, dataOG){

  propData <- propData %>%
    as_tibble(rownames = "id") %>%
    pivot_longer(-id) %>%
    mutate(name = map_chr(str_extract_all(name, paste(c(colnames(dataOG), ":"), collapse = "|")),
                          paste0, collapse = "")) %>%
    group_by(id, name) %>%
    summarise(value = sum(value)) %>%
    pivot_wider()


  propData <- propData[,-1]
  propData <- as.matrix(propData)


  return(propData)

}

