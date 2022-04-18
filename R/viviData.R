#' vimpBart
#'
#' @description A matrix with nMCMC rows with each variable as a column.
#' Each row represents an MCMC iteration. For each variable, the total count
#' of the number of times that variable is used in a tree is given.
#'
#' @param treeData A data frame created by treeData function.
#' @param type What value to return. Either the raw count 'val', the proportion 'prop',
#' the column means of the proportions 'propMean', or the median of the proportions 'propMedian'.
#'
#' @return A matrix
#'
#' @importFrom dplyr ungroup
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr select
#'
#' @export
#'


vimpBart <- function(treeData, type = 'prop'){

  if (!(type %in% c("val", "prop", "propMean", "propMedian"))) {
    stop("type must be \"val\", \"prop\", \"propMedian\", or \"propMean\"")
  }

  df <- treeData$structure

  # get count of vars
  vCount <- df %>%
    ungroup() %>%
    select(iteration, var) %>%
    group_by(iteration) %>%
    table() %>%
    as.data.frame.matrix()

  # turn into matrix (with all values)
  nam <- treeData$varName
  mat <- matrix(0, nrow = treeData$nMCMC, ncol = length(nam))
  colnames(mat) <- nam

  namesDat <- colnames(vCount)
  matchCol <- which(nam %in% namesDat)

  for (i in matchCol) {
    mat[, nam[i]] <- vCount[nam[i]][[nam[i]]]
  }

  if(type == 'prop'){
    mat <- proportions(mat, 1)
  }else if(type == 'propMean'){
    mat <- proportions(mat, 1)
    mat <- colMeans(mat)
  }else if(type == 'propMedian'){
    mat <- proportions(mat, 1)
    mat <- apply(mat, 2, median)
  }

  return(mat)
}


#' vintBart
#'
#' @description A data frame displaying the count of each pair of variables used
#' for splitting within a tree.
#'
#'
#' @param treeData A data frame created by treeData function.
#'
#' @return A dataframe containing the count of each pair of variables and the mean of the proportions.
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr select
#' @importFrom dplyr count
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr rename
#' @importFrom dplyr arrange
#' @importFrom stringr str_detect
#' @importFrom stringr str_extract_all
#' @importFrom stringr str_glue
#'
#' @export
#'
#'


vintBart <- function(treeData){

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

  propMatVint <- proportions(dfVint, 1)
  propMatVintMean <- colMeans(propMatVint)
  propMatVintMedian <- apply(propMatVint, 2, median)

  # get quantiles of proportions
  vint25 <- apply(propMatVint, 2, function(x) quantile(x, c(.25)))
  vint50 <- apply(propMatVint, 2, function(x) quantile(x, c(.50)))
  vint75 <- apply(propMatVint, 2, function(x) quantile(x, c(.75)))

  propMM <- reshape::melt(propMatVintMean)
  propMM <- tibble::rownames_to_column(propMM, "var")

  countVint <- colSums(dfVint)
  countM <- reshape::melt(countVint)
  countM <- tibble::rownames_to_column(countM, "var")

  propMM$count <- countM$value

  propMM <- propMM %>%
    mutate(var = map(
      stringr::str_split(var, pattern = ":"),
      ~ sort(.x) %>% trimws(.) %>% paste0(., collapse = ':')
    ))

  propMM$var <- as.character(propMM$var)

  # add in quantile columns

propMM$lowerQ <- vint25
propMM$median <- vint50
propMM$upperQ <- vint75

  propFinal <- propMM %>%
    group_by(var) %>%
    mutate(propMean = mean(value))  %>%
    mutate(countFinal = sum(count)) %>%
    select(var, countFinal, propMean, lowerQ, median, upperQ) %>%
    distinct() %>%
    arrange(-countFinal)

  # add adjustment
  vimps <- bartMan::vimpBart(treeData, type = 'propMean')

  splitN <- t(simplify2array(strsplit(as.character(propFinal[["var"]]), ":")))
  values <- t(apply(splitN, 1, function(x){c(vimps[x[1]], vimps[x[2]])}))
  suppressMessages(
    res <- cbind(propFinal, values)
  )
  names(res) <- c('var', 'countFinal', 'meanProp','lowerQ', 'median',
                  'upperQ', 'propVimp1', 'propVimp2')

  trans <- apply(res[,c('meanProp', 'propVimp1', 'propVimp2')], 1, function(x){
    (x[1]-x[2]*x[3])/sqrt(x[2]*x[3])
  })

  propFinal$adjusted <- trans

  return(propFinal)
}


# vintPlot ----------------------------------------------------------------


vintPlot <- function(treeData, data, type = 'val', transform  = FALSE, vimps = NULL){

  if (!(type %in% c("val", "propMean"))) {
    stop("type must be \"val\"  or \"propMean\"")
  }


  res <- data
  ovars <- treeData$varName
  if(is.null(vimps)){
    vimps <- bartMan::vimpBart(treeData, type = 'propMean')
  }

  if(type == 'val'){

    res[["var"]] <- reorder(res[["var"]], res[["count"]])
    # create matrix of values
    vars2 <- t(simplify2array(strsplit(as.character(res[["var"]]), ":"))) # split/get feature names
    mat <- matrix(0, length(ovars), length(ovars)) # create matrix
    rownames(mat) <- colnames(mat) <- ovars # set names
    mat[vars2] <- res[["count"]] # set values

  }else if(type == 'propMean'){

    if(transform){
      splitN <- t(simplify2array(strsplit(as.character(res[["var"]]), ":")))
      values <- t(apply(splitN, 1, function(x){c(vimps[x[1]], vimps[x[2]])}))
      suppressMessages(
        res <- cbind(res, values)
      )
      names(res) <- c('var', 'meanProp', 'prop1', 'prop2')

      trans <- apply(res[,c('meanProp', 'prop1', 'prop2')], 1, function(x){
        (x[1]-x[2]*x[3])/sqrt(x[2]*x[3])
      })

      res$transformed <- trans

      res[["var"]] <- reorder(res[["var"]], res[["transformed"]])
      # create matrix of values
      vars2 <- t(simplify2array(strsplit(as.character(res[["var"]]), ":"))) # split/get feature names
      mat <- matrix(0, length(ovars), length(ovars)) # create matrix
      rownames(mat) <- colnames(mat) <- ovars # set names
      mat[vars2] <- res[["transformed"]] # set values
    }else{
      res[["var"]] <- reorder(res[["var"]], res[["propMean"]])
      # create matrix of values
      vars2 <- t(simplify2array(strsplit(as.character(res[["var"]]), ":"))) # split/get feature names
      mat <- matrix(0, length(ovars), length(ovars)) # create matrix
      rownames(mat) <- colnames(mat) <- ovars # set names
      mat[vars2] <- res[["propMean"]] # set values
    }
  }


  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  diag(mat) <- vimps
  return(mat)
}
