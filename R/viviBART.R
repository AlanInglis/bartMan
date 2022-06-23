#' viviBart
#'
#' @description Returns a list containing a dataframe of variable importance summaries
#' and a dataframe of variable interaction summaries.
#'
#' @param treeData A data frame created by extractTreeData function.
#'
#' @return A list of dataframes of VIVI summaries.
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom dplyr mutate
#' @importFrom dplyr distinct
#' @importFrom dplyr select
#' @importFrom dplyr rename
#' @importFrom tibble rownames_to_column
#' @importFrom stringr str_extract_all
#' @importFrom stringr str_glue
#' @importFrom purrr map_chr
#'
#' @export
#'

viviBart <- function(treeData){

  # Vimps -------------------------------------------------------------------

  # get vimps
  vimps <- bartMan::vimpBart(treeData, type = 'prop')
  vimpsVal <- bartMan::vimpBart(treeData, type = 'val')
  #propVimp <- proportions(vimps, 1)
  vImp <- colMeans(vimps)
  vimpsVal <- colSums(vimpsVal)

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


  # add in quantile, sd, and se columns
  propMM$lowerQ <- vint25
  propMM$median <- vint50
  propMM$upperQ <- vint75
  propMM$SEofCI <- SEvint
  propMM$SD <- vintSD

  # make symmetrical
  propMM <- propMM %>%
    mutate(var = map(
      stringr::str_split(var, pattern = ":"),
      ~ sort(.x) %>% trimws(.) %>% paste0(., collapse = ':')
    ))

  propMM$var <- as.character(propMM$var)

  propFinal <- propMM %>%
    group_by(var) %>%
    mutate(propMean = sum(value))  %>%
    mutate(count = sum(count)) %>%
    mutate(SD = mean(SD)) %>%
    mutate(lowerQ = mean(lowerQ)) %>%
    mutate(median = mean(median)) %>%
    mutate(upperQ = mean(upperQ)) %>%
    mutate(SEofCI = mean(SEofCI)) %>%
    select(var, count, propMean, SD, lowerQ, median, upperQ, SEofCI) %>%
    ungroup %>%
    distinct() %>%
    arrange(-count)

  # add adjustment
  vimps <- bartMan::vimpBart(treeData, type = 'propMean')

  splitN <- t(simplify2array(strsplit(as.character(propFinal[["var"]]), ":")))
  values <- t(apply(splitN, 1, function(x){c(vimps[x[1]], vimps[x[2]])}))
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
