#' viviBart
#'
#' @description Returns a list containing a dataframe of variable importance summaries
#' and a dataframe of variable interaction summaries.
#'
#' @param trees A data frame created by `extractTreeData` function.
#' @param out Choose to either output just the variable importance ('vimp'),
#' the variable interaction ('vint'), or both ('vivi') (default).
#'
#' @return A list of dataframes of VIVI summaries.
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
#'  viviBart(trees = trees_data, out = 'vivi')
#'  }
#'
#'
#' @export
#'

viviBart <- function(trees, out = 'vivi'){

  vivis <- viviBartInternal(trees = trees)

  if (!(out %in% c("vimp", "vint", "vivi"))) {
    stop("out must be \"vimp\", \"vint\"  or \"vivi\"")
  }

    if(out == 'vimp'){
      out <-  vivis$Vimp
    }else if(out == 'vint'){
      out <- vivis$Vint
    }else if(out == 'vivi'){
      out <- vivis
    }


  return(out)



}
