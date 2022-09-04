#' viviBart
#'
#' @description Returns a list containing a dataframe of variable importance summaries
#' and a dataframe of variable interaction summaries.
#'
#' @param treeData A data frame created by extractTreeData function.
#' @param combineFact If a variable is a factor in a data frame, when building the BART model it is replaced with dummies.
#' Note that q dummies are created if q>2 and one dummy is created if q=2, where q is the number of levels of the factor.
#' If combineFact = TRUE, then both the importance and interactions are calculated for the entire factor by aggregating the dummy variables’
#' inclusion proportions.
#' @param out Choose to either output just the variable importance ('vimp'), t
#' he variable interaction ('vint'), or both ('vivi') (default).
#'
#' @return A list of dataframes of VIVI summaries.
#'
#'
#' @export
#'

viviBart <- function(treeData, combineFact = FALSE, out = 'vivi'){

  vivis <- viviBartInternal(treeData = treeData, combineFact = combineFact)


    if(out == 'vimp'){
      out <-  vivis$Vimp
    }else if(out == 'vint'){
      out <- vivis$Vint
    }else if(out == 'vivi'){
      out <- vivis
    }


  return(out)



}
