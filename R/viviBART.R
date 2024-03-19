#' viviBart
#'
#' @description Returns a list containing a dataframe of variable importance summaries
#' and a dataframe of variable interaction summaries.
#'
#' @param treeData A data frame created by extractTreeData function.
#' @param out Choose to either output just the variable importance ('vimp'),
#' the variable interaction ('vint'), or both ('vivi') (default).
#'
#' @return A list of dataframes of VIVI summaries.
#'
#'
#' @export
#'

viviBart <- function(treeData, out = 'vivi'){

  vivis <- viviBartInternal(treeData = treeData)

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
