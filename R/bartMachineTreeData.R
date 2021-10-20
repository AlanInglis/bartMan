#' bartMachineTreeData
#'
#' @description Creates a data frame of all tree attributes for a model
#' created by the bartMachine package.
#'
#' @param model Model created from the bartMachine package.
#'
#' @return A data frame of every tree and its attributes.
#'
#'
#' @importFrom bartMachine extract_raw_node_data
#' @importFrom rrapply rrapply
#' @importFrom purrr transpose
#' @importFrom purrr map
#' @importFrom purrr flatten
#' @importFrom stringr str_c
#' @importFrom dplyr arrange
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#'
#'
#' @export


bartMachineTreeData <- function(model){

  # extract the raw node data
  nodeData <- bartMachine::extract_raw_node_data(bM)
  listNodesBM <- list(nodeData)

  # bind node data together
  out <- rrapply::rrapply(listNodesBM, how = 'bind')

  # get position of grep
  i1 <- grep('string_location', names(out))

  # extract relevant info:
  resNodeID <- out[i1] %>%
    transpose %>%
    purrr::map( ~ keep(.x, lengths(.x) > 0)) %>%
    flatten %>%
    setNames(str_c('NodeID.', seq_along(.)))

  resLeaf <- out[i1-5]%>%
    transpose %>%
    purrr::map( ~ keep(.x, lengths(.x) > 0)) %>%
    flatten %>%
    setNames(str_c('isLeaf.', seq_along(.)))

  resSplitValue <- out[i1+2]%>%
    transpose %>%
    purrr::map( ~ keep(.x, lengths(.x) > 0)) %>%
    flatten %>%
    setNames(str_c('splitValue.', seq_along(.)))

  resStump <- out[i1-1]%>%
    transpose %>%
    purrr::map( ~ keep(.x, lengths(.x) > 0)) %>%
    flatten %>%
    setNames(str_c('isStump.', seq_along(.)))

  resYPred <- out[i1+3] %>%
    transpose %>%
    purrr::map( ~ keep(.x, lengths(.x) > 0)) %>%
    flatten %>%
    setNames(str_c('isStump.', seq_along(.)))

  myLists <- list(resNodeID, resSplitValue, resLeaf, resYPred, resStump)

  # Turn into dataframe
  df <- data.frame(matrix(unlist(myLists), nrow = length(myLists[[1]]), byrow=FALSE))
  names(df) <- c("nodeID", "splitValue", "isLeaf", "leafValue", "isStump")
  df <- df[rowSums(is.na(df)) != ncol(df), ] # remove NA rows

  # Add tree number
  df$treeNum <-cumsum(df$nodeID=="P")

  # rearrange the rows (ie nodes) so the left nodes come before
  # the right nodes
  df <- df %>%
    group_by(treeNum) %>%
    mutate(level = nchar(nodeID)) %>%
    group_by(treeNum) %>%
    arrange(level, ifelse(nodeID=="P", 1, match(LETTERS, nodeID)+1),
            .by_group = TRUE)


  # define node number sequentially
  df <- df %>%
    group_by(treeNum) %>%
    mutate(node = row_number())

  # get the parent node
  df <- df %>%
    group_by(treeNum) %>%
    mutate(parentNode = substr(nodeID, 0, nchar(nodeID)-1))

  # define parent node of P as NA
  df$parentNode[df$nodeID == "P"] <- NA

  # define parent nodes nodeID with still empty parent node name as "P"
  df$parentNode[df$parentNode==""] <- "P"

  # Match parent node names to node numbers
  df <- df %>%
    group_by(treeNum) %>%
    mutate(parentNodeNo = match(parentNode, nodeID))

  # Add iteration column
  df$iteration <- 1

  # add value column
  df <-  df %>%
    mutate(value = coalesce(splitValue, leafValue))

  # add label column
  df <- transform(df, label = ifelse(is.na(splitValue), value, paste(nodeID, value, sep = " ≤ ")))


  return(df)

}
