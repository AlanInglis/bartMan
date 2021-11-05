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

  # Get variable names
  varNames <- colnames(model$X)

  # Get No of iterations after burn in
  iter <- model$num_iterations_after_burn_in

  # extract the raw node data for all iterations
  nodeData <- NULL
  if (iter > 1) {
    nodeData <- NULL
    for (i in 1:iter) {
      nodeData[[i]] <- bartMachine::extract_raw_node_data(model, g = i)
    }
  } else {
    nodeData <- bartMachine::extract_raw_node_data(model)
  }


  # bind node data together
  out <- rrapply::rrapply(nodeData, how = 'bind')

  # get position of grep
  i1 <- grep('string_location', names(out))

  # extract relevant info:
  resVar <- out[i1+1] %>%
    transpose %>%
    purrr::map( ~ keep(.x, lengths(.x) > 0)) %>%
    flatten %>%
    setNames(str_c('var', seq_along(.)))

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

  myLists <- list(resVar, resNodeID, resSplitValue, resLeaf, resYPred, resStump)

  # Turn into dataframe
  df <- data.frame(matrix(unlist(myLists), nrow = length(myLists[[1]]), byrow=FALSE))
  names(df) <- c("varID", "nodeID", "splitValue", "isLeaf", "leafValue", "isStump")
  df <- df[rowSums(is.na(df)) != ncol(df), ] # remove NA rows

  # match var number to varName
  names(varNames) <- c(0:(length(varNames)-1))
  df$var <- varNames[as.character(df$varID)]

  # Add tree number (ignoring iterration)
  df$treeNumID <-cumsum(df$nodeID=="P")

  # rearrange the rows (ie nodes) so the left nodes come before
  # the right nodes
  df <- df %>%
    group_by(treeNumID) %>%
    mutate(level = nchar(nodeID)) %>%
    group_by(treeNumID) %>%
    arrange(level, ifelse(nodeID=="P", 1, match(LETTERS, nodeID)+1),
            .by_group = TRUE)


  # define node number sequentially
  df <- df %>%
    group_by(treeNumID) %>%
    mutate(node = row_number())

  # get the parent node
  df <- df %>%
    group_by(treeNumID) %>%
    mutate(parentNode = substr(nodeID, 0, nchar(nodeID)-1))

  # define parent node of P as NA
  df$parentNode[df$nodeID == "P"] <- NA

  # define parent nodes nodeID with still empty parent node name as "P"
  df$parentNode[df$parentNode==""] <- "P"

  # Match parent node names to node numbers
  df <- df %>%
    group_by(treeNumID) %>%
    mutate(parentNodeNo = match(parentNode, nodeID))

  # Add iteration column
  noTrees <- model$num_trees
  dfIteration <- data.frame(nodeID = df$nodeID)
  dfIteration <- dfIteration %>%
    mutate(iter = replace(rep(NA_integer_, n()), nodeID == 'P',
                          as.integer(gl(sum(nodeID == 'P'), noTrees, sum(nodeID == 'P'))))) %>%
    fill(iter)

  df$iteration <- dfIteration$iter

  # round values
  df$splitValue <- round(as.numeric(df$splitValue),4)
  df$leafValue <- round(as.numeric(df$leafValue),  4)

  # add value column
  df <-  df %>%
    mutate(value = coalesce(splitValue, leafValue))

  # add label column
  df <- transform(df, label = ifelse(is.na(splitValue), value, paste(var, value, sep = " ≤ ")))

  # Change treeNum to reflect the iteration
  df <- df %>%
    group_by(iteration) %>%
    mutate(treeNum = cumsum(nodeID == "P"))

  # turn into tibble
  df <- as_tibble(df)

  return(df)

}
