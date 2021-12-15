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
#' @importFrom tidyr fill
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr arrange
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr as_tibble
#' @importFrom data.table rowid
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

  # Melting the tree data into useable format
  df <- rrapply::rrapply(nodeData, how = 'melt')
  nCol <- ncol(df)

  df <- df %>%
    pivot_longer(cols = 3:(nCol-1), values_drop_na = TRUE, names_repair = "unique") %>%
    filter(grepl('depth|isLeaf|is_stump|string_location|y_pred|splitValue|splitAttributeM|n_eta', value...5)) %>%
    select(-name) %>%
    mutate(rn = rowid(L1, L2, value...5)) %>%
    pivot_wider(names_from = value...5, values_from = value...3) %>%
    select(-rn) %>% as_tibble()

  # convert to correct types
  df$depth    <- as.numeric(df$depth)
  df$isLeaf   <- as.logical(df$isLeaf)
  df$n_eta    <- as.numeric(df$n_eta)
  df$is_stump <- as.logical(df$is_stump)
  df$string_location <- as.character(df$string_location)
  df$splitAttributeM <- as.numeric(df$splitAttributeM)
  df$splitValue <- as.numeric(df$splitValue)
  df$y_pred <- as.numeric(df$y_pred)


  # match var number to varName
  names(varNames) <- c(0:(length(varNames)-1))
  df$var <- varNames[as.character(df$splitAttributeM)]

  # Add tree number (ignoring iterration)
  df$treeNumID <-cumsum(df$string_location=="P")

  # define node number sequentially
  df <- df %>%
    group_by(treeNumID) %>%
    mutate(node = row_number())

  # get the parent node
   df <- df %>%
     group_by(treeNumID) %>%
     mutate(parentNode = substr(string_location, 0, nchar(string_location)-1))
  #
  # # define parent node of P as NA
   df$parentNode[df$string_location == "P"] <- NA
  #
  # # define parent nodes nodeID with still empty parent node name as "P"
   df$parentNode[df$parentNode==""] <- "P"

  # Match parent node names to node numbers
   df <- df %>%
     group_by(treeNumID) %>%
     mutate(parentNodeNo = match(parentNode, string_location))

  # round values
   df$splitValue <- round(as.numeric(df$splitValue),4)
   df$y_pred <- round(as.numeric(df$y_pred),  4)

  # add value column
  df <-  df %>%
    mutate(value = coalesce(splitValue, y_pred))

  # add label column
  df <- transform(df, label = ifelse(is.na(splitValue), value, paste(var, value, sep = " ≤ ")))

  # add new column defining the 'to', for the nodes 'from-to'
  df <- df %>%
    mutate(to = node)
  df <- transform(df, to = ifelse(is.na(parentNode), NA, to))


  # turn into tibble
  df <- as_tibble(df)

  # rename columns and reorder/remove cols
  names(df) <- c("iteration", "treeNum", "depth", "isLeaf", "noObs", "isStump", "direction",
                 "splitAtt", "splitValue", "leafValue", "var", "treeNumID",
                 "node", "parentNode", "from", "value", "label", "to")
  df <- df %>%
    select(
      "iteration",
      "treeNum",
      "var",
      "isStump",
      "isLeaf",
      "splitValue",
      "depth",
      "noObs",
      "direction",
      "leafValue",
      "node",
      "parentNode",
      "treeNumID",
      "label",
      "value",
      'from',
      "to")


  class(df) <- c("bartMach", "tbl_df", "tbl", "data.frame")

  return(df)

}
