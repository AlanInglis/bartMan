#' extractTreeData
#'
#' @description Creates a list of all tree attributes for a model
#' created by either the BART, dbarts or bartMachine packages.
#'
#' @param model Model created from either the BART, dbarts or bartMachine packages.
#'
#' @return A list of every tree and its attributes.
#'
#'
#' @importFrom readr read_table
#' @importFrom readr cols
#' @importFrom readr col_integer
#' @importFrom readr col_double
#' @importFrom purrr map_df
#' @importFrom dplyr tibble
#' @importFrom dplyr left_join
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @importFrom dplyr coalesce
#' @importFrom stats complete.cases
#' @importFrom dplyr one_of
#' @importFrom dplyr group_split
#' @importFrom dplyr filter
#' @importFrom dplyr summarize
#'
#' @importFrom dplyr rename
#'
#' @importFrom bartMachine extract_raw_node_data
#' @importFrom rrapply rrapply
#' @importFrom purrr transpose
#' @importFrom purrr flatten
#' @importFrom stringr str_c
#' @importFrom tidyr fill
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr arrange
#' @importFrom dplyr as_tibble
#' @importFrom data.table rowid
#'
#' @export



extractTreeData <- function(model){
  trees <- extractTrees(model)
  return(trees)
}


# -------------------------------------------------------------------------

# Main function:
extractTrees <- function(model) {
  UseMethod("extractTrees")
}


# BART --------------------------------------------------------------------

extractTrees.wbart <- function(model){

  # variable names:
  varNames <- names(model$varcount.mean)
  # get trees from model
  modelTrees <- model$treedraws$trees

  # extracting tree structure
  trees <- list()
  trees$structure <- suppressWarnings(
    readr::read_table(
      file = modelTrees,
      col_names = c("node", "var", "splitValue", "leafValue"),
      col_types =
        readr::cols(
          node = readr::col_integer(),
          var = readr::col_integer(),
          splitValue = readr::col_integer(),
          leaf = readr::col_double()
        ),
      skip = 1,
      na = c("")
    )
  )

  # Adding in columns
  trees$structure$var <- varNames[trees$structure$var + 1] # as vars are indexed at 0
  trees$structure$splitID <- trees$structure$splitValue + 1
  trees$structure$tier <- as.integer(floor(log2(trees$structure$node)))

  # getting split points
  splitPoints <- purrr::map_df(
    .x = model$treedraws$cutpoints,
    .f = ~ dplyr::tibble(splitValue = ., splitID = 1:length(.)),
    .id = "var"
  )

  # adding split points into tree structure
  trees$structure <- dplyr::left_join(
    dplyr::select(trees$structure, -splitValue),
    splitPoints,
    by = c("var", "splitID")
  )

  # Add in model fit info
  modelInfo <- unlist(strsplit(modelTrees, " "))[1:3]
  modelInfo <- gsub("(^\\d+)([\a-zA-Z0-9]*)", "\\1", modelInfo)

  trees$nMCMC <- as.integer(modelInfo[1])
  trees$nTree <- as.integer(modelInfo[2])
  trees$nVar  <- as.integer(modelInfo[3])

  trees$structure$uniqueTreeID <- cumsum(is.na(trees$structure$var) & is.na(trees$structure$splitValue) & is.na(trees$structure$leafValue))
  trees$structure$iteration <- ((trees$structure$uniqueTreeID - 1) %/% trees$nTree) + 1
  trees$structure$treeNum <- ((trees$structure$uniqueTreeID - 1) %% trees$nTree) + 1
  trees$structure$uniqueTreeID <- NULL

  # remove information about tree groups (i.e., rows with missing data)
  trees$structure$missingData <- complete.cases(trees$structure)
  missingIndex <- which(trees$structure$missingData == F)
  if (length(missingIndex == 0)) {
    trees$structure <- trees$structure[-missingIndex, ]
  }
  trees$structure$missingData <- NULL

  # Functions to get the left and right children nodes
  # and the parent nodes
  childLeft <- function(nodes) {
    childL <- nodes * 2
    childL[!childL %in% nodes] <- NA_integer_

    return(childL)
  }

  childRight <- function(nodes) {
    childR <- nodes * 2 + 1
    childR[!childR %in% nodes] <- NA_integer_

    return(childR)
  }

  parent <- function(nodes) {
    parents <- nodes %/% 2
    parents[parents == 0] <- NA_integer_

    return(parents)
  }

  trees$structure <- dplyr::group_by(trees$structure, iteration, treeNum)
  trees$structure <- dplyr::mutate(
    trees$structure,
    childLeft = childLeft(node),
    childRight = childRight(node)
  )

  trees$structure <- dplyr::ungroup(trees$structure)


  # Add is leaf column
  trees$structure$isLeaf <- is.na(trees$structure$childLeft) & is.na(trees$structure$childRight)
  # Remove leaf values for non-leaves
  trees$structure$leafValue <- ifelse(trees$structure$isLeaf, trees$structure$leafValue, NA_real_)
  # Remove split values for leaves
  trees$structure$splitValue <- ifelse(trees$structure$isLeaf, NA_real_, trees$structure$splitValue)
  # Remove var names for leaves
  trees$structure$var <- ifelse(trees$structure$isLeaf, NA_character_, trees$structure$var)
  # Add a label column
  trees$structure$label <- ifelse(trees$structure$isLeaf,
                                  as.character(round(trees$structure$leafValue, digits = 2)),
                                  paste(trees$structure$var, " ≤ ", round(trees$structure$splitValue, digits = 2))
  )
  # Add parent column
  trees$structure$parent <- parent(trees$structure$node)
  # Add value column
  trees$structure <- trees$structure %>%
    dplyr::mutate(value = dplyr::coalesce(splitValue, leafValue))

  # reordering the data and removing unnecessary columns
  trees$structure <- dplyr::select(
    dplyr::group_by(trees$structure, iteration, treeNum),
    var,
    splitValue,
    node,
    isLeaf,
    leafValue,
    childLeft,
    childRight,
    parent,
    iteration,
    treeNum,
    label,
    value,
    -splitID,
    -tier
  )

  trees$varName <- colnames(model$varcount)

  # add tree depth
  treeDepth <- function(x) {
    vals <- !is.na(x)
    l1_vals <- !is.na(lag(x))
    l2_vals <- !is.na(lag(x, 2L))
    vals & (l1_vals | (l2_vals & cumsum(vals) %% 2L == 0L))
  }

  tDepth <- trees$structure %>%
    group_by(iteration, treeNum) %>%
    summarize(Depth = sum(treeDepth(var))+1)

  lengthDepth <- trees$structure %>%
    group_by(iteration, treeNum) %>%
    summarise(depth = n())

  trees$structure$depth <- rep(tDepth$Depth, times = lengthDepth$depth)



  # add class
  class(trees) <- c("list", "bart")
  return(trees)
}



# dbarts ------------------------------------------------------------------

extractTrees.bart <- function(model){
  # get all trees
  treesTotal <- model$call$ntree
  iteration  <- model$call$ndpost

  trees <- list()
  trees$structure <- model$fit$getTrees(treeNums = 1:treesTotal, sampleNums = 1:iteration)

  # add other info
  trees$nMCMC <- as.integer(iteration)
  trees$nTree <- as.integer(treesTotal)
  trees$nVar  <- as.integer(length(model$varcount))

  # Get variable names
  varNames <- colnames(model$fit$data@x)

  # set up data frame
  trees$structure$node <- c(1:(nrow(trees$structure)))
  trees$structure$value <- round(trees$structure$value, 4)
  trees$structure <- transform(trees$structure, isLeaf = ifelse(var < 0, TRUE, FALSE))
  trees$structure <- transform(trees$structure, leafValue = ifelse(isLeaf == TRUE, value, NA_integer_))
  trees$structure <- transform(trees$structure, splitValue = ifelse(isLeaf == FALSE, value, NA_integer_))
  trees$structure <- transform(trees$structure, varName = ifelse(var < 0, NA, var))
  trees$structure$varName <- varNames[trees$structure$varName]
  trees$structure <- transform(trees$structure, label = ifelse(is.na(varName), value, paste(varName, value, sep = " ≤ ")))
  trees$structure <-  trees$structure %>%
    mutate(value = coalesce(splitValue, leafValue))


  trees$structure <- trees$structure %>%
    group_by(tree, sample) %>%
    mutate(node = row_number()) %>%
    ungroup() %>%
    mutate(var = varName) %>%
    rename(iteration = sample, treeNum = tree) %>%
    select( - varName)

  # reorder columns
  trees$structure <- trees$structure %>%
    select(
      var,
      splitValue,
      node,
      isLeaf,
      leafValue,
      iteration,
      treeNum,
      label,
      value,
      n
    )

  trees$varName <- colnames(model$varcount)

  # add tree depth
  treeDepth <- function(x) {
    vals <- !is.na(x)
    l1_vals <- !is.na(lag(x))
    l2_vals <- !is.na(lag(x, 2L))
    vals & (l1_vals | (l2_vals & cumsum(vals) %% 2L == 0L))
  }

  tDepth <- trees$structure %>%
    group_by(iteration, treeNum) %>%
    summarize(Depth = sum(treeDepth(var))+1)

  lengthDepth <- trees$structure %>%
    group_by(iteration, treeNum) %>%
    summarise(depth = n())

  trees$structure$depth <- rep(tDepth$Depth, times = lengthDepth$depth)

  # add class
  class(trees) <- c("list", "dbarts")
  return(trees)
}


# bartMachine -------------------------------------------------------------

extractTrees.bartMachine <- function(model){
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

  suppressMessages(
  df <- df %>%
    pivot_longer(cols = 3:(nCol-1), values_drop_na = TRUE, names_repair = "unique") %>%
    filter(grepl('depth|isLeaf|is_stump|string_location|y_pred|splitValue|splitAttributeM|n_eta', value...5)) %>%
    select(-name) %>%
    mutate(rn = rowid(L1, L2, value...5)) %>%
    pivot_wider(names_from = value...5, values_from = value...3) %>%
    select(-rn) %>% as_tibble()
  )

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

  # Add tree number (ignoring iteration)
  df$treeNumID <-cumsum(df$string_location=="P")

  # define node number sequentially
  df <- df %>%
    group_by(treeNumID) %>%
    mutate(node = row_number())

  # get the parent node
  df <- df %>%
    group_by(treeNumID) %>%
    mutate(parentNode = substr(string_location, 0, nchar(string_location)-1))

  # # define parent node of P as NA
  df$parentNode[df$string_location == "P"] <- NA

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

  df$iteration <- as.numeric(df$iteration)
  df$treeNum <- as.numeric(df$treeNum)

  trees <- list()
  trees$structure <- df
  trees$MCMC <- iter
  trees$nTree <- model$num_trees
  trees$nVar <- model$p
  trees$varName <- colnames(model$X)


  class(trees) <- c("bartMach", "list")

  return(trees)
}