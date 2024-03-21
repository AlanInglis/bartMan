#' extractTreeData
#'
#'
#' @description Creates a list of all tree attributes for a model
#' created by either the BART, dbarts or bartMachine packages.
#'
#' @param model Model created from either the BART, dbarts or bartMachine packages.
#' @param data a data frame used to build the BART model.
#' @return A list containing the extracted and processed tree data, including dataframe of trees and attributes.
#'
#'
#' @importFrom purrr map_df
#' @importFrom dplyr tibble
#' @importFrom dplyr left_join
#' @importFrom dplyr bind_rows
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @importFrom dplyr coalesce
#' @importFrom stats complete.cases
#' @importFrom dplyr filter
#' @importFrom dplyr summarize
#' @importFrom utils read.table
#' @importFrom dplyr rename
#' @importFrom rrapply rrapply
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr as_tibble row_number
#' @importFrom dbarts extract
#' @importFrom rJava .jcall is.jnull
#'
#' @export
# Extract Tree Data Method ------------------------------------------------

extractTreeData <- function(model, data){
  trees <- extractTrees(model, data)

  hideHelper1 <- function(df){
    class(df) <- c("hideHelper1", class(df))
    df
  }
  trees <- hideHelper1(trees)
  return(trees)
}

# END

# Main Method Function -----------------------------------------------------------

extractTrees <- function(model, data) {
  UseMethod("extractTrees")
}

# END

############################
###   Package Methods    ###
############################


# BART package method -----------------------------------------------------

extractTrees.pbart <- function(model, data){
  extractTrees.wbart(model, data)
}

extractTrees.wbart <- function(model, data){
  # variable names:
  varNames <- names(model$varcount.mean)
  # get trees from model
  modelTrees <- model$treedraws$trees

  # extracting tree structure
  trees <- list()
  trees$structure <- utils::read.table(text = modelTrees,
                                       skip = 1,
                                       fill = NA,
                                       col.names = c("node", "var", "splitValue", "leafValue"))
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
  trees$varName <- varNames
  trees$data  <- data

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
  trees$structure$terminal <- is.na(trees$structure$childLeft) & is.na(trees$structure$childRight)
  # Remove leaf values for non-leaves
  trees$structure$leafValue <- ifelse(trees$structure$terminal, trees$structure$leafValue, NA_real_)
  # Remove split values for leaves
  trees$structure$splitValue <- ifelse(trees$structure$terminal, NA_real_, trees$structure$splitValue)
  # Remove var names for leaves
  trees$structure$var <- ifelse(trees$structure$terminal, NA_character_, trees$structure$var)
  # Add a label column
  trees$structure$label <- ifelse(trees$structure$terminal,
                                  as.character(round(trees$structure$leafValue, digits = 2)),
                                  paste(trees$structure$var, " \U2264 ", round(trees$structure$splitValue, digits = 2))
  )
  # Add parent column
  trees$structure$parent <- parent(trees$structure$node)
  # Add value column
  trees$structure <- trees$structure  |>
    dplyr::mutate(value = dplyr::coalesce(splitValue, leafValue))


  # renumber nodes to keep consistent with other packages
  trees$structure <- trees$structure |>
    dplyr::group_by(iteration, treeNum) |>
    dplyr::mutate(node = dplyr::row_number()) |>
    dplyr::ungroup()

  # add depth
  names(trees$structure)[names(trees$structure) == "tier"] <- "depth"

  # max depth
  trees$structure <-  trees$structure |>
    dplyr::group_by(iteration, treeNum) |>
    dplyr::mutate(depthMax = max(depth)) |>
    dplyr::ungroup()

  # get children and parent columns to be consistent with new node ordering
  trees$structure <-  trees$structure |>
    dplyr::group_by(iteration, treeNum, node)
  trees$structure <- getChildren(data = trees$structure)
  trees$structure <-  trees$structure |> dplyr::ungroup()

  cat("Extracting Observation Data...\n")
  dataNew <- as.data.frame(data)
  dat <- BART::bartModelMatrix(dataNew)
  dat <- as.data.frame(dat)

  # get observations
  trees$structure <- getObservations(data = dat, treeData = trees$structure)

  # add is stump column
  trees$structure <-  trees$structure  |>
    dplyr::mutate(isStump = is.na(childLeft) & is.na(childRight) & is.na(parent) & depth == 0)

  # reordering the data and removing unnecessary columns
  trees$structure <- dplyr::select(
    dplyr::group_by(trees$structure, iteration, treeNum),
    var,
    splitValue,
    terminal,
    leafValue,
    iteration,
    treeNum,
    node,
    childLeft,
    childRight,
    parent,
    depth,
    depthMax,
    isStump,
    label,
    value,
    obsNode,
    noObs)|>
    ungroup()

  # add class
  class(trees) <- c("list", "bart", "wbart")

  return(trees)
}

# END

# dbarts package method -----------------------------------------------------


extractTrees.bart <- function(model, data){
  # get all trees
  treesTotal <- model$call$ntree
  iteration  <- model$call$ndpost


  trees <- list()
  trees$structure <- dbarts::extract(model, "trees")

  # add other info
  trees$nMCMC <- as.integer(iteration)
  trees$nTree <- as.integer(treesTotal)
  trees$nVar  <- as.integer(length(colMeans((model$varcount))))
  trees$data  <- data
  trees$varName <- colnames(model$varcount)

  # Get variable names
  varNames <- colnames(model$fit$data@x)

  # set up data frame
  trees$structure$node <- 1:(nrow(trees$structure))
  trees$structure$value <- round(trees$structure$value, 4)
  trees$structure <- transform(trees$structure, terminal = ifelse(var < 0, TRUE, FALSE))
  trees$structure <- transform(trees$structure, leafValue = ifelse(terminal == TRUE, value, NA_integer_))
  trees$structure <- transform(trees$structure, splitValue = ifelse(terminal == FALSE, value, NA_integer_))
  trees$structure <- transform(trees$structure, varName = ifelse(var < 0, NA, var))
  trees$structure$varName <- varNames[trees$structure$varName]

  # set var column
  trees$structure <- trees$structure |>
    group_by(tree, sample) |>
    mutate(node = dplyr::row_number()) |>
    ungroup() |>
    mutate(var = varName) |>
    dplyr::rename(iteration = sample, treeNum = tree) |>
    select( - varName)

  # label
  trees$structure$label <- ifelse(trees$structure$terminal,
                                  as.character(round(trees$structure$leafValue, digits = 2)),
                                  paste(trees$structure$var, " \U2264 ", round(trees$structure$splitValue, digits = 2))
  )

  # reorder columns
  trees$structure <- trees$structure |>
    select(
      var,
      splitValue,
      node,
      terminal,
      leafValue,
      iteration,
      treeNum,
      label,
      value
    )

  # add depth
  depthList <- lapply(split(trees$structure, ~treeNum + iteration),
                   function(x) cbind(x, depth = node_depth(x)-1))

  trees$structure <- dplyr::bind_rows(depthList, .id = "list_id")

   # max depth
  trees$structure <-  trees$structure |>  # GEN
    group_by(iteration, treeNum) |>
    mutate(depthMax = max(depth)) |>
    ungroup()

  # get children and parent columns
  trees$structure <-  trees$structure |>
    group_by(iteration, treeNum, node)
  trees$structure <- getChildren(data = trees$structure)
  trees$structure <-  trees$structure |> ungroup()

  cat("Extracting Observation Data...\n")
  # get observations
  dat <- as.data.frame(model$fit$data@x)
  trees$structure <- getObservations(data = dat, treeData = trees$structure)

  # add is stump column
  trees$structure <-  trees$structure  |>
    mutate(isStump = is.na(childLeft) & is.na(childRight) & is.na(parent) & depth == 0)

  # reordering the data and removing unnecessary columns
  trees$structure <- dplyr::select(
    dplyr::group_by(trees$structure, iteration, treeNum),
    var,
    splitValue,
    terminal,
    leafValue,
    iteration,
    treeNum,
    node,
    childLeft,
    childRight,
    parent,
    depth,
    depthMax,
    isStump,
    label,
    value,
    obsNode,
    noObs)|>
    ungroup()

  # add class

  class(trees) <- c("list", "dbarts")

  return(trees)
}

# END

# bartMachine package method -----------------------------------------------------

extractTrees.bartMachine <- function(model, data){

  # Get variable names
  varNames <- model$training_data_features

  # Get No of iterations after burn in
  iter <- model$num_iterations_after_burn_in

  # function to extract tree data
  nodeData <- vector("list", iter)


  # extract node data
  #  progress bar
  cat("Extracting Node Data:\n")
  pb <- txtProgressBar(min = 0, max = iter, style = 3)

  # Define wrapper for progress bar
  wrapped_function <- function(i) {
    setTxtProgressBar(pb, i)
    extract_raw_node_dataSP(model, g = i, iter = iter)
  }

  # actually extract nodes
  nodeData <- lapply(1:iter, wrapped_function)

  # Close progress bar
  close(pb)


  # Melting the tree data into useable format
  df <- rrapply::rrapply(nodeData, how = 'melt')
  nCol <- ncol(df)

  suppressMessages(
    df <- df |>
      tidyr::pivot_longer(cols = 3:(nCol-1), values_drop_na = TRUE, names_repair = "unique") |>
      filter(grepl('depth|isLeaf|is_stump|string_location|y_pred|splitValue|splitAttributeM', value...5)) |>
      select(-name) |>
      # mutate(rn = rowid(L1, L2, value...5)) |>
      group_by(L1, L2, value...5) |>  mutate(rn = dplyr::row_number()) |>
      ungroup() |>
      tidyr::pivot_wider(names_from = value...5, values_from = value...3) |>
      select(-rn) |>  as_tibble()
  )

  # ®ename columns
  names(df) <- c("iteration", "treeNum", "depth", "terminal", "isStump", "direction",
                 "splitAtt", "splitValue", "leafValue")

  # convert to correct types
  df$iteration <- as.numeric(df$iteration)
  df$treeNum <- as.numeric(df$treeNum)
  df$depth  <- as.numeric(df$depth)
  df$terminal <- as.logical(df$terminal)
  df$isStump <- as.logical(df$isStump)
  df$direction <- as.character(df$direction)
  df$splitAtt <- as.numeric(df$splitAtt)
  df$splitValue <- as.numeric(df$splitValue)
  df$leafValue <- as.numeric(df$leafValue)
  # Remove leaf values for non-terminal nodes
  df$leafValue <- ifelse(df$terminal, df$leafValue, NA_real_)

  # match var number to varName
  names(varNames) <- c(0:(length(varNames)-1))
  df$var <- varNames[as.character(df$splitAtt)]

  # define node number sequentially
  df <- df |>
    group_by(iteration, treeNum) |>
    mutate(node = dplyr::row_number()) |>
    ungroup()


  # create paernt and children columns
  df <- getChildren(data = df)

  # fix potential overflow error with bartmachine
  df$splitValue <- ifelse(is.na(df$splitAtt), NA, df$splitValue)

  # round values
  df$splitValue <- signif(df$splitValue, digits = 4)
  df$leafValue <- signif(df$leafValue, digits = 4)

  # add value column
  df <-  df |>
    mutate(value = coalesce(splitValue, leafValue))

  # add label column
  df$label  <- ifelse(df$terminal,
                      as.character(round(df$leafValue, digits = 2)),
                      paste(df$var, " \U2264 ", round(df$splitValue, digits = 2)))

  # max depth column
  df <- df |>
    group_by(iteration, treeNum) |>
    mutate(depthMax = max(depth)) |>
    ungroup()



  cat("Extracting Observation Data...\n")
  # get which observations
  dat <- model$model_matrix_training_data
  dat <- as.data.frame(dat)
  dat <- dat[,-(length(dat))]

  df <- getObservations(treeData = df,
                        data = dat)



  # reordering the data and removing unnecessary columns
  df <- dplyr::select(
    dplyr::group_by(df, iteration, treeNum),
    var,
    splitValue,
    terminal,
    leafValue,
    iteration,
    treeNum,
    node,
    childLeft,
    childRight,
    parent,
    depth,
    depthMax,
    isStump,
    label,
    value,
    obsNode,
    noObs) |>
    ungroup()

  # turn into list object
  trees <- list()
  trees$structure <- df
  trees$nMCMC <- iter
  trees$nTree <- model$num_trees
  trees$nVar <- model$p
  trees$varName <- model$training_data_features
  trees$data  <- model$X

  # set the class
  class(trees) <- c("bartMach", "list")

  return(trees)
}

# END




# Print Helper Function ---------------------------------------------------


#' print.hideHelper
#' @description This function hides parts from the print out
#' but are still accessible via indexing.
#' @param x A data frame of trees
#' @param ... Extra parameters
#' @export

print.hideHelper1 <- function(x, ...) {
  cat("Tree dataframe:\n")
  print(x$structure)
  cat("Variable names:\n")
  print(x$varName)
  cat("nMCMC:\n")
  print(x$nMCMC)
  cat("nTree:\n")
  print(x$nTree)
  cat("nVar:\n")
  print(x$nVar)
}

# END


# Function to improve bartMachine df generation time ----------------------


# Function to improve bartMachine speeds ----------------------------------

extract_raw_node_dataSP <- function (bart_machine, g = 1, iter)
{

  raw_data_java = rJava::.jcall(bart_machine$java_bart_machine, "[LbartMachine/bartMachineTreeNode;",
                                "extractRawNodeInformation", as.integer(g - 1), simplify = TRUE)

  raw_data <- vector('list', iter)
  raw_data <- lapply(raw_data_java, bMachineNode)
  raw_data
}


# recursively go through java object
bMachineNode <- function (node_java)
{
  if (!requireNamespace("bartMachine", quietly = TRUE)) {
    stop("Package \"bartMachine\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  BAD_FLAG_INT = -2147483647
  BAD_FLAG_DOUBLE = -1.7976931348623157e+308


  node_data = list()
  #node_data = vector("list", 19)
  node_data$java_obj = node_java
  node_data$depth = node_java$depth
  node_data$isLeaf = node_java$isLeaf
  node_data$n_eta = node_java$n_eta
  node_data$is_stump = node_java$isStump()
  node_data$string_location = node_java$stringLocation()


  if (node_java$splitAttributeM == BAD_FLAG_INT) {
    node_data$splitAttributeM = NA
  }
  else {
    node_data$splitAttributeM = node_java$splitAttributeM
  }


  if (node_java$splitValue == BAD_FLAG_DOUBLE) {
    node_data$splitValue = NA
  }
  else {
    node_data$splitValue = node_java$splitValue
  }


  if (node_java$y_pred == BAD_FLAG_DOUBLE) {
    node_data$y_pred = NA
  }
  else {
    node_data$y_pred = node_java$y_pred
  }


  if (!rJava::is.jnull(node_java$left)) {
    node_data$left = bMachineNode(node_java$left)
  }
  else {
    node_data$left = NA
  }
  if (!rJava::is.jnull(node_java$right)) {
    node_data$right = bMachineNode(node_java$right)
  }
  else {
    node_data$right = NA
  }
  node_data
}

