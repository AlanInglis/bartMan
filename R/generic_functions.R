#' Get Observations Falling into Each Node
#'
#' This function determines which observations from a given dataset fall into
#' which nodes of a tree, based on a tree structure defined in `treeData`.
#' The treeData object must include `iteration`, `treeNum`, `var`, and `splitValue` columns.
#'
#' @param data A data frame used to build BART model.
#' @param treeData A data frame representing the tree structure, including the
#'        necessary columns `iteration`, `treeNum`, `var`, and `splitValue`.
#'
#' @return A modified version of `treeData` that includes two new columns: `obsNode` and
#'         `noObs`. `obsNode` lists the observations falling into each node, and
#'         `noObs` provides the count of observations for each node.
#'
#' @importFrom dplyr group_by mutate
#'
#' @examples
#'  \donttest{
#' treeDataWithObservations <- getObservations(data = my_data, treeData = my_trees)
#' }
#' @export

getObservations <- function(data, treeData){

  dat <- data

  evalNode <- function(df, x, v) {
    out <- vector("list", length(x))
    stk <- vector("list", sum(is.na(x)))
    pos <- 1L
    stk[[pos]] <- df

    for (i in seq_along(x)) {
      if (!is.na(x[[i]])) {

        subs <- pos + c(0L, 1L)
        stk[subs] <- split(stk[[pos]], stk[[pos]][[x[[i]]]] <= v[[i]])

        names(stk)[subs] <- trimws(paste0(
          names(stk[pos]), ",", x[[i]], c(">", "<="), v[[i]]
        ), "left", ",")

        out[[i]] <- rev(stk[subs])
        pos <- pos + 1L
      } else {

        out[[i]] <- stk[pos]
        stk[[pos]] <- NULL
        pos <- pos - 1L
      }
    }

    return(out)
  }

  dfObs <-  treeData |>
    dplyr::group_by(iteration, treeNum) |>
    dplyr::mutate(obsList = evalNode(dat, var, splitValue))

  obsIndex <- lapply(dfObs$obsList, function(x) {
    lapply(x, row.names)
  })

  whichObs <- lapply(obsIndex, rapply, f = c)
  whichObs <- lapply(whichObs, as.numeric)

  treeData$obsNode <- whichObs

  # get number of observation
  noObser <- NULL
  for(i in 1:nrow(dfObs)){
    noObser[[i]] <- lapply(dfObs$obsList[[i]], dim)
  }

  treeData$noObs <- sapply(noObser, function(y) sum(do.call(rbind, y)[, 1]))

  return(treeData)
}


# -------------------------------------------------------------------------

#' Generate Terminal Node Indicator
#'
#' Adds a boolean `terminal` column to the dataset indicating whether each node is terminal.
#'
#' @param data A data frame containing tree structure information with at least `treeNum`,
#'        `iteration`, and `depth` columns.
#'
#' @return The modified data frame with an additional `terminal` column.
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export


terminalFunction <- function(data){

  message("Generating Terminal Node Indicator:\n")
  # Set up the progress bar
  total_iterations <- nrow(data)
  pb <- txtProgressBar(min = 0, max = total_iterations, style = 3)

  # get terminal column
  data$terminal <- FALSE

  for(i in 1:nrow(data)){
    # progress bar
    setTxtProgressBar(pb, i)

    # Check if last row
    if(i == nrow(data)){
      data$terminal[i] <- TRUE
      break
    }

    current_treeNum <- data$treeNum[i]
    current_iteration <- data$iteration[i]
    current_depth <- data$depth[i]

    next_treeNum <- data$treeNum[i+1]
    next_iteration <- data$iteration[i+1]
    next_depth <- data$depth[i+1]

    # Determine if the current node is terminal
    if(current_treeNum != next_treeNum | current_iteration != next_iteration | next_depth <= current_depth){
      data$terminal[i] <- TRUE
    } else {
      # For nodes within the same tree and iteration, check if the next node is at a higher depth
      if(current_depth < next_depth){
        data$terminal[i] <- FALSE
      }
    }
  }
  # Close the progress bar
  close(pb)
  return(data)
}



# -------------------------------------------------------------------------

#' Generate Child and Parent Node Relationships
#'
#' Populates `childLeft`, `childRight`, and `parent` columns in the dataset to establish
#' parent-child relationships between nodes based on tree structure.
#'
#' @param data A data frame with tree structure, including `iteration`,
#'        `treeNum`, `node`, and `depth` columns, along with a `terminal` indicator.
#'
#' @return The modified data frame with `childLeft`, `childRight`, and `parent` columns added,
#'         detailing the tree's parent-child node relationships.
#'
#' @examples
#'  \dontrun{
#' # Assuming `treeData` has been prepared with terminal node indicators
#' treeData <- getChildren(data = treeData)
#' }
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export


getChildren <- function(data) {
  message("Generating Child/Parent Mappings:\n")

  # Initialize the progress bar
  total_iterations <- nrow(data)
  pb <- txtProgressBar(min = 0, max = total_iterations, style = 3)

  # Prepare the data frame: initialize columns for children and parent node numbers
  data$childLeft <- NA
  data$childRight <- NA
  data$parent <- NA

  for (i in 1:nrow(data)) {
    setTxtProgressBar(pb, i)

    if (!data$terminal[i]) { # Skip terminal nodes
      current_iteration <- data$iteration[i]
      current_treeNum <- data$treeNum[i]
      current_node <- data$node[i]

      # Find immediate children nodes within the same iteration and tree, at the next depth
      children <- which(data$iteration == current_iteration & data$treeNum == current_treeNum &
                          data$depth == data$depth[i] + 1 & data$node > current_node)

      if (length(children) >= 1) {
        # Assign childLeft as the node number of the first child
        data$childLeft[i] <- data$node[children[1]]
        # Set the parent of the left child as the current node
        data$parent[children[1]] <- current_node

        if (length(children) >= 2) {
          # Assign childRight as the node number of the second child, if it exists
          data$childRight[i] <- data$node[children[2]]
          # Set the parent of the right child as the current node
          data$parent[children[2]] <- current_node
        }
      }
    }
  }

  # Close the progress bar
  close(pb)

  return(data)
}

# END



# -------------------------------------------------------------------------

#' Calculate Node Depths in a Tree Data Frame
#'
#' Computes the depth of each node in a given tree data frame, assuming a binary tree structure.
#' Requires the tree data frame to contain a logical column `terminal` indicating terminal nodes.
#'
#' @param tree A data frame representing a tree, must contain a `terminal` column.
#'
#' @return A vector of depths corresponding to each node in the tree.
#' @export


node_depth <- function(tree) {
  stopifnot(is.data.frame(tree)) # Ensure input is a data frame
  stopifnot("terminal" %in% names(tree)) # Ensure 'terminal' column exists
  depths <- rep(0, nrow(tree)) # Initialize depths vector with zeros

  # Recursive function to calculate depth
  recurse <- function(idx, current_depth, depths) {
    if (idx > nrow(tree)) { # Base case: exceed number of rows
      return(list(idx=idx, depths=depths))
    }
    depths[idx] = current_depth # Assign current depth to node
    if (tree$terminal[idx]) { # If node is terminal, move to next node
      return(list(idx=idx+1, depths=depths))
    }
    # Recursive call for left child
    step <- recurse(idx+1, current_depth + 1, depths)
    # Recursive call for right child using updated index and depths
    step <- recurse(step$idx, current_depth + 1, depths=step$depths)
    return(list(idx=step$idx, depths=step$depths))
  }

  # Start recursion from the first node with depth 1
  recurse(1, 1, depths)$depths
}


# END

# -------------------------------------------------------------------------

#' Cluster Trees by Variable
#'
#' Reorders a list of tree structures based on the clustering of variables within each tree.
#'
#' @param tree_list A list of trees, where each tree is expected to have a 'var' column.
#'
#' @return A list of trees reordered based on the clustering of variables.
#'
#' @importFrom dplyr mutate group_by arrange ungroup
#' @importFrom tidyr replace_na
#' @importFrom purrr map
#' @importFrom tibble as_tibble
#'
#' @export

####################
# Cluster Function #
####################

clusterTrees <- function(tree_list) {

  #df <- data
  indIDS <- map(tree_list, function(x) {
    x %>%
      pull(var) %>%
      replace_na("a") %>%
      paste0(collapse = "")
  }) %>%
    unlist(use.names = F) %>%
    as_tibble() %>%
    mutate(ids = 1:n()) %>%
    group_by(value) %>%
    mutate(count = n():1) %>%
    arrange(value)

  ind <- indIDS %>%
    group_by(value) %>%
    mutate(valrank = max(count)) %>%
    ungroup() %>%
    arrange(-valrank, value, -count) %>%
    pull(ids)


  tree_list <- tree_list[ind]

  return(tree_list)
}

# END


# -------------------------------------------------------------------------
#' Sort Trees by Maximum Depth
#'
#' @param tree_list List of `tbl_graph` trees.
#' @return Sorted list of `tbl_graph` trees by decreasing maximum depth.
#' @export
sort_trees_by_depthMax <- function(tree_list) {
  # Extract the maximum depthMax value from each tree
  max_depths <- sapply(tree_list, function(tree) {
    max(as_tibble(tree)$depthMax, na.rm = TRUE)
  })

  # Sort the list of trees by the extracted maximum depthMax values
  sorted_trees <- tree_list[order(max_depths,decreasing = TRUE)]

  return(sorted_trees)
}

# END

# -------------------------------------------------------------------------
#' Determines the stump color for a legend based on its mean value
#'
#' This function is internal and is used to compute the color of a stump
#' for the purpose of legend display, based on the mean value relative to specified limits.
#'
#' @param lims A numeric vector of length 2 specifying the limits within which the mean value falls.
#' @param mean_value The mean value for which the color needs to be determined.
#' @param palette A character vector of colors representing the palette from which the color is selected.
#' @return A character string specifying the color corresponding to the mean value.

get_stump_colour_for_legend <- function(lims, mean_value, palette){

  # Calculate the position of the mean value within  range
  normalized_position <- (mean_value - lims[1]) / (lims[2] - lims[1])

  # Calculate the index in the high-resolution palette that corresponds to the mean value
  index <- round(normalized_position * (length(palette) - 1)) + 1

  # Get the color corresponding to the mean value
  mean_color <- palette[index]

  return(mean_color)
}

# END



# -------------------------------------------------------------------------

#' Update Dummy Variable Names
#'
#' @description This function updates the `var` column in the `structure` component of the `trees` list,
#' replacing dummy variable names derived from factor variables with their original factor variable names.
#'
#' @param trees A list containing at least two components: `data` and `structure`.
#' `data` should be a dataframe, and `structure` a dataframe that includes a `var` column.
#'
#' @return The modified `trees` list with updated `var` column entries in `trees$structure`.
#'
#' @details The function first identifies factor variables in `trees$data`, then checks each entry
#' in `trees$structure$var` for matches with these factor variables. If a match is found, indicating
#' a dummy variable, the entry is replaced with the original factor variable name.
#'
#' @examples
#' \donttest{
#' df_trees <- extractTreeData(model = my_model, data = my_data)
#' combined_trees <- combineDummy(trees = df_trees)
#' }
#'
#' @export

combineDummy <- function(trees) {
  # Identify factor columns in trees$data
  dfOG <- trees$data
  factorColNam <- names(which(!(sapply(dfOG[colnames(dfOG)], is.numeric))))
  factorCols <- which((colnames(dfOG) %in% factorColNam))

  # Create a list of the variables split into their factors
  facLevelsList <- list()
  for (i in 1:length(factorCols)) {
    facLevels <- unique(dfOG[, factorCols[i]])
    facLevelsList[[names(dfOG)[factorCols[i]]]] <- as.character(facLevels)
  }

  # Function to find original factor name for dummy variables
  find_original_factor_name <- function(varName, factorList) {
    for (factorName in names(factorList)) {
      if (grepl(factorName, varName)) {
        return(factorName)
      }
    }
    return(varName) # Return original if no match found
  }

  # Update entries in trees$structure$var
  trees$structure$var <- sapply(trees$structure$var, find_original_factor_name, factorList = facLevelsList)
  names(trees$structure$var) <- NULL
  trees$varName <- names(trees$data)
  return(trees)
}
