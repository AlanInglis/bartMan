#' plotAllTrees
#'
#' @description Plots all the trees.
#'
#' @param treeData A list of tree attributes created by exctractTreeData function.
#' @param iter The selected iteration
#' @param treeNo The selected tree number.
#' @param cluster LOGICAL. If TRUE, then cluster by tree structures.
#' @param sizeNode Whether to size node width by the number of observations that fall into that node.
#' @param pal A palette to colour the terminal nodes.
#' @param fillBy Which parameter to colour the terminal nodes. Either 'response' or 'mu'.
#' @param removeStump LOGICAL. If TRUE, then stumps are removed from plot. If False, stumps
#' remain in plot and are coloured grey.
#' @param selectedVars A vector of selected variables to display. Either a character vector of names
#' or the variables column number.
#' @param combineFact If a variable is a factor in a data frame, when building the BART model it is replaced with dummies.
#' Note that q dummies are created if q>2 and one dummy is created if q=2, where q is the number of levels of the factor.
#' If combineFact = TRUE, then the dummy factors are combined into their original factor.
#' @return A list containing vectors of the indices of observations from leaf nodes.
#'
#' @import ggraph
#' @importFrom igraph gsize
#' @importFrom cowplot get_legend
#' @importFrom cowplot plot_grid
#' @importFrom scales hue_pal
#' @importFrom dplyr select
#' @importFrom dplyr one_of
#' @importFrom dplyr mutate
#' @importFrom dplyr group_split
#' @importFrom dplyr filter
#' @importFrom purrr map
#' @importFrom purrr map2
#' @importFrom tidygraph tbl_graph
#' @importFrom dplyr transmute
#' @importFrom dplyr row_number
#' @importFrom dplyr pull
#' @importFrom dplyr n
#' @importFrom stats setNames
#' @importFrom rlang :=
#' @importFrom grDevices colorRampPalette
#'
#' @export

plotAllTrees <- function(treeData,
                         iter = NULL,
                         treeNo = NULL,
                         cluster = NULL,
                         sizeNode = TRUE,
                         pal = rev(colorRampPalette(c('steelblue', '#f7fcfd', 'orange'))(5)),
                         fillBy = NULL,
                         selectedVars = NULL,
                         removeStump = FALSE,
                         combineFact = FALSE
                         ) {

  # ONLY FOR BART ATM:
  if(combineFact){
    treeData <-  cFactTrees(treeData)
  }



  if(length(selectedVars) > length(treeData$varName)){
    message("SelectedVars is longer than number of available variables. Selecting all variables")
    selectedVars <- c(1:length(treeData$varName))
  }

  allTrees <- plotAll(treeData, iter = iter, treeNo = treeNo, cluster = cluster)

  suppressWarnings(
    p <- plotAllTreesPlotFn(allTrees,
                            sizeNode = sizeNode,
                            pal = pal,
                            fillBy = fillBy,
                            name = treeData$varName,
                            selectedVars = selectedVars,
                            removeStump = removeStump,
                            combineFact = combineFact,
                            treeData = treeData
                            )

  )
  return(p)
}



# -------------------------------------------------------------------------

# Main plot function:
plotAll <- function(treeData, iter = NULL, treeNo = NULL, cluster = NULL) {
  UseMethod("plotAll")
}



# BART --------------------------------------------------------------------


#' @export
plotAll.bart <- function(treeData, iter = NULL, treeNo = NULL, cluster = NULL) {
  maxIter <- treeData$nMCMC

  # if (is.null(iter) & is.null(treeNo)) {
  #   df <- treeData$structure %>%
  #     filter(iteration == maxIter)
  # } else if (is.null(iter) & !is.null(treeNo)) {
  #   df <- treeData$structure %>%
  #     filter(treeNum == treeNo)
  # } else if (!is.null(iter) & is.null(treeNo)) {
  #   df <- treeData$structure %>%
  #     filter(iteration == iter)
  # } else {
  #   df <- treeData$structure %>%
  #     filter(iteration == iter, treeNum == treeNo)
  # }

  if (is.null(iter) & is.null(treeNo)) {
    df <- treeData$structure
    message("Both iter and treeNo are NULL. Recommend filtering trees")
  } else if (is.null(iter) & !is.null(treeNo)) {
    df <- treeData$structure %>%
      filter(treeNum == treeNo)
  } else if (!is.null(iter) & is.null(treeNo)) {
    df <- treeData$structure %>%
      filter(iteration == iter)
  } else {
    df <- treeData$structure %>%
      filter(iteration == iter, treeNum == treeNo)
  }


  # add mean response per node:
  respNode <- df$obsNode #apply(df, 1, function(x) {y[x$obsNode]})
  respNode <- lapply(respNode, mean)
  df$respNode <- unlist(respNode)

  if (!is.null(cluster)) {
    if (cluster == "depth") {
      df <- df[with(df, order(-depthMax)), ]
    }
  }

  # Which columns to display
  keeps <- c("var", "node", "parent", "iteration", "treeNum", "label", "value", "depthMax", 'noObs', 'respNode', 'obsNode')

  res <- dplyr::select(
    df,
    dplyr::one_of(keeps)
  )

  # Create edge and node list
  res <- dplyr::mutate(res,
                       newNode = seq_along(node),
                       newParent = newNode[match(parent, node)],
                       node = newNode,
                       parent = newParent
  )
  res <- dplyr::select(res, -newNode, -newParent)

  nodeList <- dplyr::group_split(dplyr::select(res, -parent), .keep = TRUE)

  edgeList <- purrr::map(
    dplyr::group_split(dplyr::select(
      res,
      iteration,
      treeNum,
      parent,
      node
    ), .keep = FALSE),
    ~ dplyr::filter(., !is.na(parent))
  )

  # Turn into data structure for tidy graph manipulation
  tblgList <- purrr::map2(
    .x = nodeList,
    .y = edgeList,
    .f = ~ tidygraph::tbl_graph(
      nodes = .x,
      edges = .y,
      directed = TRUE
    )
  )


  if (!is.null(cluster)) {
    if (cluster == "var") {
      tblgList <- clusterTrees(tblgList)
    }
  }

  return(tblgList)
}


# dbarts ------------------------------------------------------------------
# -------------------------------------------------------------------------
#' @export
plotAll.dbarts <- function(treeData, iter = NULL, treeNo = NULL, cluster = NULL) {

  maxIter <- treeData$nMCMC

  if (is.null(iter) & is.null(treeNo)) {
    treeData$structure <- treeData$structure
    message("Both iter and treeNo are NULL. Recommend filtering trees")
  } else if (is.null(iter) & !is.null(treeNo)) {
    treeData$structure <- treeData$structure %>%
      filter(treeNum == treeNo)
  } else if (!is.null(iter) & is.null(treeNo)) {
    treeData$structure <- treeData$structure %>%
      filter(iteration == iter)
  } else {
    treeData$structure <- treeData$structure %>%
      filter(iteration == iter, treeNum == treeNo)
  }

  noObservations <- max(treeData$structure$noObs)

  treeData$structure <- treeData$structure %>%
    group_by(iteration, treeNum)

  # -------------------------------------------------------------------------

  # add mean response per node:
  respNode <- treeData$structure$noObs #apply(treeData$structure, 1, function(x) {y[x$obsNode]})
  respNode <- lapply(respNode, mean)
  treeData$structure$respNode <- unlist(respNode)


  # cluster trees by depth
  if (!is.null(cluster)) {
    if (cluster == "depth") {
      treeData$structure <-  treeData$structure[with(treeData$structure, order(-depthMax)), ]
    }
  }

  # Which columns to display
  keeps <- c("var", "node", "isLeaf", "iteration", "treeNum", "label", "noObs", "value", "depthMax", "respNode")

  treeData$structure <- dplyr::select(
    treeData$structure,
    dplyr::one_of(keeps)
  )


  treeData$structure <- transform(treeData$structure, varValue = ifelse(is.na(var), -1, var))
  treeData$structure <- treeData$structure %>%
    tibble()


  # treesSplit <- treeData$structure %>%
  #   group_split(cumsum(noObs == noObservations))
  treesSplit <- treeData$structure %>%
    ungroup() %>%
    group_split(cumsum(noObs == noObservations))

  # add the depth of the tree
  treeDepth <- function(trees) {
    if (trees$varValue[1] == -1) {
      return(c(depth = 1, size = 1))
    }

    left <- treeDepth(trees[-1, , drop = FALSE])
    right <- treeDepth(trees[seq.int(2 + left[["size"]], nrow(trees)), , drop = FALSE])

    depthSize <- c(
      depth = 1 + max(left[["depth"]], right[["depth"]]),
      size = 1 + left[["size"]] + right[["size"]]
    )
    depthSize
  }

  # get every tree max depth and add to trees list
  dS <- as.vector(sapply(treesSplit, treeDepth)[1, ])
  treesSplit <- mapply(cbind, treesSplit, "depthNEW" = dS, SIMPLIFY = F)



  # Get Nodes and Edges -----------------------------------------------------

  # Get all the edges
  getEdges <- function(trees) {
    maxDepth <- trees$depthNEW[1]

    if (maxDepth > 2) {
      # Set the edge value:
      # this gets the 1st group of edges
      edgeSet1 <- do.call(rbind, lapply(
        split(trees, ~ cumsum(!isLeaf)),
        function(x) {
          with(x, expand.grid(from = node[!isLeaf], to = node[isLeaf]))
        }
      ))

      # get second group of edges
      edgeSet2 <- setNames(rev(data.frame(stats::embed(trees$node[!trees$isLeaf], 2))), c("from", "to"))

      # bind them together
      edges <- rbind(edgeSet1, edgeSet2)

      # If any number in the to column appears more than twice then
      # replace with n-(number of times it appears-2)
      newFrom <- edges %>%
        group_by(from) %>%
        transmute(from := from - c(rep(0, 2), row_number())[row_number()]) %>%
        ungroup()

      # add to edges df
      edges$from <- pull(newFrom)

      # Change to right side of tree
      backToRoot <- edges %>%
        group_by(from) %>%
        filter(n() == 1)

      backToRoot$from <- 1
      edges <- rbind(edges, backToRoot)

      # remove duplicated rows and single row entries
      edges <- edges[!duplicated(edges), ]
      dups <- edges$from[duplicated(edges$from)]
      edges$unique <- edges$from %in% dups
      edgeIndex <- which(edges$unique == F)
      if (length(edgeIndex == 0)) {
        edges <- edges[-edgeIndex, ]
      }
      drops <- c("unique")
      edges <- edges[, !(names(edges) %in% drops)]

      # reorder
      "edges" <- `row.names<-`(edges[with(edges, order(from, to)), ], NULL)
      edges <- edges[order(edges$to), ]
    } else {
      edges <- do.call(rbind, lapply(
        split(trees, ~ cumsum(!isLeaf)),
        function(x) {
          with(x, expand.grid(from = node[!isLeaf], to = node[isLeaf]))
        }
      ))
    }
  }

  allEdges <- sapply(treesSplit, getEdges)

  #  for single type of tree
  needLe <- paste0("\\b", paste(c(1, 1, 1, 4, 4, 6, 6, 1), collapse = ","), "\\b")
  changeTo <- paste(c(1, 1, 3, 4, 4, 6, 6, 3), collapse = ",")
  allEdges <- lapply(
    strsplit(sapply(allEdges, function(L) gsub(needLe, changeTo, paste(L, collapse = ","))), ","),
    as.integer)
  allEdges <- matrix(allEdges, nrow = 2)
  rownames(allEdges) <- c('from', 'to')


  # remove unnecessary columns
  treesSplit <- lapply(treesSplit, function(x) {
    x["varValue"] <- x["isLeaf"] <- NULL
    x
  })


  # Turn into a table graph object
  eachTree <- list()
  for (i in 1:length(treesSplit)) {
    eachTree[[i]] <- tidygraph::tbl_graph(treesSplit[[i]], as.data.frame(allEdges[, i]))
  }

  if (!is.null(cluster)) {
    if (cluster == "var") {
      eachTree <- clusterTrees(eachTree)
    }
  }

  return(eachTree)
}



# bartMachine -------------------------------------------------------------
# -------------------------------------------------------------------------
#' @export
plotAll.bartMach <- function(treeData, iter = NULL, treeNo = NULL, cluster = NULL) {

  df <- treeData$structure
  maxIter <- treeData$nMCMC
  noObservations <- treeData$structure$noObs[1]

  if (is.null(iter) & is.null(treeNo)) {
    df <- df
    message("Both iter and treeNo are NULL. Recommend filtering trees")
    #%>%
     # filter(iteration == maxIter)
  } else if (is.null(iter) & !is.null(treeNo)) {
    df <- df %>%
      filter(treeNum == treeNo)
  } else if (!is.null(iter) & is.null(treeNo)) {
    df <- df %>%
      filter(iteration == iter)
  } else {
    df <- df %>%
      filter(iteration == iter, treeNum == treeNo)
  }

  df <- df %>%
    group_by(iteration, treeNum)

  # -------------------------------------------------------------------------

  # add mean response per node:
  respNode <- df$obsNode #apply(df, 1, function(x) {y[x$obsNode]})
  respNode <- lapply(respNode, mean)
  df$respNode <- unlist(respNode)



  # cluster trees
  suppressWarnings(
  if (!is.null(cluster)) {
    if (cluster == "depth") {
      df <- df[with(df, order(-depthMax)), ]

      # split the dataframe into a list of dfs, one for each tree
      list_edges <- df %>%
        ungroup() %>%
        group_split(cumsum(noObs == noObservations))
    }
  }
  )

  suppressWarnings(
  if(cluster == 'var' || is.null(cluster)){
    # split the dataframe into a list of dfs, one for each tree
    list_edges <- df %>%
      group_split(cumsum(noObs == noObservations))
  }
)


  # remove unnecessary columns
  treesSplit <- lapply(list_edges, function(x) {
    x["isStump"] <- x["to"] <- x["from"] <- x["node"] <- x["parentNode"] <- x["treeNumID"] <- NULL
    x
  })

  # create dataframe of edges
  dfOfEdges <- lapply(list_edges, function(df_tree) {
    res <- data.frame(
      from = df_tree$from,
      to = df_tree$to
    )
    # delete NAs from result
    res <- stats::na.omit(res)
    return(res)
  })


  # Turn into a table graph object
  eachTree <- list()
  for (i in 1:length(treesSplit)) {
    eachTree[[i]] <- tidygraph::tbl_graph(treesSplit[[i]], dfOfEdges[[i]])
  }

  if (!is.null(cluster)) {
    if (cluster == "var") {
      eachTree <- clusterTrees(eachTree)
    }
  }

  return(eachTree)
}


# -------------------------------------------------------------------------
# cluster function by variable --------------------------------------------
# -------------------------------------------------------------------------

clusterTrees <- function(treeList) {

  #df <- data
    indIDS <- map(treeList, function(x) {
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


    treeList <- treeList[ind]

  return(treeList)
}

# -------------------------------------------------------------------------

#' plotAllTreesplotFn
#'
#' @description This function is used to creat the plots.
#'
#' @param treeList A list of trees created by treeList function.
#' @param sizeNode Whether to size node width by the number of observations that fall into that node.
#' @param pal A palette to colour the terminal nodes.
#' @param selectedVars A vector of selected variables to display.
#' @param fillBy Which parameter to colour the terminal nodes. Either 'response' or 'mu'.
#' @param removeStump LOGICAL. If TRUE, then stumps are removed from plot. If False, stumps
#' remain in plot and are coloured grey.
#' @param combineFact If a variable is a factor in a data frame, when building the BART model it is replaced with dummies.
#' Note that q dummies are created if q>2 and one dummy is created if q=2, where q is the number of levels of the factor.
#' If combineFact = TRUE, then the dummy factors are combined into their original factor.
#' @param name Variable names
#' @param treeData A list of tree attributes created by exctractTreeData function.
#'
#'
#' @return A list containing vectors of the indices of observations from leaf nodes.
#'
#'

plotAllTreesPlotFn <- function(treeList,
                               sizeNode = TRUE,
                               pal =  rev(colorRampPalette(c('steelblue', '#f7fcfd', 'orange'))(5)),
                               fillBy = NULL,
                               selectedVars = NULL,
                               removeStump = FALSE,
                               combineFact = FALSE,
                               treeData,
                               name) {


  # remove stumps
  if(removeStump){
    treeList <- Filter(function(x) igraph::gsize(x) > 0, treeList)
    stumpIdx <- NULL
  }else{
    # get the stump index
    whichStump = NULL
    for(i in 1:length(treeList)){
      whichStump[[i]] <-  which(igraph::gsize(treeList[[i]]) == 0)
    }
    stumpIdx <- which(whichStump == 1)

    if(length(stumpIdx >=1)){
    # create new tree list
    newTrees <- treeList[stumpIdx]

    # create df of tree stumps
    newTreesDF <- NULL
    for(i in 1:length(newTrees)){
      newTreesDF[[i]] <- newTrees[[i]] %>%
        activate(nodes) %>%
        data.frame()
      if(is.null(fillBy)){
        newTreesDF[[i]]$var <- "Stump/Leaf"
      }else{
        newTreesDF[[i]]$var <- "Stump"
      }
    }


    # create edge data for stumps
    newDF_Nodes <- newDF_Edges <- NULL
    for(i in 1:length(newTreesDF)){
      newDF_Nodes[[i]] <- rbind(newTreesDF[[i]], newTreesDF[[i]][rep(1), ])
      newDF_Edges[[i]] <- data.frame(from = c(1,1), to = c(1,1))
    }

    # turn into tidygraph trees
    newTree <- NULL
    for(i in 1:length(newDF_Nodes)){
      newTree[[i]] <- tbl_graph(nodes = newDF_Nodes[[i]], edges = newDF_Edges[[i]])
    }
    # replace stumps with new stumps
    treeList[stumpIdx] <- newTree
    }
  }


  # get limits
  if(!is.null(fillBy)){
  if(fillBy == 'response'){
    lims <- range(unlist(lapply(treeList, . %>% activate(nodes) %>% pull(respNode))))
    lims <- pretty(c(lims[1], lims[2]))
    lims <- c(min(lims), max(lims))
    nam <- 'Mean \nResponse'
  } else if(fillBy == "mu"){
    lims <- range(unlist(lapply(treeList, . %>% activate(nodes) %>% filter(is.na(var) | var == "Stump") %>% pull(value))))
    lims <- c(-max(abs(lims)), max(abs(lims)))
    nam <- 'Mu'
  }
  }else{
    nam <- 'Variable'
    }


  # set node colours
  if(!is.null(selectedVars)){

    if(is.numeric(selectedVars)){
      nodeNamesImp <- name[selectedVars]
      nodeNamesOthers <- name[-selectedVars]
    }else if(is.character(selectedVars)){
      nodeNamesImp <- name[name %in% selectedVars]
      nodeNamesOthers <- name[!(name %in% selectedVars)]
    }

    if(combineFact){
      nameChange <- cFactTreesSelVars(treeData)
      for (i in 1:length(nameChange$dfnew)) {
        nodeNamesImp[which(nodeNamesImp %in% nameChange$dfnew[[i]])] <- nameChange$factorColNam[i]
        nodeNamesOthers[which(nodeNamesOthers %in% nameChange$dfnew[[i]])] <- nameChange$factorColNam[i]
      }
      nodeNamesImp <- unique(nodeNamesImp)
      nodeNamesOthers <- unique(nodeNamesOthers)
    }


    nodeColorsImp <- setNames(scales::hue_pal(c(0, 360) + 15, 100, 64, 0, 1)(length(nodeNamesImp)), nodeNamesImp)
    nodeColorsOth <- setNames(rep('#e6e6e6', length(nodeNamesOthers)), nodeNamesOthers)
    nodecolors <- c(nodeColorsImp, nodeColorsOth)

    namedOthers <- setNames('#e6e6e6',  "Others")
    legColours <- c(nodeColorsImp, namedOthers)
    # dfLegend <- reshape::melt(legColours) %>%
    #   tibble::rownames_to_column(var = 'varName')

    dfLegend <- utils::stack(legColours)
    colnames(dfLegend) <- c('value', 'varName')
    dfLegend$varName <- as.character(dfLegend$varName)

    dfLegend$val <- rep(1, times = length(dfLegend$varName))
    dfLegend$varName <- factor(dfLegend$varName, levels = names(legColours))
    pLeg <- ggplot(dfLegend, aes(x = varName, y = val, fill = varName)) +
      geom_bar(stat = 'identity') +
      scale_fill_manual(values = dfLegend$value, name = 'Variable')
    if(is.null(fillBy)){
      if(length(stumpIdx) >=  1){
        nodecolors[["Stump/Leaf"]] <- '#e6e6e6'
      }
    }else{
      if(length(stumpIdx) >=  1){
        nodecolors[["Stump"]] <- '#e6e6e6'
      }
    }
  }else{
    nodeNames <- unique(stats::na.omit(unlist(lapply(treeList, . %>% activate(nodes) %>% pull(var)))))
    nodeNames <- sort(nodeNames)
    nodecolors <- setNames(scales::hue_pal(c(0, 360) + 15, 100, 64, 0, 1)(length(nodeNames)), nodeNames)

    # colour stumps grey
    if(length(stumpIdx) >=  1){
      if(is.null(fillBy)){
        nodecolors[["Stump/Leaf"]] <- '#808080'
      }else if(fillBy == 'response'){
        nodecolors[["Stump"]] <- pal[ceiling(length(pal)/2)]
      }else if(fillBy == 'mu'){
        # stumpVal <- lapply(treeList, . %>% activate(nodes) %>% filter(is.na(var) | var == "Stump") %>% pull(value))
        # sv  <- as.data.frame(stumpVal[stumpIdx])
        # cdfLims <- ecdf(lims)
        # cdfVal <- cdfLims(sv[1,1])
        # nodecolors[["Stump"]] <-  pal[abs(length(pal) * cdfVal)]
        nodecolors[["Stump"]] <- pal[ceiling(length(pal)/2)]
      }
    }
  }

  suppressMessages(
  allPlots <- lapply(treeList,
                     plotFun,
                     n = length(treeList),
                     color = nodecolors,
                     sizeNode = sizeNode,
                     pal = pal,
                     range = lims,
                     name = nam,
                     fill = fillBy)
  )


  # get legend
  if(!is.null(selectedVars) & !is.null(fillBy)){
    themeMargin <- theme(legend.box.margin = margin(100, 15, 80, 20))
    legend1 <- cowplot::get_legend(pLeg + themeMargin)
    legend <- cowplot::get_legend(allPlots[[1]] + themeMargin)
    legend$grobs[[2]] <- legend1
  }else if(!is.null(selectedVars) & is.null(fillBy)){
    themeMargin <- theme(legend.box.margin = margin(100, 15, 80, 20))
    legend <- cowplot::get_legend(pLeg + themeMargin)
  }else{
    ggdf <- data.frame(x = names(nodecolors), y = c(1:length(nodecolors)))
    ggLegend <- ggplot(ggdf, aes(x=x, y=y))+
                 geom_point(aes(color = nodecolors), shape = 15, size = 10) +
                 scale_color_manual(values = unname(nodecolors),
                                    label = names(nodecolors),
                                    name = 'Variable') +
                 theme_void() +
                 theme(legend.position = 'right')
    legend <- cowplot::get_legend(ggLegend)
    #legend <- cowplot::get_legend(allPlots[[1]])
  }


  # remove legends from individual plots
  allPlots <- lapply(allPlots, function(x) x + theme(legend.position = "none"))
  if(removeStump == FALSE){
    for(i in stumpIdx){
      allPlots[[i]]$data <- allPlots[[i]]$data[-2, ]
    }
  }
  # n <- length(allPlots)
  # nRow <- floor(sqrt(n))
  # allTreesPlot <- arrangeGrob(grobs=allPlots, nrow=nRow)
  # cowplot::plot_grid(allTreesPlot, legend, rel_widths = c(0.9, 0.13), ncol = 2)

  treesPlot <- cowplot::plot_grid(plotlist = allPlots)
  legendPlot <- cowplot::plot_grid(legend)
  cowplot::plot_grid(treesPlot, legendPlot,  rel_widths = c(0.5, 0.13))


}



# Plot function -----------------------------------------------------------

plotFun <- function(List,
                    color = NULL,
                    n,
                    sizeNode = TRUE,
                    pal = rev(colorRampPalette(c('steelblue', '#f7fcfd', 'orange'))(5)),
                    range,
                    name,
                    fill) {

  if(is.null(pal)){
    pal = "grey"
  }

  if(!is.null(fill)){
  if (fill == "response") {
  if (sizeNode) {
    plot <- ggraph(List, "partition", weight = noObs) +
      geom_node_tile(aes(fill = var), size = 0.25) +
      geom_node_text(aes(label = ""), size = 4) +
      scale_y_reverse() +
      theme_void()
    if (!is.null(colors)) {
      plot <- plot + scale_fill_manual(values = color, name = "Variable") +
        ggnewscale::new_scale_fill() +
        ggnewscale::new_scale_color() +
        geom_node_tile(
          size = 0.15,
          data = . %>% filter(is.na(var)),
          aes(fill = respNode)
        ) +
        scale_fill_gradientn(
          colours = pal,
          limits = range,
          name = name,
          guide = guide_colorbar(
            frame.colour = "black",
            ticks.colour = "black",
            order = 2
          )
        )
    }
  } else {
    plot <- ggraph(List, "partition") +
      geom_node_tile(aes(fill = var), size = 0.25) +
      geom_node_text(aes(label = ""), size = 4) +
      scale_y_reverse() +
      theme_void()
    if (!is.null(colors)) {
      plot <- plot + scale_fill_manual(values = color, name = "Variable") +
        ggnewscale::new_scale_fill() +
        ggnewscale::new_scale_color() +
        geom_node_tile(
          size = 0.15,
          data = . %>% filter(is.na(var)),
          aes(fill = respNode)
        ) +
        scale_fill_gradientn(
          colours = pal,
          limits = range,
          name = name,
          guide = guide_colorbar(
            frame.colour = "black",
            ticks.colour = "black",
            order = 2
          )
        )
    }
  }
  }
  }

  if(!is.null(fill)){
  if (fill == "mu") {
    if (sizeNode) {
      plot <- ggraph(List, "partition", weight = noObs) +
        geom_node_tile(aes(fill = var), size = 0.25) +
        geom_node_text(aes(label = ""), size = 4) +
        scale_y_reverse() +
        theme_void()
      if (!is.null(colors)) {
        plot <- plot + scale_fill_manual(values = color, name = "Variable") +
          ggnewscale::new_scale_fill() +
          ggnewscale::new_scale_color() +
          geom_node_tile(
            size = 0.15,
            data = . %>% filter(is.na(var)),
            aes(fill = value)
          ) +
          scale_fill_gradientn(
            colours = pal,
            limits = range,
            name = name,
            guide = guide_colorbar(
              frame.colour = "black",
              ticks.colour = "black",
              order = 2
            )
          )
      }
    } else {
      plot <- ggraph(List, "partition") +
        geom_node_tile(aes(fill = var), size = 0.25) +
        geom_node_text(aes(label = ""), size = 4) +
        scale_y_reverse() +
        theme_void()
      if (!is.null(colors)) {
        plot <- plot + scale_fill_manual(values = color, name = "Variable") +
          ggnewscale::new_scale_fill() +
          ggnewscale::new_scale_color() +
          geom_node_tile(
            size = 0.15,
            data = . %>% filter(is.na(var)),
            aes(fill = value)
          ) +
          scale_fill_gradientn(
            colours = pal,
            limits = range,
            name = name,
            guide = guide_colorbar(
              frame.colour = "black",
              ticks.colour = "black",
              order = 2
            )
          )
      }
    }
  }
  }

  if(is.null(fill)){
    if (sizeNode) {
      plot <- ggraph(List, "partition", weight = noObs) +
        geom_node_tile(aes(fill = var), size = 0.25) +
        geom_node_text(aes(label = ""), size = 4) +
        scale_y_reverse() +
        theme_void()
      if (!is.null(colors)) {
        plot <- plot + scale_fill_manual(values = color, name = "Variable")
      }
    } else {
      plot <- ggraph(List, "partition") +
        geom_node_tile(aes(fill = var), size = 0.25) +
        geom_node_text(aes(label = ""), size = 4) +
        scale_y_reverse() +
        theme_void()
      if (!is.null(colors)) {
        plot <- plot + scale_fill_manual(values = color, name = "Variable")
    }
    }
  }

  return(plot)
}



# -------------------------------------------------------------------------


# Combine factor function for trees ---------------------------------------

cFactTrees <- function(treeData){
  dfOG <- treeData$data
  # find out which columns in my original data are factors
  factorColNam <- names(which(!(sapply(dfOG[colnames(dfOG)], is.numeric))))
  factorCols <- which((colnames(dfOG) %in% factorColNam))
  dfnew <- list()

  # create a list of the variables split into their factors
  if(any(class(treeData) == 'bart')){
    for (i in 1:length(factorCols)) {
      facLevels <- unique(dfOG[ ,factorCols[i]])
      dfnew[[i]] <- paste0(factorColNam[i],  as.numeric(facLevels))
    }
  }else if(any(class(treeData) ==  'bartMach')){

    for (i in 1:length(factorCols)) {
      facLevels <- unique(dfOG[,factorCols[i]])
      dfnew[[i]] <- paste0(factorColNam[i], "_", facLevels)
    }
  }else if(any(class(treeData) == 'dbarts')){

    for (i in 1:length(factorCols)) {
      facLevels <- unique(dfOG[,factorCols[i]])
      dfnew[[i]] <- paste0(factorColNam[i], ".", facLevels)
    }
  }


  # rename each element it's original name
  for (i in 1:length(dfnew)) {
    treeData$structure$var[which(treeData$structure$var %in% dfnew[[i]])] <- factorColNam[i]
  }

  return(treeData)
}


cFactTreesSelVars <- function(treeData){
  dfOG <- treeData$data
  # find out which columns in my original data are factors
  factorColNam <- names(which(!(sapply(dfOG[colnames(dfOG)], is.numeric))))
  factorCols <- which((colnames(dfOG) %in% factorColNam))


  # create a list of the variables split into their factors
  dfnew <- list()
  if(any(class(treeData) == 'bart')){
    for (i in 1:length(factorCols)) {
      facLevels <- unique(dfOG[ ,factorCols[i]])
      dfnew[[i]] <- paste0(factorColNam[i],  as.numeric(facLevels))
    }
  }else if(any(class(treeData) ==  'bartMach')){

    for (i in 1:length(factorCols)) {
      facLevels <- unique(dfOG[,factorCols[i]])
      dfnew[[i]] <- paste0(factorColNam[i], "_", facLevels)
    }
  }else if(any(class(treeData) == 'dbarts')){

    for (i in 1:length(factorCols)) {
      facLevels <- unique(dfOG[,factorCols[i]])
      dfnew[[i]] <- paste0(factorColNam[i], ".", facLevels)
    }
  }

  myList <- list(dfnew = dfnew, factorColNam = factorColNam)
  return(myList)
}




