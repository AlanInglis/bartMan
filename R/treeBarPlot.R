#' treeBarPlot
#'
#' @description Creates a barplot displaying the frequency of different tree structures.
#'
#' @param treeList A list of trees created using the treeList function.
#' @param topTrees integer value to show the top x variables.
#' @param iter The selected iteration
#' @param treeNo The selected tree number.
#' @param removeStump LOGICAL. If TRUE, then stumps are removed from plot. If False, stumps
#' remain in plot and are coloured grey.
#'
#' @return A barplot plot.
#'
#'
#' @import ggplot2
#' @importFrom purrr map_chr
#' @importFrom purrr map2
#' @importFrom tidygraph activate
#' @importFrom tidygraph tbl_graph
#' @importFrom tidyr separate
#' @importFrom tidyr replace_na
#' @importFrom dplyr as_tibble
#' @importFrom dplyr arrange
#' @importFrom dplyr mutate
#' @importFrom dplyr slice
#' @importFrom igraph gsize
#' @importFrom ggraph ggraph
#' @importFrom ggraph geom_node_tile
#' @importFrom ggraph geom_node_text
#' @importFrom patchwork area
#' @importFrom patchwork plot_layout
#' @importFrom cowplot get_legend
#' @importFrom cowplot plot_grid
#'
#'
#' @export

treeBarPlot <- function(treeData,
                        iter = NULL,
                        treeNo = NULL,
                        topTrees = NULL,
                        removeStump = FALSE){


  treeList <- plotAll(treeData, iter = iter, treeNo = treeNo, cluster = NULL)

  # remove stumps
  #treeList <- Filter(function(x) igraph::gsize(x) > 0, treeList)

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
        newTreesDF[[i]]$var <- "Stump"
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



  # get the frequency of similar trees:
  freqs <- map(treeList, function(x){
    x %>%
      pull(var) %>%
      tidyr::replace_na("..") %>%
      paste0(collapse = "")
  }) %>%
    unlist(use.names = F) %>%
    as_tibble() %>%
    mutate(ids = 1:n()) %>%
    group_by(value) %>%
    mutate(val = n():1)


  freqDf <-  freqs %>%
    slice(1) %>%
    arrange(-val) %>%
    rename(frequency = val)  # frequency tibble
  freqDf$treeNum <- seq(1:nrow(freqDf)) # add tree number


  if(!is.null(topTrees)){
    freqDf <- freqDf[1:topTrees,]
  }

  lengthFreq <- length(freqDf$value)
  ids <- freqDf$ids
  #ids <- freqs %>% slice(1) %>% pull(ids) # remove duplicates
  freqs <- freqs[ids,] %>% pull(val) # get frequencies


  treeList <- treeList[ids]
  treeListNew <- purrr::imap(treeList, ~.x %>%
                               mutate(frequency = freqs[.y]) %>%
                               select(var, frequency))

  # # return new list of trees
  # treeList <- treeListNew[sort(ids)]
    treeList <- treeListNew

  # add plot name as number
  # for(i in 1:(length(treeList))){
  #   treeList[[i]] <- treeList[[i]] %>%
  #     activate(nodes) %>%
  #     mutate(name = c(i, rep("", length.out = igraph::gsize(treeList[[i]]))))
  # }


  # Create barplot of frequency ---------------------------------------------
  names <- factor(freqDf$value, levels = freqDf$value)

  bp <- freqDf %>%
    ggplot() +
    geom_bar(aes(x = value, y = frequency), fill = 'steelblue', stat = "identity") +
    scale_x_discrete(limits = rev(levels(names))) +
    ggtitle("") +
    ylab("Count") +
    xlab("") +
    theme_bw() +
    theme(legend.position = "none") +
    coord_flip()


  # ggraph plotting funtion -------------------------------------------------


  # set node colours
  nodenames <- unique(na.omit(unlist(lapply(treeList, .%>%activate(nodes) %>% pull(var)))))
  nodenames <- sort(nodenames)
  nodecolors <- setNames(scales::hue_pal(c(0,360)+15, 100, 64, 0, 1)(length(nodenames)), nodenames)

  if(length(stumpIdx) >=  1){
      nodecolors[["Stump"]] <- '#808080'
  }



  # plotting function
  plotFun <- function(List, colors = NULL, n) {

    plot <- ggraph(List, "partition") +
      geom_node_tile(aes(fill = var), size = 0.25) +
      geom_node_text(aes(label = ''), size = 4) +
      theme(legend.position = "bottom") +
      scale_y_reverse() +
      theme_void()
    if (!is.null(colors)) {
      plot <- plot + scale_fill_manual(values = colors, name = "Variable") +
        scale_color_manual(values = colors, na.value = "grey")  +
        theme(legend.position = "bottom")
    }
  }


  allPlots <- lapply(treeList, plotFun, n = length(treeList), color = nodecolors)

  # get legend
  legend <- cowplot::get_legend(allPlots[[1]])

  # remove legends from individual plots
  allPlots <- lapply(allPlots, function(x) x + theme(legend.position = "none"))

   if(removeStump == FALSE){

    whichStumpData = NULL
    for(i in 1:length(allPlots)){
      whichStumpData[[i]] <-  which(allPlots[[i]]$data$leaf[1] == TRUE)
    }

    stumpIdxNew <- which(whichStumpData == 1)
    for(i in stumpIdxNew){
      allPlots[[i]]$data <- allPlots[[i]]$data[-2, ]
    }
   }

  # filter top X% of plots
  if(!is.null(topTrees)){
    allPlots <- allPlots[1:topTrees]
  }

  # Create final barplot ----------------------------------------------------
  width = 1
  p_axis <- ggplot(freqDf) +
    geom_blank(aes(y = value)) +
    purrr::map2(allPlots,
                rev(seq_along(allPlots)),
                ~ annotation_custom(ggplotGrob(.x),
                                    ymin = .y - width / 2,
                                    ymax = .y + width / 2,
                                    )) +
    theme_void()

  bp1 <- bp + theme(axis.text.y = element_blank())
  ppp <- p_axis + theme(aspect.ratio = 10)
  px <- ppp|bp1

  bpFinal <- cowplot::plot_grid(px, legend, rel_heights = c(1, .1), ncol = 1)

  return(bpFinal)

}



