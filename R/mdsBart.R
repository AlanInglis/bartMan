#' mdsBart
#'
#' @description Multi-dimensional Scaling Plot of proximity matrix from a BART model.
#'
#' @param treeData A data frame created by treeData function.
#' @param type What value to return. Either the raw count 'val', the proportion 'prop'
#' or the column means of the proportions 'propMean'
#'
#' @return For this function, the MDS coordinates are calculated for each iteration.
#' Procrustes method is then applied to align each of the coordinates to a target set
#' of coordinates. The returning result is then a clustered average of each point.
#'
#' @import ggplot2
#' @importFrom smacof Procrustes
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr bind_rows
#' @importFrom dplyr arrange
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr tibble
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggforce geom_mark_ellipse
#'
#' @export
#'

mdsBart <- function(treeData,
                    data,
                    target,
                    type = "mean",
                    plotType = "rows",
                    showGroup = TRUE){

  if (!(type %in% c("mean", "kmeans"))) {
    stop("type must be \"mean\" or \"kmeans\"")
  }

  if (!(plotType %in% c("point", "rows", "all", 'none'))) {
    stop("plotType must be \"point\", \"rows\", \"all\", or \"none\"")
  }

  # remove stumps
  whichStumps <- which(treeData$structure$isLeaf == TRUE & treeData$structure$node == 1)
  df <- NULL
  df$structure <- treeData$structure[-whichStumps,]

  # set target matrix
  targetFit <- cmdscale(1 - target, eig = TRUE, k = 2)
  iter <- treeData$nMCMC

  # get all rotation matrices
  rotationMatrix <- list()
  for(i in 1:iter){
    rotationMatrix[[i]] <- bartMan::proximityMatrix(df,
                                           data,
                                           reorder = F,
                                           normalize = T,
                                           iter = i)
  }

  # get all MDS fits
  fitRot <- list()
  suppressWarnings(
    for(i in 1:length(rotationMatrix)){
      fitRot[[i]] <- cmdscale(1 - rotationMatrix[[i]], eig = TRUE, k = 2)
    }
  )

  # perform procrustes on all
  allProc <- list()
  suppressWarnings(
    for(i in 1:length(fitRot)){
      allProc[[i]] <- smacof::Procrustes(targetFit$points, fitRot[[i]]$points)
    }
  )

  # turn into data frames
  dfRot <- list()
  for(i in 1:length(allProc)){
    dfRot[[i]] <- data.frame(x = allProc[[i]]$Yhat[,1],
                             y = allProc[[i]]$Yhat[,2])
  }

  # add in response & make a group variable for each df and combine
  responseNum <- which(!(names(data) %in% treeData$varName))
  response <- names(data[responseNum])

  addResponse <- function(x){
     d <- x %>%
      mutate(response = data[response][,1])
     return(d)
  }

  dfRot <- lapply(dfRot, addResponse)

  # make a group variable for each df and combine
  addRowID <- function(x) {
    d <- x
    return(d %>% mutate(rowNo = rownames(d)))
  }

  dfRot <- lapply(dfRot, addRowID)

  # put together
  dfRotateAll <- bind_rows(lapply(dfRot,
                                  function(x){
                                    addRowID(x)
                                  }),
                           .id = "dfID")

  # turn response into something useable
  dfRotateAll <- dfRotateAll %>%
    mutate(factResponse = ifelse(response ==  0 ,1,2))

  if(type == 'mean'){

    dfRotateAll$rowNo <- as.numeric(dfRotateAll$rowNo)
    dfM <- dfRotateAll %>%
      arrange(rowNo) %>%
      group_by(rowNo) %>%
      dplyr::summarize(mX = mean(x), mY = mean(y))
    dfRotateAll$rowNo <-  as.character(dfRotateAll$rowNo)

  }else if(type == "kmeans"){

   dfM <-  dfRotateAll %>%
      arrange(rowNo) %>%
      group_by(rowNo) %>%
      summarise(kmX = bind_rows(kmeans(x, 1)[2]),
                kmY = bind_rows(kmeans(y, 1)[2]))
   dfM <- tibble(rowNo = dfM$rowNo,
                 mX = dfM$kmX$centers[,1],
                 mY = dfM$kmY$centers[,1])
  }

  # set up plot
  factorLevel <- as.factor(dfRotateAll$factResponse)
  nlevs <- nlevels(factorLevel)
  suppressWarnings(
    pal <-  RColorBrewer::brewer.pal(nlevs, "Set1")
    )
  pal <- pal[as.numeric(dfRotateAll$factResponse)]
  pal <- pal[as.numeric(dfM$rowNo)]

  # plot

  if(showGroup){
    p <- ggplot(data = dfRotateAll, aes(x = x, y = y)) +
      ggforce::geom_mark_ellipse(data = dfRotateAll,
                                 aes(x = x,
                                     y = y,
                                     fill = factor(rowNo),
                                     col = rowNo),
                                 alpha = 0.1)
    if(plotType == 'point'){
      p <- p + #ggplot(dfM, aes(x = mX, y = mY)) +
        geom_point(data = dfM, aes(x = mX, y = mY), col = pal)
    }else if(plotType == 'rows'){
      p <- p + #ggplot(dfM, aes(x = mX, y = mY)) +
        geom_text(data = dfM, aes(x = mX, y = mY, label = rowNo), col = pal)
    }else if(plotType == "all"){
      p <- p +
        geom_text(label = dfRotateAll$rowNo, col = dfRotateAll$rowNo)
    }else if(plotType == 'none'){
      p <- p
    }
  }else{
    if(plotType == 'point'){
      p <-  ggplot(dfM, aes(x = mX, y = mY)) +
        geom_point(col = pal)
    }else if(plotType == 'rows'){
      p <-  ggplot(dfM, aes(x = mX, y = mY)) +
        geom_text(label = dfM$rowNo, col = pal)
    }else if(plotType == "all"){
      p <- ggplot(dfRotateAll, aes(x = x, y = y)) +
        geom_text(label = dfRotateAll$rowNo, col = dfRotateAll$rowNo)
    }else if(plotType == 'none'){
      stop("plotType cannot be 'none' if showGroup is FALSE")
    }
  }


  p <- p +
    xlab("Dimension 1") +
    ylab("Dimension 2") +
    theme_bw() +
    theme(legend.position = 'none')

  return(p)
}









