#' mdsBart
#'
#' @description Multi-dimensional Scaling Plot of proximity matrix from a BART model.
#'
#' @param treeData A data frame created by treeData function.
#' @param plotType xxx
#' @param data  a dataframe
#' @param target A target proximity matrix to
#' @param plotType Type of plot to show. Either 'interactive' - showing interactive confidence ellipses.
#' 'point' - a point plot showing the average position of a observation.
#' 'rows' - displaying the average position of a observation number instead of points.
#' 'all' - show all observations (not averaged).
#' @param showGroup Logical. Show confidence ellipses.
#' @param level The confidence level to show. Default is 95\% confidence level.
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
#' @importFrom ggiraph geom_polygon_interactive
#' @importFrom ggiraph ggiraph
#' @importFrom ggiraph opts_hover_inv
#' @importFrom ggiraph opts_hover
#'
#' @export
#'

mdsBart <- function(
    treeData,
    data,
    target,
    plotType = "rows",
    showGroup = TRUE,
    level = 0.95){

  if (!(plotType %in% c("interactive", "point", "rows", "all"))) {
    stop("plotType must be  \"interactive\", \"point\", \"rows\", or \"all\"")
  }

  # remove stumps
  whichStumps <- which(treeData$structure$isLeaf == TRUE & treeData$structure$node == 1)
  df <- NULL
  if(length(whichStumps > 0)){
    df$structure <- treeData$structure[-whichStumps,]
  }else{
    df$structure <- treeData$structure
  }


  # set target matrix
  targetFit <- cmdscale(1 - target, eig = TRUE, k = 2)
  iter <- treeData$nMCMC

  # get all rotation matrices
  message('Getting proximites...')
  rotationMatrix <- list()
  for(i in 1:iter){
    rotationMatrix[[i]] <- bartMan::proximityMatrix(df,
                                                    data,
                                                    reorder = F,
                                                    normalize = T,
                                                    iter = i)
  }

  # get all MDS fits
  message('Getting MDS...')
  fitRot <- list()
  suppressWarnings(
    for(i in 1:length(rotationMatrix)){
      fitRot[[i]] <- cmdscale(1 - rotationMatrix[[i]], eig = TRUE, k = 2)
    }
  )

  # perform procrustes on all
  message('Performing procrustes...')
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

  # turn response into something usable
  dfRotateAll$response <- as.numeric(as.factor(dfRotateAll$response))
  #dfRotateAll$palette <- pal1[dfRotateAll$response]
  if(max(dfRotateAll$response) < 2){
    dfRotateAll <- dfRotateAll %>%
      mutate(factResponse = ifelse(response ==  0 ,1,2))
  }else{
    dfRotateAll$factResponse <- dfRotateAll$response
  }

  # create data frame to plot
  dfM <-  dfRotateAll %>%
    arrange(rowNo) %>%
    group_by(rowNo) %>%
    summarise(kmX = bind_rows(kmeans(x, 1)[2]),
              kmY = bind_rows(kmeans(y, 1)[2]))
  dfM <- tibble(rowNo = dfM$rowNo,
                #pals = dfM$palette,
                mX = dfM$kmX$centers[,1],
                mY = dfM$kmY$centers[,1])


  # set up plot
  factorLevel <- as.factor(dfRotateAll$factResponse)
  nlevs <- nlevels(factorLevel)
  suppressWarnings(
    pal <-  RColorBrewer::brewer.pal(nlevs, "Set1")
  )
  pal <- pal[as.numeric(dfRotateAll$factResponse)]
  pal <- pal[as.numeric(dfM$rowNo)]

  # plot

  # set the limits
  rangeRot <- rbind(dfM$mX, dfM$mY)
  limitsRot <- range(rangeRot)
  limitsRot <- range(labeling::rpretty(limitsRot[1], limitsRot[2]))


  suppressMessages(
    suppressWarnings(
      if(plotType == 'interactive'){
        # set the limits
        rangeRot <- rbind(dfRotateAll$x, dfRotateAll$y)
        limitsRotInt <- range(rangeRot)
        limitsRotInt <- range(labeling::rpretty(limitsRotInt[1], limitsRotInt[2]))

        p <- ggplot(dfRotateAll, aes(x = x, y = y)) +
          #geom_text(label = dfRotateAll$rowNo) +
          geom_point(data = dfM, aes(x = mX, y = mY), col = 'steelblue') +
          #geom_text(data = dfM, aes(x = mX, y = mY, label = rowNo)) +
          stat_ellipse(alpha = 0, aes(col = rowNo), level = level)+
          xlab("Dimension 1") +
          ylab("Dimension 2") +
          theme_bw() +
          theme(legend.position = 'none') +
          scale_x_continuous(limits = limitsRotInt) +
          scale_y_continuous(limits = limitsRotInt)


        build <- ggplot_build(p)$data
        ell <- build[[2]]
        ell$group <- as.factor(ell$group)
        levels(ell$group) <- levels(as.factor(dfM$rowNo))

        pp <- p + ggiraph::geom_polygon_interactive(data = ell,
                                           aes(x, y,
                                               group = group,
                                               fill  = group,
                                               tooltip = group,
                                               data_id = group),
                                           fill = ell$group,
                                           alpha = 0.3)


        pFinal <- ggiraph::ggiraph(ggobj = pp,
                          options = list(
                            ggiraph::opts_hover_inv(css = "opacity:0.05;"),
                            ggiraph::opts_hover(ggiraph::girafe_css(
                              css = "fill:blue;stroke:gray;",
                              line = "fill:blue",
                              area = "stroke-width:10px",
                              point = "stroke-width:10px",
                              image = "outline:2px red"
                            ))
                          )
        )
      }else{
        if(showGroup){
          rangeRot <- rbind(dfRotateAll$x, dfRotateAll$y)
          limitsRotInt <- range(rangeRot)
          limitsRotInt <- range(labeling::rpretty(limitsRotInt[1], limitsRotInt[2]))


          p <- ggplot(dfRotateAll, aes(x = x, y = y)) +
            stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill = rowNo), level = level)+
            xlab("Dimension 1") +
            ylab("Dimension 2") +
            theme_bw() +
            theme(legend.position = 'none') +
            scale_x_continuous(limits = limitsRotInt) +
            scale_y_continuous(limits = limitsRotInt)

          p
          # p <- ggplot(data = dfRotateAll, aes(x = x, y = y)) +
          #   ggforce::geom_mark_ellipse(data = dfRotateAll,
          #                              aes(x = x,
          #                                  y = y,
          #                                  fill = factor(rowNo),
          #                                  col = rowNo),
          #                              alpha = 0.1)
          if(plotType == 'point'){
            p <- p +
              geom_point(data = dfM, aes(x = mX, y = mY), col = pal)
          }else if(plotType == 'rows'){
            p <- p +
              geom_text(data = dfM, aes(x = mX, y = mY, label = rowNo), col = pal)
          }else if(plotType == "all"){
            p <- p +
              geom_text(label = dfRotateAll$rowNo, col = dfRotateAll$rowNo)
          }
        }else{
          rangeRot1 <- rbind(dfRotateAll$x, dfRotateAll$y)
          limitsRotInt <- range(rangeRot1)
          limitsRotInt <- range(labeling::rpretty(limitsRotInt[1], limitsRotInt[2]))

          p <- ggplot(dfRotateAll, aes(x = x, y = y)) +
            xlab("Dimension 1") +
            ylab("Dimension 2") +
            theme_bw() +
            theme(legend.position = 'none') +
            scale_x_continuous(limits = limitsRot) +
            scale_y_continuous(limits = limitsRot)

          if(plotType == 'point'){
            p <- p +
              geom_point(data = dfM, aes(x = mX, y = mY), col = pal)
          }else if(plotType == 'rows'){
            p <-   p +
              geom_text(data = dfM, aes(x = mX, y = mY, label = rowNo), col = pal)
          }else if(plotType == "all"){
            p <- p +
              geom_text(label = dfRotateAll$rowNo, col = dfRotateAll$rowNo)  +
              scale_x_continuous(limits = limitsRotInt) +
              scale_y_continuous(limits = limitsRotInt)
          }
        }


        pFinal <- p +
          theme_bw() +
          theme(legend.position = 'none')
      }
    )
  )

  #return(pFinal)
  suppressWarnings(print(pFinal))
}

