#' mdsBart
#'
#' @description Multi-dimensional Scaling Plot of proximity matrix from a BART model.
#'
#' @param trees A data frame created by `extractTreeData` function.
#' @param plotType Type of plot. One of "interactive", "point", "rows", or "all".
#' @param data  a dataframe used in building the model.
#' @param response The name of the response for the fit.
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
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr bind_rows
#' @importFrom dplyr arrange
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr tibble
#' @importFrom ggiraph geom_polygon_interactive
#' @importFrom ggiraph girafe
#' @importFrom ggiraph opts_hover_inv
#' @importFrom ggiraph opts_hover
#' @importFrom grDevices rainbow
#'
#' @examples
#'if (requireNamespace("dbarts", quietly = TRUE)) {
#'  # Load the dbarts package to access the bart function
#'  library(dbarts)
#'  # Get Data
#'  df <- na.omit(airquality)
#'  # Create Simple dbarts Model For Regression:
#'  set.seed(1701)
#'  dbartModel <- bart(df[2:6],
#'    df[, 1],
#'    ntree = 5,
#'    keeptrees = TRUE,
#'    nskip = 10,
#'    ndpost = 10
#'  )
#'  # Tree Data
#'  trees_data <- extractTreeData(model = dbartModel, data = df)
#'  # Cretae Porximity Matrix
#'  bmProx <- proximityMatrix(
#'    trees = trees_data,
#'    reorder = TRUE,
#'    normalize = TRUE,
#'    iter = 1
#'  )
#'  # MDS plot
#'  mdsBart(
#'    trees = trees_data, data = df, target = bmProx,
#'    plotType = "interactive", level = 0.25, response = "Ozone"
#'  )
#'}
#'
#' @export

mdsBart <- function(
    trees,
    data,
    target,
    response,
    plotType = "rows",
    showGroup = TRUE,
    level = 0.95){

  if (!(plotType %in% c("interactive", "point", "rows", "all"))) {
    stop("plotType must be  \"interactive\", \"point\", \"rows\", or \"all\"")
  }

  # remove stumps
  whichStumps <- which(trees$structure$terminal == TRUE & trees$structure$node == 1)
  df <- NULL
  if(length(whichStumps > 0)){
    df$structure <- trees$structure[-whichStumps,]
  }else{
    df$structure <- trees$structure
  }


  # set target matrix
  targetFit <- stats::cmdscale(1 - target, eig = TRUE, k = 2)
  iter <- trees$nMCMC

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
      fitRot[[i]] <- stats::cmdscale(1 - rotationMatrix[[i]], eig = TRUE, k = 2)
    }
  )

  # perform procrustes on all
  message('Performing procrustes...')
  allProc <- list()
  suppressWarnings(
    for(i in 1:length(fitRot)){
      allProc[[i]] <- ProcrustesFn(targetFit$points, fitRot[[i]]$points)
    }
  )

  # turn into data frames
  dfRot <- list()
  for(i in 1:length(allProc)){
    dfRot[[i]] <- data.frame(x = allProc[[i]]$Yhat[,1],
                             y = allProc[[i]]$Yhat[,2])
  }

  # add in response & make a group variable for each df and combine
  # responseNum <- which(!(names(data) %in% trees$varName))
  # response <- names(data[responseNum])
  # responseNum <- which(names(data) == response)
  # response <- names(data[responseNum])

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
    summarise(kmX = bind_rows(stats::kmeans(x, 1)[2]),
              kmY = bind_rows(stats::kmeans(y, 1)[2]))
  dfM <- tibble(rowNo = dfM$rowNo,
                #pals = dfM$palette,
                mX = dfM$kmX$centers[,1],
                mY = dfM$kmY$centers[,1])


  # set up plot
  factorLevel <- as.factor(dfRotateAll$factResponse)
  nlevs <- nlevels(factorLevel)
  suppressWarnings(
    pal <-  grDevices::rainbow(nlevs)#RColorBrewer::brewer.pal(nlevs, "Set1")
  )
  pal <- pal[as.numeric(dfRotateAll$factResponse)]
  pal <- pal[as.numeric(dfM$rowNo)]

  # plot

  # set the limits
  rangeRot <- rbind(dfM$mX, dfM$mY)
  limitsRot <- range(rangeRot)
  limitsRot <- range(pretty(c(limitsRot[1], limitsRot[2])))
    #range(labeling::rpretty(limitsRot[1], limitsRot[2]))


  suppressMessages(
    suppressWarnings(
      if(plotType == 'interactive'){

        p <- ggplot(dfRotateAll, aes(x = x, y = y)) +
          #geom_text(label = dfRotateAll$rowNo) +
          geom_point(data = dfM, aes(x = mX, y = mY), col = pal) +
          #geom_text(data = dfM, aes(x = mX, y = mY, label = rowNo)) +
          stat_ellipse(alpha = 0, aes(col = rowNo), level = level)+
          xlab("Dimension 1") +
          ylab("Dimension 2") +
          theme_bw() +
          theme(legend.position = 'none')


        yLims <- layer_scales(p)$y$range$range
        xLims <- layer_scales(p)$x$range$range

        rangeRot <- rbind(yLims, xLims)
        limitsRotInt <- range(rangeRot)
        limitsRotInt <- range(pretty(c(limitsRotInt[1], limitsRotInt[2])))

        p <- p  +
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


        pFinal <- ggiraph::girafe(ggobj = pp,
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
          limitsRotInt <- range(pretty(c(limitsRotInt[1], limitsRotInt[2])))
            #range(labeling::rpretty(limitsRotInt[1], limitsRotInt[2]))


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
          limitsRotInt <- range(pretty(c(limitsRotInt[1], limitsRotInt[2])))
            #range(labeling::rpretty(limitsRotInt[1], limitsRotInt[2]))

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


ProcrustesFn <- function (X, Y) {
    n <- dim(X)[1]
    E <- diag(1, nrow = n)
    eins <- rep(1, n)
    k <- 1/n
    Z <- E - k * eins %*% t(eins)
    C <- t(X) %*% Z %*% Y
    s <- svd(C)
    f <- diag(s$d)
    P <- s$u
    Q <- s$v
    T <- Q %*% t(P)
    streck <- sum(diag((t(X) %*% Z %*% Y %*% T)))/sum(diag((t(Y) %*%
                                                              Z %*% Y)))
    trans <- as.vector(k * t(X - streck * Y %*% T) %*% eins)
    Yhut <- streck * Y %*% T + eins %*% t(trans)
    colnames(Yhut) <- rownames(T) <- colnames(T) <- names(trans) <- colnames(Y)
    dX <- stats::dist(X)
    dY <- stats::dist(Y)
    dYhat <- stats::dist(Yhut)
    cong <- sum(dX * dY)/(sqrt(sum(dX^2)) * sqrt(sum(dY^2)))
    alien <- sqrt(1 - cong^2)
    pairdist <- sort(sqrt(rowSums((X - Yhut)^2)))
    res <- list(X = X, Y = Y, Yhat = Yhut, translation = trans,
                dilation = streck, rotation = T, confdistX = dX, confdistY = dY,
                confdistYhat = dYhat, congcoef = cong, aliencoef = alien,
                pairdist = pairdist, call = match.call())
    class(res) <- "procr"
    return(res)
  }







