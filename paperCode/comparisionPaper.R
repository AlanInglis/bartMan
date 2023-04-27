
## Load libraries:
library(dbarts) # for model
library(bartMan) # for visualizations
library(vivid) # for agnostic visualizations
library(ggplot2) # for visualizations
library(scales) # for selecting colours
library(RColorBrewer) # for selecting colours


## Read in and setup data:
# Create some data
f <- function(x) {
  10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 +
    10 * x[, 4] + 5 * x[, 5]
}


set.seed(1701)
sigma <- 1.0
n <- 250
x <- matrix(runif(n * 10), n, 10)
colnames(x) <- paste0("x", 1:10)
Ey <- f(x)
y <- rnorm(n, Ey, sigma)
fData <- as.data.frame(cbind(x, y))

x <- fData[, 1:10]
y <- fData$y



## Build models
set.seed(1701)
dB20 <- bart(x.train = x,
             y.train = y,
             ntree = 20,
             keeptrees = TRUE,
             nskip = 100,
             ndpost = 1000
)

set.seed(1701)
dB100 <- bart(x.train = x,
              y.train = y,
              ntree = 100,
              keeptrees = TRUE,
              nskip = 100,
              ndpost = 1000
)

set.seed(1701)
dB200 <- bart(x.train = x,
              y.train = y,
              ntree = 200,
              keeptrees = TRUE,
              nskip = 100,
              ndpost = 1000
)

## Create dataframe of trees


dbT20 <- extractTreeData(model  = dB20, data = fData)
dbT100 <- extractTreeData(model = dB100, data = fData)
dbT200 <- extractTreeData(model = dB200, data = fData)






# -------------------------------------------------------------------------
# 20 trees ----------------------------------------------------------------
# -------------------------------------------------------------------------

# plot standard heatmap
myMatStd <- viviBartMatrix(dbT20,
                           type = 'standard',
                           metric = 'propMean',
                           reorder = F)

viviBartPlot(myMatStd)

# plot vsup
myMat <- viviBartMatrix(dbT20,
                        type = 'vsup',
                        metric = 'propMean',
                        metricError = "CV",
                        reorder = F)


viviBartPlot(myMat,label = 'CV')


# -------------------------------------------------------------------------


## vivid for 20 trees
response = 'y'
data = fData
# create predict function
responseIdx <- which(colnames(data) == response)
pFun <- function(fit, data, prob=TRUE) apply(predict(fit, data[,-responseIdx]), 2, mean)


# run vivid
set.seed(1701)
mat <- vivid::vivi(fit = dB20,
                   data = fData,
                   response = 'y',
                   reorder = F,
                   gridSize = 10,
                   nmax = 500,
                   normalized = FALSE,
                   class = 1,
                   predictFun = pFun)

colors <- scales::colour_ramp(
  colors = c(blue = '#FFFFCC', red = '#800026')
)((0:7)/7)
newCols <- RColorBrewer::brewer.pal(9, 'GnBu')
colors2 <- newCols[-1]


viviHeatmap(mat,
            intPal = colors,
            impPal = colors2)



#rm(data, fData,x, colors, colors2, Ey, n, newOrder, response, responseIdx,sigma, y, f, pFun)
#save.image(file='compare20.RData')


# -------------------------------------------------------------------------
# 100 trees ---------------------------------------------------------------
# -------------------------------------------------------------------------


# plot standard heatmap
myMatStd <- viviBartMatrix(dbT100,
                           type = 'standard',
                           metric = 'propMean',
                           reorder = F)

viviBartPlot(myMatStd, impLims = c(0, 0.16))

# plot vsup
myMat <- viviBartMatrix(dbT100,
                        type = 'vsup',
                        metric = 'propMean',
                        metricError = "CV",
                        reorder = F)


viviBartPlot(myMat, label = 'CV', impLims = c(0, 0.16))





# -------------------------------------------------------------------------
## vivid for 100 trees


response = 'y'
data = fData
# create predict function
responseIdx <- which(colnames(data) == response)
pFun <- function(fit, data, prob=TRUE) apply(predict(fit, data[,-responseIdx]), 2, mean)


# run vivid
set.seed(1701)
mat <- vivid::vivi(fit = dB100,
                   data = fData,
                   response = 'y',
                   reorder = F,
                   gridSize = 10,
                   nmax = 500,
                   normalized = FALSE,
                   class = 1,
                   predictFun = pFun)

colors <- scales::colour_ramp(
  colors = c(blue = '#FFFFCC', red = '#800026')
)((0:7)/7)
newCols <- RColorBrewer::brewer.pal(9, 'GnBu')
colors2 <- newCols[-1]

viviHeatmap(mat,
            intPal = colors,
            impPal = colors2)




#rm(data, fData,x, colors, colors2, Ey, n, newOrder, response, responseIdx,sigma, y, f, pFun)
#save.image(file='compare100.RData')

# -------------------------------------------------------------------------
# 200 trees ---------------------------------------------------------------
# -------------------------------------------------------------------------


# plot standard heatmap
myMatStd <- viviBartMatrix(dbT200,
                           type = 'standard',
                           metric = 'propMean',
                           reorder = F)

viviBartPlot(myMatStd, impLims = c(0, 0.14))

# plot vsup
myMat <- viviBartMatrix(dbT200,
                        type = 'vsup',
                        metric = 'propMean',
                        metricError = "CV",
                        reorder = F)


viviBartPlot(myMat, label = 'CV',impLims = c(0, 0.14))

# -------------------------------------------------------------------------

## vivid for 200 trees
response = 'y'
data = fData
# create predict function
responseIdx <- which(colnames(data) == response)
pFun <- function(fit, data, prob=TRUE) apply(predict(fit, data[,-responseIdx]), 2, mean)


# run vivid
set.seed(1701)
mat <- vivid::vivi(fit = dB200,
                   data = fData,
                   response = 'y',
                   reorder = F,
                   gridSize = 10,
                   nmax = 500,
                   normalized = FALSE,
                   class = 1,
                   predictFun = pFun)

colors <- scales::colour_ramp(
  colors = c(blue = '#FFFFCC', red = '#800026')
)((0:7)/7)
newCols <- RColorBrewer::brewer.pal(9, 'GnBu')
colors2 <- newCols[-1]

viviHeatmap(mat,
            intPal = colors,
            impPal = colors2)


#rm(data, fData,x, colors, colors2, Ey, n, newOrder, response, responseIdx,sigma, y, f, pFun)
#save.image(file='compare200.RData')






