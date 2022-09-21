
library(bartMachine)
library(dbarts)
library(bartMan)
library(ggplot2)

# load data
data(iris)
iris2 = iris[1:100,]
iris2$Species <- factor(iris2$Species)
iris2$Species <- ifelse(iris2$Species == "setosa", 0, 1)

# bartMachine
set.seed(100)
bm <-  bartMachine(X = iris2[,1:4],
                   y = iris2[,5],
                   num_trees = 20,
                   seed = 100)



# dbarts
set.seed(100)
dB <- bart(x.train = iris2[,1:4],
           y.train = iris2[,5],
           ntree = 20,
           keeptrees = TRUE,
           nskip = 250,
           ndpost = 1000
)


# -------------------------------------------------------------------------
# Create data frame of trees ---------------------------------------------
bmDF <- extractTreeData(bm, iris2)
dbDF <- extractTreeData(dB, iris2)

# -------------------------------------------------------------------------
# VIVI bart ---------------------------------------------------------------
vimpPlot(bmDF, plotType = 'point', metric = 'median')



# plot standard heatmap
myMatStd <- viviBartMatrix(bmDF,
                           type = 'standard',
                           metric = 'propMean', reorder = T)

viviBartPlot(myMatStd, impLims = c(0, 0.45))

# plot vsup
myMat <- viviBartMatrix(bmDF,
                        type = 'vsup',
                        metric = 'propMean',
                        metricError = "CV", reorder = T)


viviBartPlot(myMat,
             impLims = c(0,0.45),
             label = 'CV')

# -------------------------------------------------------------------------
# tree plots --------------------------------------------------------------

# finding iteration with lowest residual sd
bmPost <- bartMachine::bart_machine_get_posterior(bm, iris2[,1:4])

resid = NULL
for(i in 1:1000){
  resid[[i]] <- iris2[,5] - bmPost$y_hat_posterior_samples[,i]
}

finalRes <- lapply(resid, sd)
sort(unlist(finalRes))
which.min(finalRes)

# plot all trees
plotAllTrees(bmDF, treeNo = 20)

plotAllTrees(bmDF, iter = 736, sizeNode = T, fillBy = 'mu')
plotAllTrees(bmDF, iter = 736, sizeNode = T, fillBy = 'response')
plotAllTrees(bmDF, iter = 736, cluster = "depth")
plotAllTrees(bmDF, iter = 736, cluster = "var")


treeBarPlot(bmDF, topTrees = 10)




# -------------------------------------------------------------------------
# Diagnostic plots --------------------------------------------------------

# acceptance rate
acceptRate(bmDF) + ylim(c(0, 1))
acceptRate(dbDF) + ylim(c(0, 1))


# tree depth per iteration
treeDepth(bmDF) + ylim(c(0.7, 1.6))
treeDepth(dbDF) + ylim(c(0.7, 1.6))

# tree nodes per iteration
treeNodes(bmDF) + ylim(c(2.3, 4.3))
treeNodes(dbDF) + ylim(c(2.3, 4.3))

# split density
splitDensity(bmDF, data = iris2, display = 'dataSplit')
splitDensity(dbDF, data = iris2, display = 'both2')




# -------------------------------------------------------------------------

# mds plot
bmProx <- proximityMatrix(bmDF, iris2, reorder = T, normalize = T, iter = 736)
mdsBart(treeData = bmDF, data = iris2, target =  bmProx,
        plotType = 'interactive', showGroup = F)


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------


#save.image(file = 'iris.Rdata')







