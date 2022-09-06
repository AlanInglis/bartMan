library(lubridate)
library(BART)
library(bartMan)
library(dplyr)
library(ggplot2)


# read data
sbd <- read.csv('/Users/alaninglis/Desktop/SeoulBikeData.csv', sep = ',', check.names = F)


# data manipulation
sbd <- sbd |>
  mutate(dates = dmy(Date))

dd <- sbd |>
  group_by(dates, Season, Holiday) |>
  summarise(Count = sum(Count),
            Temp = mean(Temp),
            Humidity = mean(Humidity),
            Wind.Spd = mean(Wind.Spd),
            Visibility = mean(Visibility),
            Dew.Pt = mean(Dew.Pt),
            Solar.R = mean(Solar.R),
            Rainfall = mean(Rainfall),
            Snowfall = mean(Snowfall)
  ) |>
  mutate(Day = day(dates), Month = month(dates), Year = year(dates))

dd$Month <- month.abb[dd$Month]
dd$Day <- weekdays(as.Date(dd$dates,'%d-%m-%Y'))

# create weekend column
wkday <- c("Monday", "Tuesday" ,  "Wednesday" ,"Thursday", 'Friday' )
wke <- c( "Saturday",  "Sunday" )

dd <- transform(dd, Wkend = ifelse(dd$Day %in% wke, 'wkend', 'wkday'))

# make factor variables
dd$Season  <- as.factor(dd$Season)
dd$Holiday  <- as.factor(dd$Holiday)
dd$Month <- as.factor(dd$Month)
dd$Day <- as.factor(dd$Day)
dd$Year <- as.factor(dd$Year)
dd$Wkend <- as.factor(dd$Wkend)

# remove unnessesery columns
dd <- dd[,-c(1)]
dd <- dd[,-12]
# -------------------------------------------------------------------------


# remove zero counts
zeroC <- which(dd$Count == 0)
dd <- dd[-zeroC, ]

# transform skewed variables
logCols <- c('Count')
dd[logCols] <- (dd[logCols] + 1)^(1/3)




# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

dd <- as.data.frame(dd)
# BART model
xData <- dd[,-3]
yData <- dd[, 3]

#246
set.seed(8642)
bt <- wbart(x.train = xData,
            y.train = yData,
            nskip = 100,
            ndpost = 1000, # MCMC iters
            nkeeptreedraws = 1000,
            ntree = 100
)



# diagnostic plots
bartDiag(model = bt, response = dd$Count, burnIn = 100, combineFact = T, data = dd)

# tree df
btDF <- extractTreeData(model = bt, data = dd)
toneR::picR(1)
toneR::tone(2)



# -------------------------------------------------------------------------
# -------------------------------------------------------------------------


# plot standard heatmap
myMatStd <- viviBartMatrix(btDF,
                           type = 'standard',
                           metric = 'propMean', reorder = T,
                           combineFact = TRUE)

viviBartPlot(myMatStd, angle = 45)

myMat <- viviBartMatrix(btDF,
                        type = 'vsup',
                        metric = 'propMean',
                        metricError = "CV",
                        combineFact = T,
                        reorder = T)

# vsup
viviBartPlot(myMat, label = 'CV')  +
  ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 0, angle = 45))

# -------------------------------------------------------------------------
# tree plot
# -------------------------------------------------------------------------

# finding iteration with lowest residual sd
btPost <- bt$yhat.train

resid = NULL
for(i in 1:1000){
  resid[[i]] <- dd$Count - btPost[i,]
}

finalRes <- lapply(resid, sd)
which.min(finalRes)



# plot
plotAllTrees(btDF,
             iter = 996,
             selectedVars = c(15:26, 1:4, 7, 12, 13),
             cluster = 'var',
             removeStump = F,
             combineFact = T)




# -------------------------------------------------------------------------
# MDS
# -------------------------------------------------------------------------


# proximity matrix and MDS plots
bmProx <- proximityMatrix(btDF, dd, reorder = T, normalize = T, iter = 996)

mdsBart(treeData = btDF, data = dd, target = bmProx, plotType = 'interactive')





targetFit <- cmdscale(1 - bmProx, eig = TRUE, k = 2)

treeData = btDF
data= dd

save.image(file='sbdDailys.RData')


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
