library(reshape2)
library(ggplot2)
library(multiscales)

library(scales)
library(dplyr)
library(purrr)

set.seed(100)
# create some data
x <- LETTERS[1:5]
s <- matrix(runif(25, 0, 5),5,5)
s[lower.tri(s)] <- t(s)[lower.tri(s)]

rownames(s) <- colnames(s) <- x
s <- melt(s)

# create fake uncertainty vals
s$uncertainty <- runif(25)

s$Var1 <- factor(s$Var1, levels = x)
pal =  rev(colorspace::sequential_hcl(palette = "Purples 3", n = 8))

#lim <- data.frame(a = c(0, 5), b = c(0, 1))

# works
pp <- ggplot(s) +
  geom_raster(aes(x = Var1, y = Var2, fill = zip(value, uncertainty))) +
  bivariate_scale(
    name = c("Value", "Uncertainty"),
    aesthetics = "fill",
    limits = list(c(0, 5), c(0, 1)),
    palette = pal_vsup(
      values = rev(colorspace::sequential_hcl(palette = "Purples 3", n = 8)),
      ),
    guide = "colorfan"
  )

pp


lims <- tibble(
  x = c(0.5,5),
  y = c(0,1)
)


# fails
p <- ggplot(s) +
  geom_raster(aes(x = Var1, y = Var2, fill = zip(value, uncertainty))) +
  bivariate_scale(
    name = c("Value", "Uncertainty"),
    aesthetics = "fill",
    limits = list(c(1, 5), c(0, 1)), # changed 0 to 0.5 for the limits of value
    breaks = list(1:5, c(0,0.5,1)),
    palette = pal_vsup(
      values = rev(colorspace::sequential_hcl(palette = "Purples 3", n = 8)),
    ),
    #guide = "colorbox" # changing to colorbox fixes error
    guide = "colorfan"
  )
p

# -------------------------------------------------------------------------

px <- ggplot(s) +
  geom_raster(aes(x = Var1, y = Var2, fill = zip(value, uncertainty)))

