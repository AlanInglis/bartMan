% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treeNodes.R
\name{treeNodes}
\alias{treeNodes}
\title{treeNodes}
\usage{
treeNodes(trees)
}
\arguments{
\item{trees}{A list of tree attributes created using the extractTreeData function.}
}
\value{
A plot of tree number of nodes over iterations.
}
\description{
A plot of number of nodes over iterations.
}
\examples{
if(requireNamespace("dbarts", quietly = TRUE)){
 # Load the dbarts package to access the bart function
 library(dbarts)
 # Get Data
 df <- na.omit(airquality)
 # Create Simple dbarts Model For Regression:
 set.seed(1701)
 dbartModel <- bart(df[2:6], df[, 1], ntree = 5, keeptrees = TRUE, nskip = 10, ndpost = 10)

 # Tree Data
 trees_data <- extractTreeData(model = dbartModel, data = df)
 treeNodes(trees = trees_data)
}

}
