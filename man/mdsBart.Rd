% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mdsBart.R
\name{mdsBart}
\alias{mdsBart}
\title{mdsBart}
\usage{
mdsBart(
  treeData,
  data,
  target,
  response,
  plotType = "rows",
  showGroup = TRUE,
  level = 0.95
)
}
\arguments{
\item{treeData}{A data frame created by treeData function.}

\item{data}{a dataframe used in building the model.}

\item{target}{A target proximity matrix to}

\item{response}{The name of the response for the fit.}

\item{plotType}{Type of plot to show. Either 'interactive' - showing interactive confidence ellipses.
'point' - a point plot showing the average position of a observation.
'rows' - displaying the average position of a observation number instead of points.
'all' - show all observations (not averaged).}

\item{showGroup}{Logical. Show confidence ellipses.}

\item{level}{The confidence level to show. Default is 95\% confidence level.}
}
\value{
For this function, the MDS coordinates are calculated for each iteration.
Procrustes method is then applied to align each of the coordinates to a target set
of coordinates. The returning result is then a clustered average of each point.
}
\description{
Multi-dimensional Scaling Plot of proximity matrix from a BART model.
}
