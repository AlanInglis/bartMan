% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VSUP.R
\name{pal_vsup}
\alias{pal_vsup}
\title{Variance suppressing uncertainty palette}
\usage{
pal_vsup(
  values,
  unc_levels = 4,
  max_light = 0.9,
  max_desat = 0,
  pow_light = 0.8,
  pow_desat = 1
)
}
\arguments{
\item{values}{Color values to be used at minimum uncertainty. Needs to be a vector of
length `2^unc_levels`.}

\item{unc_levels}{Number of discrete uncertainty levels. The number of discrete colors
at each level doubles.}

\item{max_light}{Maximum amount of lightening}

\item{max_desat}{Maximum amount of desaturation}

\item{pow_light}{Power exponent of lightening}

\item{pow_desat}{Power exponent of desaturation}
}
\description{
Returns a palette function that turns `v` (value) and `u` (uncertainty) (both between 0 and 1) into
colors.
}
