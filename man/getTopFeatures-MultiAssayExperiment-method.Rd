% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods-getTopFeatures.R
\name{getTopFeatures,MultiAssayExperiment-method}
\alias{getTopFeatures,MultiAssayExperiment-method}
\title{Get the Top Ranked Features in a \code{MultiAssayExperiment} object}
\usage{
\S4method{getTopFeatures}{MultiAssayExperiment}(object, numFeats, ...)
}
\arguments{
\item{object}{A \code{SummarizedExperiment} to extract top features from}

\item{numFeats}{An \code{integer} number of top ranked features to extract.}

\item{...}{Fall through arguments to \code{rankFeatures}.}
}
\value{
A \code{character} vector of top ranked features, with the features in
order of rank ascending.
}
\description{
Get the Top Ranked Features in a \code{MultiAssayExperiment} object
}
\seealso{
\code{\link{rankFeatures,MultiAssayExperiment-method}}
}
