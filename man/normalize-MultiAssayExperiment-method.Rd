% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods-normalize.R
\name{normalize,MultiAssayExperiment-method}
\alias{normalize,MultiAssayExperiment-method}
\title{Normalize the \code{assays} of a \code{MutliAssayExperiment} Object}
\usage{
\S4method{normalize}{MultiAssayExperiment}(
  object,
  MARGIN = 2,
  FUN = "preprocessCaret",
  ...,
  whichAssays = seq_along(assays(object))
)
}
\arguments{
\item{object}{A \code{SummarizedExperiment} object with assays to normalize.}

\item{MARGIN}{An \code{integer} indicating if rows (1) or columns (2) should be
normalized. Defaults to 2 for columns.}

\item{FUN}{A function to normalize your data with. Should accept a
rectangular object such as a \code{matrix}, \code{data.frame}, or \code{data.table} and
return an object of the same class with the data normalized using FUN.}

\item{...}{Fall through parameters to FUN. For the default FUN, these are
passed to \code{caret::preProcess} to allow configuration of the normalization
method. Omitting any arguments with the default FUN will scale and center
the data.}

\item{whichAssays}{A \code{numeric} or \code{character} vector specifying the indices
of the assays to normalize. Defaults to all assays.}
}
\value{
The \code{MultiAssayExperiment} with one or more of the assays normalized
and information about the normalization method in the \code{normalization} item
of the object \code{metadata}.
}
\description{
For this method to work, there must be a \code{normalize} method defined
for all classes of experiments in the \code{MultiAssayExperiment}
}
\details{
When using the default FUN, it is also possible to impute missing values.
See \code{?caret::preProcess} for information on available methods.
}
\seealso{
\code{\link{preprocessCaret}}, \code{\link[caret:preProcess]{caret::preProcess}}
}
