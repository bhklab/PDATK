% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{hasColDataColumns}
\alias{hasColDataColumns}
\title{Check for Column Names in the colData Slot of a SummarizedExperiment}
\usage{
hasColDataColumns(SE, values)
}
\arguments{
\item{SE}{A \code{SummarizedExperiment} object to check for the existence
of colData columns.}

\item{values}{A \code{character} vector with one or more column name to
check for in the column data.}
}
\value{
\code{logical} True if all of values are column names in the
\code{SummarizedExperiment} object, FALSE otherwise.
}
\description{
Check for Column Names in the colData Slot of a SummarizedExperiment
}
\examples{
SE <- SummarizedExperiment(matrix(rnorm(100), 10, 10),
    colData=data.frame(test=rep(1, 10)))
hasColDataColumns(SE, 'test')

}
