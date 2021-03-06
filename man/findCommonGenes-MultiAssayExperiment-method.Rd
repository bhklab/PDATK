% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods-findCommonGenes.R
\name{findCommonGenes,MultiAssayExperiment-method}
\alias{findCommonGenes,MultiAssayExperiment-method}
\title{Intersect Gene Names for All \code{experiments} in a \code{MultiAssayExperiment}}
\usage{
\S4method{findCommonGenes}{MultiAssayExperiment}(object)
}
\arguments{
\item{object}{A \code{MultiAssayExperiment} where rownames represent genes.}
}
\value{
A \code{character} vector of genes common to all \code{experiments}
in the \code{MutliAssayExperiment}.
}
\description{
Intersect Gene Names for All \code{experiments} in a \code{MultiAssayExperiment}
}
