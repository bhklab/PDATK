% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class-SurvivalModel.R
\name{validationStats,SurvivalModel-method}
\alias{validationStats,SurvivalModel-method}
\title{Accessor for the \code{validationStats} slot of a \code{SurvivalModel} object.}
\usage{
\S4method{validationStats}{SurvivalModel}(object)
}
\arguments{
\item{object}{A \code{SurvivalModel} object to get validation statistics from.}
}
\value{
A \code{data.table} of validation statistics for the \code{SurvivalModel}
object.
}
\description{
Accessor for the \code{validationStats} slot of a \code{SurvivalModel} object.
}
\examples{
data(samplePCOSPmodel)
validationStats(samplePCOSPmodel)

}
