% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class-S4Model.R
\name{models}
\alias{models}
\title{Accessor for the models slot of an \code{S4} object}
\usage{
models(object, ...)
}
\arguments{
\item{object}{An \code{S4} object to retrieve models from.}

\item{...}{Allow new parameters to be defined for this generic.}
}
\value{
An \code{S4} object representing a model.
}
\description{
Accessor for the models slot of an \code{S4} object
}
\examples{
data(CSPC_MAE)
set.seed(1987)
metaclustModel <- ConMetaclustModel(CSPC_MAE, randomSeed=1987)
models(metaclustModel)

}
