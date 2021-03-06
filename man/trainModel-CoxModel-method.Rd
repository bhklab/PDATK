% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods-trainModel.R
\name{trainModel,CoxModel-method}
\alias{trainModel,CoxModel-method}
\title{Fit Models to the trainData in a CoxModel Object}
\usage{
\S4method{trainModel}{CoxModel}(object)
}
\arguments{
\item{object}{A \code{CoxModel} object to fit models for.}
}
\value{
A \code{CoxModel} object with the results of \code{coxph}, \code{survfit} and
\code{survdiff} in the models slot as lists where each item corresponds to
the data in \code{modelData}. For convenience, all the model data has also
been merged into a single \code{data.table} in the \code{modelDT} item of \code{models}.
}
\description{
Computes models with the \code{survival} package for \code{coxph}, \code{survfit},
\code{survdiff} as well as computes the fit p-values using \code{pchisq} with the
chisq values from \code{survdiff}. Modelling data is stored in \code{modelData},
as well as a \code{data.table} with all model data merged in \code{modelDT}. These
items are all assigned to the \code{models} slot.
}
