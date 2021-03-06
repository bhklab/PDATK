% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods-predictClasses.R
\name{predictClasses,SurvivalExperiment,ClinicalModel-method}
\alias{predictClasses,SurvivalExperiment,ClinicalModel-method}
\title{Predict Survival Prognosis Classes and Risk Scores for A \code{SurvivalModel} Using
a \code{ClinicalModel} Object.}
\usage{
\S4method{predictClasses}{SurvivalExperiment,ClinicalModel}(object, model, ..., na.action = "na.exclude", type = "response")
}
\arguments{
\item{object}{A \code{SurvivalExperiment} object with the correct columns in
\code{colData} to match the formula for the \code{ClinicalModel} object.}

\item{model}{A trained \code{ClinicalModel} object, as return by \code{trainModel}.}

\item{...}{Fall through parameters to \code{\link[stats:predict]{stats::predict}}.}

\item{na.action}{The \code{na.action} paramter passed to \code{\link[stats:predict.glm]{stats::predict.glm}}.}

\item{type}{The \code{type} parameter passed to \code{\link[stats:predict.glm]{stats::predict.glm}}}
}
\value{
A \code{SurvivalExperiment} with the model predictions in the colData
slot as clinical_prob_good.
}
\description{
Predict Survival Prognosis Classes and Risk Scores for A \code{SurvivalModel} Using
a \code{ClinicalModel} Object.
}
\examples{
data(sampleClinicalModel)
data(samplePCSIsurvExp)

# Set parallelization settings
BiocParallel::register(BiocParallel::SerialParam())

# Train Model
trainedClinicalModel <- trainModel(sampleClinicalModel)

# Make predictions
ClinicalPredSurvExp <- predictClasses(samplePCSIsurvExp,
  model=trainedClinicalModel)
head(colData(ClinicalPredSurvExp))

}
