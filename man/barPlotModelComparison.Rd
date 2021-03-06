% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods-barPlotModelComparison.R
\name{barPlotModelComparison}
\alias{barPlotModelComparison}
\title{Make A Bar Plot Comparing Performance Between Two \code{S4} Objects Representing
Mathematical Models.}
\usage{
barPlotModelComparison(model1, model2, ...)
}
\arguments{
\item{model1}{An\code{S4} object containing results of a mathematical model}

\item{model2}{An \code{S4} object containing results of a different mathematical
model, but with the same or overlapping samples.}

\item{...}{Allow new parameters to be defined for this generic.}
}
\value{
A bar plot comparing some aspect of model1 and model2
}
\description{
Make A Bar Plot Comparing Performance Between Two \code{S4} Objects Representing
Mathematical Models.
}
\examples{
data(sampleCohortList)
data(sampleValPCOSPmodel)
data(sampleICGCmicro)

# Set parallelization settings
BiocParallel::register(BiocParallel::SerialParam())

# Setup the models
set.seed(1987)
clinicalModel <- ClinicalModel(sampleICGCmicro,
  formula='prognosis ~ sex + age + T + N + M + grade',
  randomSeed=1987)

# Train the models
trainedClinicalModel <- trainModel(clinicalModel)

# Make predctions
clinicalPredCohortList <- predictClasses(sampleCohortList[c('PCSI', 'TCGA')],
  model=trainedClinicalModel)

# Validate the models
validatedClinicalModel <- validateModel(trainedClinicalModel,
  valData=clinicalPredCohortList)

# Plot the comparison
modelCompBarPlot <- barPlotModelComparison(validatedClinicalModel,
 sampleValPCOSPmodel, stat='AUC')

}
