#' @title A Set of Example Patient Cohorts
#'
#' @name sampleCohortList
#'
#' @examples
#' data(sampleCohortList)
#'
#' @description
#' A `CohortList` object containing sample data for the PCOSP vignette. This
#' data is a subset of the the Pancreas datasets available in `MetaGxPancreas`.
#'
#' @seealso `MetaGxPancreas::loadPancreasDatasets`
#'
#' @md
NULL


#' @title A Sample SurvivalExperiment Containing Data from the ICGC micro-array
#'   cohort from `MetaGxPancreas`
#'
#' @examples
#' data(sampleICGCmicro)
#'
#' @name sampleICGCmicro
#'
#' @seealso `MetaGxPancreas::loadPancreasDatasets`
#'
#' @md
NULL


#' @title A Sample PCOSP Model Containing the ICGC micro-array cohort from
#'   `MetaGxPancreas` as training data.
#'
#' @examples
#' data(samplePCOSPmodel)
#'
#' @name samplePCOSPmodel
#'
#' @seealso `MetaGxPancreas::loadPancreasDatasets`
#'
#' @md
NULL


#' @title Sample RLS Model Containing the ICGC micro-array cohort from
#'   `MetaGxPancreas` as training data.
#'
#' @name sampleRLSmodel
#'
#' @examples
#' data(sampleRLSmodel)
#'
#' @seealso `MetaGxPancreas::loadPancreasDatasets`
#'
#' @md
NULL


#' @title Sample RGA Model Containing the ICGC micro-array cohort from
#'   `MetaGxPancreas` as training data.
#'
#' @name sampleRGAmodel
#'
#' @examples
#' data(sampleRGAmodel)
#'
#' @seealso `MetaGxPancreas::loadPancreasDatasets`
#'
#' @md
NULL


#' @title Sample ClinicalModel Containing the ICGC micro-array cohort from
#'   `MetaGxPancreas` as training data.
#'
#' @name sampleClinicalModel
#'
#' @examples
#' data(sampleClinicalModel)
#'
#' @seealso `MetaGxPancreas::loadPancreasDatasets`
#'
#' @md
NULL

#' @title Sample SurvivalExperiment Containing the PCSI rna-sequencing cohort
#'   from `MetaGxPancreas`.
#'
#' @name samplePCSIsurvExp
#'
#' @description
#' Used as validation data for modelling examples
#'
#' @seealso `MetaGxPancreas::loadPancreasDatasets`
#'
#' @examples
#' data(samplePCSIsurvExp)
#'
#' @md
NULL

#' existingClassifierData
#'
#' @examples
#' # Loads chen, birnbaum and haiderSigScores objects
#' data(existingClassifierData)
#'
#' Three objects:
#' - `chen`: The genes and coefficients for the gene signature from
#' Chen *et al.* (2015)
#' - `birnbaum`: The genes and coefficients for the gene signature from
#' Birnbaum *et al.* (2017)
#' - `haiderSigScores`: The classifier risk scores from Haider *et al.* (2014)
#'
#' @md
NULL


#' @title Published classifier gene signature for Chen
#'
#' @name chen
#'
#' @examples
#' # Loads chen, birnbaum and haiderSigScores objects
#' data(existingClassifierData)
#'
#' @description
#' The genes and coefficients for the gene signature from
#' Chen *et al.* (2015)
#'
#' @md
NULL


#' @title Published classifier gene signature for Birnbaum
#'
#' @name birnbaum
#'
#' @examples
#' # Loads chen, birnbaum and haiderSigScores objects
#' data(existingClassifierData)
#'
#' @description
#' The genes and coefficients for the gene signature from
#' Birnbaum *et al.* (2017)
#'
#' @md
NULL


#' @title Classifier survival scores for Haider
#'
#' @name haiderSigScores
#'
#' @examples
#' # Loads chen, birnbaum and haiderSigScores objects
#' data(existingClassifierData)
#'
#' @description
#' The classifier risk scores from Haider *et al.* (2014)
#'
#' @md
NULL


#' @title A list of sample subtypes for the data in sampleCohortList
#'
#' @name sampleCohortClasses
#'
#' @examples
#' data(sampleCohortClasses)
#'
#' @description
#' A `list` of `data.frames` containing clinical subtypes for the data
#'   in the `sampleCohortList` object.
#'
#' @md
NULL