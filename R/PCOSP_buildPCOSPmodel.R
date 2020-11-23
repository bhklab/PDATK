#' Build PCOSP Model from input data
#'
#' Building Pancreatic cancer overall survival predictor using unique 89
#'   samples profiled using microarray and sequencing platform.
#'
#' @examples
#' # To ensure reproducible results
#' set.seed(1987)
#'
#' # Load the data
#' data(trainingCohorts)
#'
#' # Return the object
#' selectedModels <- buildPCOSPmodels(trainingCohorts, numModels=10, nthread=1)
#'
#' # OR Save the object to disk
#' buildPCOSPmodel(trainingCohort, numModels=10, saveDir=tempdir())
#'
#' @param trainingCohorts A named [`list`] of training cohorts to fit and select
#'   PCOSP models for.
#' @param seqCohort [`character`] The names of the cohorts containing sequencing
#'   data in trainingCohorts. All other cohorts are assumped to contain
#'   microarray data.
#' @param numModels [`integer`] The number of models to fit
#' @param saveDir [`character`] A path to a directory to save the model. If you
#'   exclude this the function will return the model object instead.
#' @param nthread [`integer`] The number of threads to parallelize across
#'
#' @return [`list`] Either returns the model object or, if \code{saveDir} is
#'   specified it saves to disk instead and returns the path.
#'
#' @section Warning: This function uses random numbers; remember to
#'   `set.seed()` before running to ensure reproducible results.
#'
#' @md
#' @export
buildPCOSPmodels <- function(trainingCohorts, seqCohort, numModels, nthread,
    saveDir)
{
    # Set number of threads to parallelize over
    if (!missing(nthread)) {
        ops <- options()
        options("mc.cores"=nthread)
        on.exit(options(ops))
    }

    # Extract cohorts from trainingCohorts
    sequenceCohort <- trainingCohorts[[seqCohort]]
    arrayCohort <- trainingCohorts[[which(!(names(trainingCohorts) %in%
        seqCohort))]]

    # Merged common ICGC seq and array trainingCohorts
    commonData <- mergeCommonData(sequenceCohort, arrayCohort)

    # Training the model on ICGC seq/array common samples cohort
    cohortMatrix <- convertCohortToMatrix(commonData)
    cohortMatrixGroups <- ifelse(commonData$OS >= 365, 1, 0)

    selectedModels <- .generateTSPmodels(cohortMatrix, cohortMatrixGroups,
                                         numModels)

    # Save to disk or return
    if (!missing(saveDir)) {
        saveRDS(selectedModels,
                file=paste0(file.path(saveDir, 'PCOSPmodels'), '.rds'))
        return(paste0('Saved model to ', saveDir))
    } else {
        return(selectedModels)
    }
}



