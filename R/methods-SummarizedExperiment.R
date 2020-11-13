#' Generic function for converting objects to the PCOSP cohort format
#'
#' @param object [any] An object to reformat as a PCOSP cohort
#' @param ... [pairlist] Allow additional parameters to be dedfined on this generic
#'
#' @return Sample x feature [data.frame] where values are measurements for each feature in each sample, two additional
#'     columns, OS and OS_Status, store survival in days and survival status (1 is deceased), respectively.
#'
#' @export
setGeneric('toPCOSP', function(object, ...) standardGeneric('toPCOSP'))

#' Method for converting SummarizedExperiment objects from MetaGxPancreas into
#'   the PCOSP cohort format
#'
#' @param object A [SummarizedExperiment] object from the MetaGxPancreas datasets; assumes the assay of interest
#'     is the first assay in the SE and that the colData of the SE has columns `days_to_death` and `vital_status`
#' @return Sample x feature [data.frame] where values are measurements for each feature in each sample, two additional
#'     columns, OS and OS_Status, store survival in days and survival status (1 is deceased), respectively.
#'
#' @export
setMethod('toPCOSP',
          signature(object='SummarizedExperiment'),
          function(object){
              print(metadata(object)$experimentData@name)
              data <- assays(object)[[1]]
              DF <- data.frame(t(data), row.names=colnames(data))
              DF$OS <- as.numeric(colData(object)$days_to_death)
              DF$OS_Status <- ifelse(colData(object)$vital_status == 'deceased', 1, 0)
              return(DF)
          })