#' Merge two `SurvivalExperiments`, subsetting to shared rows and columns
#'
#' @param x A `SurvivalExperiment`.
#' @param y A `SurvivalExperiment`.
#' @param cohortNames An optional `character` vector specifying the a name for
#'   each `SurvivalExperiment`.
#'
#' @return A `SurvivalExperiment` object with merge data from x and y, and
#'   the assay from each in the `assays` slot.
#'
#' @examples
#' data(sampleICGCmicro)
#' survExp2 <- sampleICGCmicro
#' mergedSurvExp <- merge(survExp2, sampleICGCmicro,
#'   cohortNames=c('copyICGCmicro', 'ICGCmicro'))
#' mergedSurvExp
#'
#' @md
#' @importFrom SummarizedExperiment colData colData<- rowData rowData<-
# @importFrom CoreGx .warnMsg .errorMsg
#' @importFrom data.table data.table as.data.table merge.data.table rbindlist
#'   `:=` copy .N .SD fifelse merge.data.table transpose setcolorder
#' @export
setMethod('merge', signature('SurvivalExperiment', 'SurvivalExperiment'),
    function(x, y, cohortNames)
{
    cohortL <- CohortList(list(x, y))
    names(cohortL) <- cohortNames

    # ensure the SurivalExperiments have common genes
    commonGenes <- findCommonGenes(cohortL)
    actualGenes <- Reduce(intersect, lapply(cohortL, rownames))

    # subset to common genes
    if (!all(actualGenes %in% commonGenes)) {
        cohortL <- subset(cohortL, subset=commonGenes)
        warning(.warnMsg(.context(), 'The training cohorts did not have only ',
            'common genes. Subsetting to common genes...'))
    }

    # ensure the SurvivalExperiments have common samples
    commonSamples <- findCommonSamples(cohortL)
    actualSamples <- unique(Reduce(c, lapply(cohortL, colnames)))

    # subset to common samples
    if (!all(actualSamples %in% commonSamples)) {
        cohortL <- subset(cohortL, select=commonSamples)
        warning(.warnMsg(.context(), 'The training cohorts did not have only ',
                         'common samples. Subsetting to common samples...'))
    }

    # merge the rowData
    rowDataL <- lapply(cohortL, rowData)
    rowDTL <- lapply(rowDataL, as.data.table)
    if (!all(rowDTL[[1]] == rowDTL[[2]], na.rm=TRUE)) {
        sharedRowData <- Reduce(intersect, lapply(rowDataL, colnames))
        mergeBy <- function(x, y) merge.data.table(x, y, by=sharedRowData)
        rowDT <- Reduce(mergeBy, rowDTL)
    } else {
        rowDT <- rowDTL[[1]]
    }

    # merge the colData
    ## TODO:: NA in join columns break the merge, how can I fix this?
    colDataL <- lapply(cohortL, colData)
    colDTL <- lapply(colDataL, as.data.table)
    if (!all(colDTL[[1]] == colDTL[[2]], na.rm=TRUE)) {
        sharedColData <- Reduce(intersect, lapply(colDataL, colnames))
        mergeBy <- function(x, y) merge.data.table(x, y, by=sharedColData)
        colDT <- Reduce(mergeBy, colDTL)
    } else {
        colDT <- colDTL[[1]]
    }
    metadata <- lapply(cohortL, metadata)
    metadataDT <- data.table(metadata[[1]], metadata[[2]])
    colnames(metadataDT) <- names(metadata)

    assaysL <- lapply(cohortL, assays)
    .lengthGt1 <- function(x) length(x) > 1
    multiAssays <- unlist(lapply(assaysL, .lengthGt1))
    if (any(multiAssays)) {
        namesOfAssays <- mapply(paste, names(assaysL), lapply(assaysL, names),
            sep='.', SIMPLIFY=FALSE)
    } else {
        namesOfAssays <- list(metadata[[1]]$mDataType, metadata[[2]]$mDataType)
    }
    for (i in seq_along(assays)) {
        names(assaysL[[i]]) <- namesOfAssays[[i]]
    }

    assays <- Reduce(c, assaysL)

    SurvivalExperiment(assays=assays, colData=colDT, rowData=rowDT,
        metadata=metadata)

})