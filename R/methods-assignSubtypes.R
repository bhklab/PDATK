#' Assign Sample Subtypes to an S4 Object
#'
#' @param object An `S4` object containing a slot representing samples or patients.
#' @param subtypes A mapping to assign subtypes to the samples or patients in
#'   the object.
#'
#' @return object with subtypes assigned to the sample metadata and the
#'   isSubtyped metadata item set to TRUE.
#'
#' @md
#' @export
setGeneric('assignSubtypes',
    function(object, subtypes, ...) standardGeneric('assignSubtypes'))

#' Assign Subtype Annotations to a SurvivalExperiment Object
#'
#' @param object A `SurvivalExperiment` object where the subtype annotations
#'   will be added to the `colData` slot as the `subtype` column.
#' @param subtypes A `data.frame` with
#' @param ... Force subsequent arguments to be named. Not used.
#' @param sampleCol A `character` vector specifying the name of the column
#'   containing the sample names. These must match the colnames of the
#'   `SurvivalExperiment`. If the sample names are the rownames of the `subtypes`
#'   `data.frame` then set this parameter to 'rownames'. Defaults to 'sample_name'.
#' @param subtypeCol A `character` vector specifying the name of the subtype
#'   column in the `subtypes` `data.frame` object. Defaults to 'subtype'.
#'
#' @return The `SurvivalExperiment` with the subtypes in the `subtypes` column
#'   of the colData slot and a metadata item, `hasSubtypes`, set to TRUE.
#'
#' @md
#' @export
setMethod('assignSubtypes', signature(object='SurvivalExperiment',
    subtypes='data.frame'), function(object, subtypes, ...,
        sampleCol='sample_name', subtypeCol='subtype')
{
    subtypes$rownames <- rownames(subtypes)
    columnData <- colData(object)

    subtypeSamples <- subtypes[[sampleCol]]
    columnSamples <- colnames(object)

    sampleHasSubtypes <- columnSamples %in% subtypeSamples

    if (all(!sampleHasSubtypes)) warning(.warnMsg(.context(), 'No samples in the',
        ' column names of the SurvivalExperiment match the sampleCol of the ',
        'subtypes data.frame.'))

    if (!all(sampleHasSubtypes)) message(.warnMsg(.context(), 'The samples ',
        paste0(columnSamples[!sampleHasSubtypes], collapse=', '), ' are not ',
        'present in the subtypes data.frame. Their subtype will be NA.'))

    columnData <- merge(columnData, subtypes[, c(sampleCol, subtypeCol)],
        by.x='sample_name', by.y=sampleCol, all.x=TRUE)

    # Lose rownames in join; reassign them
    rownames(columnData) <- rownames(colData(object))

    colData(object) <- columnData
    metadata(object)$hasSubtypes <- TRUE

    return(object)
})

#' Assign Subtype Annotations to a SurvivalExperiment Object
#'
#' @inherit assignSubtypes,SurvivalExperiment,data.frame-method
#'
#' @param object A `CohortList`
#' @param subtypes A `list` of `data.frame` objects, one per cohort, with
#'   to subtypes to assign to the colData slot of cohorts with a matching name.
#'
#' @return The `CohortList` with the subtypes in the `subtypes` column
#'   of the colData slot and a metadata item, `hasSubtypes`, set to TRUE for
#'   each `SurvivalExperiment`.
#'
#' @md
#' @importFrom S4Vectors mendoapply
#' @importFrom IRanges DataFrameList
#' @export
setMethod('assignSubtypes', signature(object='CohortList',
    subtypes='list'), function(object, subtypes, ..., sampleCol='sample',
        subtypeCol='subtype')
{
    if (!all(names(object) %in% names(subtypes))) stop(.errorMsg(.context(),
        'The names of the subtypes list must match the names of the CohortList',
        'passed as the object argument.'))

    subtypedCohortList <- mendoapply(FUN=assignSubtypes,
        object=object, subtypes=subtypes[names(object)],
        MoreArgs=list(sampleCol=sampleCol, subtypeCol=subtypeCol))

    return(subtypedCohortList)
})