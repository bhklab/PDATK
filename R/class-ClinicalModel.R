#' ClinicalModel Class Definition
#'
#' @inherit SurvivalModel-class
#'
#' @md
#' @include class-SurvivalModel.R
#' @export
.ClinicalModel <- setClass('ClinicalModel', contains='SurvivalModel')


#' Constructor for the ClinicalModel Class
#'
#' @param trainData A `SurvivalExperiment` or `CohortList` object to construct
#'   a clinical model using
#' @param formula A `formula` object or a `character` vector coercible to one.
#'   All columns specified in the formula must be in the colData slot of the
#'   all `SurvivalExperiment`s in trainData.
#' @param minDaysSurvived An `integer` specifying the minimum number of days
#'   required to be 'good' prognosis. Default is 365.
#' @param ... Force all subsequent parameters to be named. Not used.
#' @param randomSeed An `integer` randomSeed that was used to train the model.
#'   Users should specify this when initializing a model to ensure
#'   reproducibilty.
#'
#' @return A `ClinicalModel` object.
#'
#' @examples
#' data(sampleICGCmicro)
#' set.seed(1987)
#' clinicalModel <- ClinicalModel(sampleICGCmicro,
#'   formula='prognosis ~ sex + age + T + N + M + grade', randomSeed=1987)
#'
#' @md
#' @importFrom plyr is.formula
#' @importFrom CoreGx .errorMsg .warnMsg
#' @export
ClinicalModel <- function(trainData, formula, minDaysSurvived=365, ...,
    randomSeed)
{
    funContext <- .context(1)

    if (missing(randomSeed)) stop(.errorMsg(funContext, 'No random seed was ',
        'specied for your model. Please include the value used for set.seed ',
        'when training this model! This ensures other can reproduce your ',
        'results.'))

    survModel <- SurvivalModel(trainData, minDaysSurived=minDaysSurvived,
            randomSeed=randomSeed)

    # ensure the formula is formatted correctly and all columns are in the data
    if (!(is.formula(formula) || is.character(formula)))
        stop(.errorMsg(funContext, "The formula for a clinical model must ",
            "either be a formula object or a character vector coercible to ",
            "a formula object"))

    # add the formula to the metadata
    metadata(survModel)$modelParams <-
        c(metadata(survModel)$modelParams, list(formula=formula))

    # make the clinical model
    clinicalModel <- .ClinicalModel(survModel)
    return(clinicalModel)
}

#' @importFrom CoreGx .warnMsg
#' @export
setValidity('ClinicalModel', function(object) {

    formula <- as.formula(metadata(object)$modelParams$formula)
    formulaCols <- as.character(formula[seq(2, 3)])
    formulaCols <- unlist(strsplit(formulaCols,
        split='[\\s]*[\\+\\-\\~\\=\\*][\\s]*', perl=TRUE))
    hasFormulaCols <- formulaCols %in% colnames(colData(object))
    if (!all(hasFormulaCols)) {
        warning(.warnMsg(funContext, 'The columns ', formulaCols[!hasFormulaCols],
            ' are missing from the colData slot of the training data',
            'Please only specify valid column names in colData to the formula!'))
        FALSE
    } else {
        TRUE
    }
})
