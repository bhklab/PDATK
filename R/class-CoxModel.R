#' @inherit S4Model
#' 
#' @title CoxModel-class
#' 
#' @description Fit a cox proportional hazards model (`survival::coxph`) to one 
#'   or more experiments inside a `MultiAssayExperiment` object.
#' 
#' @md
#' @keywords internal
#' @export
.CoxModel <- setClass('CoxModel', contains='S4Model')

#' @title CoxModel constructor
#' 
#' @description Build a `CoxModel` object from a `MultiAssayExperiment` of
#'   `SurvivalExperiment`s. Allows easy application of the 
#'   `survival::coxph` function to many `SurvivalExperiment` at once,
#'   assuming they share the `survivalPredictor` column.
#' 
#' @param object A `MultiAssayExperiment` containing only `SurvivalExperiment`
#'   objects.
#' @param survivalPredictor A `character` vector indicating the name of the
#'   one or more columns in the colData slot of each `SurvivalExperiment`
#'   to use for testing survival differences between different groups. Must
#'   be a valid column in the colData of ALL experiments.
#' 
#' @return A `CoxModel` object, with `object` in the trainData slot.
#' 
#' @md
#' @export
CoxModel <- function(object, survivalPredictor='metacluster_labels') {

    ## TODO:: Finish error handling
    funContext <- .context(1)

    if (!is(object, 'MultiAssayExperiment'))
        stop(.errorMsg(funContext, 'The object argument must be a ',
            'MultiAssayExperiment!'))

    isSurvivalExperiment <- vapply(experiments(object), FUN=is, 
        'SurvivalExperiment', FUN.VALUE=logical(1))
    if (!all(isSurvivalExperiment)) stop(.errorMsg(funContext, 'All ',
        'experiments in a CoxModel must be SurvivalExperiments. Please ',
        'call the SurvivalExperiment constructor on your SummarizedExperiments',
        ' before trying to create a CoxModel.'))

    .CoxModel(trainData=object,
        modelParams=SimpleList(list(
            survivalPredictor=survivalPredictor
        )),
        models=SimpleList(),
        validationStats=data.table(),
        validationData=SimpleList(),
        mcols=DataFrame(),
        metadata=list()
    )
}