#' Extract the Training and/or Testing Data as Subsets From An S4 Object.
#'
#' @param object An `S4` object to extract training and testing data from.
#' @param extract A `character` vector specifying what to extract.
#'
#' @md
#' @export
setGeneric('extractTrainTest',
    function(object, extract, ...) standardGeneric('extractTrainTest'))
#'
#' Extract the Training and/or Testing Data as Subsets From a SurvivalExperiment
#'
#' @param extract A `character` vector specifying whether to extract the
#'   data labeled 'test', 'train' or to extract both as a named list.
#'
#' @return
#'
#' @md
#' @importFrom S4Vectors metadata
#' @export
setMethod('extractTrainTest', signature('SurvivalExperiment', 'character'),
    function(object, extract='both', ...)
{
    extract <- match.arg(extract, c('train', 'test', 'both'))

    train <- metadata(object)$trainSamples

    if (is.null(train)) object <- splitTrainTest(object, ...)

    test <- setdiff(rownames(object), train)

    switch(extract,
           'both'={
               modelData <- vector(mode='list', 2L)
               names(modelData) <- c('train', 'test')
               modelData$train <- object[, train]
               modelData$test <- object[, test]
           },
           'train'={ modelData <- object[, train]  },
           'test'={ modelData <- object[, test] })
    return(modelData)
})
#'
#' @export
setMethod('extractTrainTest', signature('CohortList', 'character'),
    function(object, extract='both', ...)
{
    extract <- match.arg(extract, c('train', 'test', 'both'))

    train <- metadata(object)$trainSamples

    if (is.null(train)) object <- splitTrainTest(object, ...)

    test <- setdiff(rownames(object), train)

    switch(extract,
           'both'={
               modelData <- vector(mode='list', 2L)
               names(modelData) <- c('train', 'test')
               modelData$train <- subset(object, select=train)
               modelData$test <- subset(object, select=test)
           },
           'train'={ modelData <- subset(object, select=train)  },
           'test'={ modelData <- subset(object, select=test) })

    return(modelData)
})