#' Split an S4 Object into Training and Testing Data
#'
#' Divide an object into training and testing data using random sammpling.
#'
#' @param object `S4` An S4 object to split into training and testing groups
#'   for use building machine learning models.
#' @param ... Allow additional parmaters to be defined for this generic.
#'
#' @export
setGeneric('splitTrainTest',
    function(object, ...) standardGeneric('splitTrainTest'))
#'
#' @param object A `SurvivalExperiment` object to split into training and
#'   testing data. This function adds items `trainSamples` and `randomSeed`
#'   to the `metadata` slot of the object.
#' @param proportionTrain `float` What proportion of samples should be used for
#'   the training data. Default is 0.7.
#' @param ... Force subsequent arguments to be named, preventing from accidental
#'   seed setting.
#' @param randomSeed `numeric` A random seed which will be set locally and reset
#'   at the end of function exectuion. If missing, no seed is set; use this
#'   if you have a global random seed already set.
#'
#' @return
#'
#' @importFrom S4Vectors metadata metadata<-
#' @export
setMethod('splitTrainTest', signature('SurvivalExperiment'),
    function(object, proportionTrain=0.7, ..., randomSeed)
{
    # set random seed locally
    if (!missing(randomSeed)) {
        oldSeed <- .Random.seed
        set.seed(randomSeed)
        on.exit({ .Random.seed <- oldSeed })
    } else {
        warning(.warnMsg(.context(), 'The randomSeed parameter was not ',
            'specified. Defaulting to 1234.'))
    }

    # get sample information
    numSamples <- ncol(object)
    sampleNames <- colnames(object)

    # sample training samples
    trainSamples <- sampleNames[sample(seq_len(numSamples),
        floor(numSamples * proportionTrain), replace=FALSE)]

    # set object metadata
    metadata(object)$trainSamples <- trainSamples
    metadata(object)$modelParams$randomSeed <-
        if (!missing(randomSeed)) randomSeed else 1234

    return(object)
})
#' @param object A `CohortList` object to split into training and testing data.
#'
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom S4Vectors mcols mcols<-
#' @export
setMethod('splitTrainTest', signature('CohortList'),
    function(object, proportionTrain=0.7, ..., randomSeed)
{
    # set random seed locally
    if (!missing(randomSeed)) {
        oldSeed <- .Random.seed
        set.seed(randomSeed)
        on.exit({ .Random.seed <- oldSeed })
    } else {
        warning(.warnMsg(.context(), 'The randomSeed parameter was not ',
            'specified. Capturing the current enviroments .Random.seed instead.'))
    }

    # Do the train test split on each SurvivalExperiment
    object <- endoapply(object, splitTrainTest,
        proportionTrain=proportionTrain, randomSeed=randomSeed)

    # Assign the training sample names to elementMetadata
    .getTrainingSamples <- function(x) list(metadata(x)$trainSamples)
    mcols(object) <-
        DataFrame(do.call(rbind, lapply(object, .getTrainingSamples)))
    colnames(mcols(object)) <- 'trainSamples'

    # Add randomSeed to metadtaa
    metadata(object)$randomSeed <- .Random.seed
    return(object)
})