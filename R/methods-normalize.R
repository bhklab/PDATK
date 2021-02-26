#' Normalize the `assays` in a `SummarizedExperiment` Object
#'
#' @param object A `SummarizedExperiment` object with assays to normalize.
#' @param MARGIN An `integer` indicating if rows (1) or columns (2) should be 
#'   normalized. Defaults to columns. Defaults to 2.
#' @param FUN A function to normalize your data with. The function should accept
#'  a matrix, normalize the columns then return a matrix. The data will be 
#'  transposed before applyying FUN if MARGIN=1.
#' @param ... Fall through parameters to FUN. When using the default FUN, the 
#'   data is scaled and centered when no additional arguments are specified.
#' @param whichAssays A `numeric` or `character` vector specifying the indices
#'   or names of the assays to normalize. Defaults to all assays.
#' 
#' @return The `SummarizedExperiment` with one or more of the matrices in 
#'   `assays` normalized and the normalization details in the `normalization` 
#'   item of the object `metadata`.
#' 
#' @details 
#' When using the default FUN, it is also possible to impute missing values. 
#'   See `?caret::preProcess` for information on available methods.
#' 
#' @seealso [`PDATK::preprocessCaret`], [`caret::preProcess`]
#' 
#' @importMethodsFrom BiocGenerics normalize
#' 
#' @md
#' @export
setMethod('normalize', signature(object='SummarizedExperiment'), 
    function(object, MARGIN=2, FUN='preprocessCaret', ..., 
        whichAssays=seq_along(assays(object)))
{
    funContext <- .context(1)

    FUN_NAME <- as.symbol(as.character(substitute(FUN)))

    if (is.character(FUN)) FUN <- get(FUN)
    if (!is.function(FUN)) stop(.errorMsg(funContext, 'The argument to the FUN ',
        'parameter is not a function... Please ensure you pass a function or
        the name of a function to calculate row summaries with!'))
    
    .FUN <- if (MARGIN == 1) function(x, ...) { t(FUN(t(x), ...)) } else FUN

    assays(object)[whichAssays] <- endoapply(assays(object)[whichAssays], .FUN)
    metadata(object)$normalization <- list(FUN=substitute(FUN_NAME(x, ...)), 
        MARGIN=MARGIN, assays=whichAssays)
    return(object)
})

#' Normalize the `assays` of a `MutliAssayExperiment` Object
#' 
#' For this method to work, there must be a `normalize` method defined
#'   for all classes of experiments in the `MultiAssayExperiment`
#' 
#' @param object A `SummarizedExperiment` object with assays to normalize.
#' @param MARGIN An `integer` indicating if rows (1) or columns (2) should be 
#'   normalized. Defaults to 2 for columns.
#' @param FUN A function to normalize your data with. Should accept a 
#'   rectangular object such as a `matrix`, `data.frame`, or `data.table` and 
#'   return an object of the same class with the data normalized using FUN.
#' @param ... Fall through parameters to FUN. For the default FUN, these are 
#'   passed to `caret::preProcess` to allow configuration of the normalization 
#'   method. Omitting any arguments with the default FUN will scale and center 
#'   the data.
#' @param whichAssays A `numeric` or `character` vector specifying the indices
#'   of the assays to normalize. Defaults to all assays.
#' 
#' @return The `MultiAssayExperiment` with one or more of the assays normalized 
#'   and information about the normalization method in the `normalization` item 
#'   of the object `metadata`.
#' 
#' @details 
#' When using the default FUN, it is also possible to impute missing values. 
#'   See `?caret::preProcess` for information on available methods.
#' 
#' @seealso [`PDATK::preprocessCaret`], [`caret::preProcess`]
#' 
#' @importMethodsFrom BiocGenerics normalize
#' @importFrom MultiAssayExperiment ExperimentList
#' 
#' @md
#' @export
setMethod('normalize', signature(object='MultiAssayExperiment'),
    function(object, MARGIN=2, FUN='preprocessCaret', ..., 
        whichAssays=seq_along(assays(object)))
{
    funContext <- .context(1)

    FUN_NAME <- as.symbol(as.character(substitute(FUN)))

    if (is.character(FUN)) FUN <- get(FUN)
    if (!is.function(FUN)) stop(.errorMsg(funContext, 'The argument to the FUN ',
        'parameter is not a function... Please ensure you pass a function or
        the name of a function to calculate row summaries with!'))

    experiments(object)[whichAssays] <- 
        endoapply(experiments(object)[whichAssays], FUN=normalize, 
            MARGIN, FUN, ...)
    metadata(object)$normalization <- list(FUN=substitute(FUN_NAME(x, ...)), 
        MARGIN=MARGIN, assays=whichAssays)
    return(object)
})

setClassUnion('data.frame_or_matrix', c('data.frame', 'matrix'))

#' Normalize a `data.frame` Object
#' 
#' @param object A `data.frame` object to normalize.
#' @param MARGIN An `integer` indicating if rows (1) or columns (2) should be 
#'   normalized. Defaults to 2 for columns.
#' @param FUN A function to normalize your data with. Should accept a 
#'   rectangular object such as a `matrix`, `data.frame`, or `data.table` and 
#'   return an object of the same class with the data normalized using FUN.
#' @param ... Fall through parameters to FUN. For the default FUN, these are 
#'   passed to `caret::preProcess` to allow configuration of the normalization 
#'   method. Omitting any arguments with the default FUN will scale and center 
#'   the data.
#' 
#' @return The `data.frame` normalized.
#' 
#' @md
#' @export
setMethod('normalize', signature(object='data.frame_or_matrix'), 
    function(object, MARGIN=2, FUN='proprocessCaret', ...)
{
    funContext <- .context(1)

    FUN_NAME <- as.symbol(as.character(substitute(FUN)))

    if (is.character(FUN)) FUN <- get(FUN)
    if (!is.function(FUN)) stop(.errorMsg(funContext, 'The argument to the FUN ',
        'parameter is not a function... Please ensure you pass a function or
        the name of a function to calculate row summaries with!'))

    .FUN <- if (MARGIN == 1) function(x, ...) { t(FUN(t(x), ...)) } else FUN

    return(.FUN(object, ...))
})

#' Normalize a S4 `DFrame` Object
#' 
#' @param object A `DFrame` or `DataFrame` object to normalize.
#' @param MARGIN An `integer` indicating if rows (1) or columns (2) should be 
#'   normalized. Defaults to 2 for columns.
#' @param FUN A function to normalize your data with. Should accept a 
#'   rectangular object such as a `matrix`, `data.frame`, or `data.table` and 
#'   return an object of the same class with the data normalized using FUN.
#' @param ... Fall through parameters to FUN. For the default FUN, these are 
#'   passed to `caret::preProcess` to allow configuration of the normalization 
#'   method. Omitting any arguments with the default FUN will scale and center 
#'   the data.
#' 
#' @return A normalized `DFrame` object. 
#' 
#' @md
#' @export
setMethod('normalize', signature(object='DFrame'),
    function(object, MARGIN=2, FUN='preprocessCaret', ...) 
{
    df <- as.data.frame(object)
    df <- normalize(df, MARGIN, FUN, ...)

    DF <- DataFrame(df)
    mcols(DF) <- mcols(object)
    metadata(DF) <- metadata(object)
    
    return(DF)
})


#' Preprocess Data Using the `caret::preProcess` method, then return the 
#'   normalized data using `predict`.
#' 
#' @param x The data to be normalized with `caret::preProcess` and 
#'   `caret::predict`.
#' @param ... Fall through parameters to `caret::preProcess`. This can be used 
#'   to apply a range of different preprocessing methods from that package.
#' 
#' @seealso [`caret::preProcess`], [`stats::predict`]
#' 
#' @importFrom caret preProcess
#' @importFrom stats predict
#' 
#' @md
#' @export
preprocessCaret <- function(x, ...) predict(preProcess(x, ...), x)