#' Return the name of the function and the name of the package that function
#'   is in when called within an R function.
#'
#' For providing context in user messages, warnings and errors
#'
#' @param n `integer` How far up the call stack to look for context. Defaults to
#'   2 since it is assumed this function will be used inside of `message`,
#'   `warning` or `stop`.
#'
#' @return A `character` vector with the name of the function
#'   `.getExecutionContext` was called from, as well as the package name,
#'   if applicable.
#'
#' @md
#' @keywords internal
#' @importFrom rlang trace_back
#' @importFrom utils packageName
#' @noRd
#' @aliases .context
.getExecutionContext <- function(n=2) {

    # name of function which called this function
    callStack <- rlang::trace_back()$calls
    context <- deparse(callStack[[length(callStack) - n]][1])

    # remove function arguments
    context <- gsub('\\(.*\\)', '', context)

    # deal with getting function names from inside an lapply statement
    ## TODO:: clean this up
    if (grepl('.*lapply.*', context)) {
        context <- deparse(callStack[[length(callStack) - (n + 1)]][3])
        context <- gsub('\\(.*\\)', '', context)
        # deal with S4 lapply calls (e.g., endoapply)
        if (grepl('.*match.fun.*', context)) {
            context <- deparse(callStack[[length(callStack) - (n + 5)]][3])
            context <- gsub('\\(.*\\)', '', context)
        }
    } else if (grepl('.*mapply.*', context)) {
        context <- deparse(callStack[[length(callStack) - (n + 1)]][1])
        context <- gsub('\\(.*\\)', '', context)
        if (grepl('.*match.fun.*', context)) {
            context <- deparse(callStack[[length(callStack) - (n + 5)]][1])
            context <- gsub('\\(.*\\)', '', context)
        }
    } else if (grepl('.*FUN.*', context)) {
        context <- deparse(callStack[[length(callStack) - (n + 2)]][3])
        context <- gsub('\\(.*\\)', '', context)
        # deal with S4 lapply calls (e.g., endoapply)
        if (grepl('.*match.fun.*', context)) {
            context <- tryCatch({ 
                deparse(callStack[[length(callStack) - (n + 6)]][3])
            }, error=function(e) 'context_failed')
            context <- gsub('\\(.*\\)', '', context)
        }
    }
    if (!grepl('::', context)) context <- paste0(packageName(), '::', context)

    return(paste0('\n[', context, '] '))
}
#' @noRd
.context <- .getExecutionContext

#' Convert factor columns in a rectangular object
#'
#' @param x A list-like rectangular object such as a `data.frame`,
#'   `data.table`, or `DataFrame`.
#'
#' @return `x` with factor columns converted to either integer or character,
#'   as appropriate.
#'
#' @examples
#' x <- data.frame(a=factor(LETTERS[1:5]), b=factor(runif(5, 0, 1)))
#' removeFactorColumns(x)
#'
#' @md
#' @export
removeFactorColumns <- function(x) {
    isFactor <- vapply(x, is.factor, logical(1))
    for (col in colnames(x)[isFactor]) {
        x[[col]] <-
            tryCatch({
                as.integer(levels(x[[col]])[x[[col]]]) },
            warning=function(w) { levels(x[[col]])[x[[col]]]
            })
    }
    return(x)
}

#' Remave any factor columns from the `colData` of an `S4` object
#'
#' @param x An `S4` object with a `colData` method defined for it.
#'
#' @return `x` with colData factor columsn converted to either integer or
#'   character, as appropriate.
#'
#' @section WARNING
#' If for some reason the factor levels are doubles, then this will coerce
#' them to integers.
#'
#' @examples
#' data(sampleICGCmicro)
#' removeColDataFactorColumns(sampleICGCmicro)
#'
#' @md
#' @export
removeColDataFactorColumns <- function(x) {
    colData(x) <- removeFactorColumns(colData(x))
    return(x)
}

#' Rename columns or do nothing if the names don't match
#'
#' @param x An object for which `colnames` is defined, probably a `data.frame`
#'   or other similar object.
#' @param value A character vector where names are the old column names
#'   and values are the new column names. Uses `gsub` internally to do
#'   the renaming.
#'
#' @return `x` with the updated column names if they are present. Does not
#'   fail if the column names are missing.
#'
#' @examples
#' x <- data.frame(a=factor(LETTERS[1:5]), b=factor(runif(5, 0, 1)))
#' renameColumns(x, c(a='c'))
#'
#' @md
#' @export
renameColumns <- function(x, values) {
    for (i in seq_along(values)) {
        colnames(x) <- gsub(names(values)[i], values[i], colnames(x))
    }
    return(x)
}

#'Rename the columns in the `colData` slot, or do nothing if they don't match
#'
#' @param x An `S4` object with a `colData` method.
#' @param values A character vector where names are the existing column
#'   names and values are the new column names.
#'
#' @return `x` with updated column names, if they match any existing columns.
#'
#' @examples
#' data(sampleICGCmicro)
#' renameColDataColumns(sampleICGCmicro, c(event_occurred='days_survived'))
#'
#' @md
#' @export
renameColDataColumns <- function(x, values) {
    colData(x) <- renameColumns(colData(x), values)
    return(x)
}

#' Check for Column Names in the colData Slot of a SummarizedExperiment
#' 
#' @param SE A `SummarizedExperiment` object to check for the existence
#'   of colData columns.
#' @param values A `character` vector with one or more column name to 
#'   check for in the column data.
#' 
#' @return `logical` True if all of values are column names in the 
#'   `SummarizedExperiment` object, FALSE otherwise.
#' 
#' @md
#' @export
hasColDataColumns <- function(SE, values) {
    all(values %in% colnames(colData(SE)))
}

#'
#' @param SE A `SummarizedExperiment` object to assign colData columns to
#' @param colname The name of the column to assign `values` to.
#' @param values The values to assign to the col
#' 
#' @md
#' @export
assignColDataColumn <- function(SE, colname, values) {
    funContext <- .context(1)
    if (length(values) != nrow(colData(SE || length(values != 1))))
        warning(.warnMsg(funContext), 'Assignment values not he same
            length as nrow(colData), values will be recycled!')
    colData(SE) <- within(colData(SE), 
        eval(str2expression(paste0(colname, ' <- ', capture.output(dput(values))))))
    return(SE)
}