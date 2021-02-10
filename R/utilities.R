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
#' @noRd
#' @aliases .context
.getExecutionContext <- function(n=2) {

    # name of function which called this function
    callStack <- rlang::trace_back()$calls
    context <- deparse(callStack[[length(callStack) - n]][1])

    # remove function arguments
    context <- gsub('\\(.*\\)', '', context)
    print(context)

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
            context <- deparse(callStack[[length(callStack) - (n + 6)]][3])
            context <- gsub('\\(.*\\)', '', context)
        }
    }
    if (!grepl('::', context)) context <- paste0(packageName(), '::', context)

    return(paste0('\n[', context, '] '))
}
#' @noRd
.context <- .getExecutionContext