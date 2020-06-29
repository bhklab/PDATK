#' Read in .tsv files of expression data and calculate then return a list with
#'   a signature score for each expression file.
#'
#' @param dataDir A \code{character} vector containing the path to the
#'    data directory. We recommend using `file.path` to make your paths
#'    platform agnostic!
#' @param fileNames A \code{character} vector of file names (with extensions)
#'    to read from and extract the gene signature scores.
#' @param geneCoefFile A \code{character} vector with the file name of the
#'    gene coefficients.
#' @param signed A \code{logical} vector indicating whether the gene coefficients
#'    read in are signed or not. Default is FALSE
#' @param ... Fall through parameters to `sig.scores` function from `genefu`.
#'    Overrides the default settings.
#'
#' @return A \code{list} of signature scores objects (as returned by `sig.score`)
#'     per specified file, named based on the file name.
#'
#' @section warning: this function assumes that the first column of expression
#'    data contains sample ID/identifiers! If it does not, please rearrange
#'    the columns before reading the files.
#'
#' @importFrom genefu sig.score
#' @export
readGenefuFilesToSigScores <- function(dataDir, fileNames, geneCoefFile, signed=FALSE, ...) {

    geneCoefficients <- read.table(file.path(dataDir, geneCoefFile),
                                   sep="\t", header=TRUE)

    expressionData <- lapply(fileNames,
                             function(file, dataDir)
                                 read.table(file.path(dataDir, file),
                                            sep="\t", header=T),
                             dataDir=dataDir)

    if (!missing(...)) {
        sigScores <- lapply(expressionData,
                            function(data, geneCoef, ...)
                                genefu::sig.score(x=geneCoef, data=data, ...),
                            geneCoef=geneCoefficients,
                            ...=...
                            )
    } else {
        sigScores <- lapply(expressionData,
                            function(data, geneCoef)
                                structure(genefu::sig.score(x=geneCoef, data=data,
                                          do.mapping=FALSE, mapping,
                                          size=0, cutoff=NA,
                                          signed=signed, verbose=FALSE
                                          )$score, .Names=data[,1]),
                            geneCoef=geneCoefficients)
    }

    names(sigScores) <- gsub("\\_.*txt|.txt", "", fileNames)
    return(sigScores)
}


# DEPRECATED --------------------------------------------------------------

#' Get a list of `sig.score` vectors named according to the probe they were
#'     calculated from.
#'
#' @param sigScores A \code{list} of signature scores caluclated with `sig.score`.
#'     As returned from the `readGenefuFilesToSigScores` function.
#'
#' @return A \code{list} of sig.score vectors with probe names as their names.
#'     The returned list has the same names as the sigScores list.
#'
#' @export
extractSigScoresProb <- function(sigScores) {
    lapply(sigScores,
           function(sigScore) structure(sigScore$score,
                                        .Names=rownames(sigScore$probe)))
}