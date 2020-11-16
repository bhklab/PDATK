#' Pancreatic Cancer Overall Survival Predictor (PCOSP) Class
#'
#'
#'
#' @slot rawdata `data.table` A sample by gene gene expression matrix
#' @slot model `` A train PCOSP model
#' @slot .intern
#'
#' @md
#' @export
setClass("PCOSP", slots=list(rawdata='matrix', model='',
    .intern='environment'))

#' Pancreatic Cancer Overall Survival Predictor (PCOSP) Constructor
#'
#' @md
#' @export
PCOSP <- function() {

}