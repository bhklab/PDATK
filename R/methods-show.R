#' Show method for a `CohortList` Object
#'
#' Implemented because the default show method was excluding list item names,
#'   even when there were valid names for the object.
#'
#' @md
#' @export
setMethod('show', signature('CohortList'), function(object)
{
    cat('< CohortList >\n')
    print(as(object, 'list'))
})