#' @inherit setAs
#' 
#' @md
#' @export
setAs('matrix', 'data.frame', function(from) as.data.frame(from, row.names=rownames(from)))

#' @inherit setAs
#' 
#' @importFrom data.table as.data.table
#' 
#' @md
#' @export
setAs(from='matrix', to='data.table', 
    function(from) as.data.table(from, keep.rownames='rownames'))