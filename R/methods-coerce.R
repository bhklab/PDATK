#' @inherit setAs
#' 
#' @md
#' @export
setAs('matrix', 'data.frame', function(x) as.data.frame(x, row.names=rownames(x)))

#' @inherit setAs
#' 
#' @importFrom data.table as.data.table
#' 
#' @md
#' @export
setAs(from='matrix', to='data.table', 
    function(x) as.data.table(x, keep.rownames='rownames'))