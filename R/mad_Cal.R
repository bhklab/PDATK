#' Parallelize calculation of the row-wise median absolute deviance
#'
#' @param dataset
#'
#' @importFrom SparkR dense_rank
#' @importFrom BiocParllel bpplapply
#' @export
calculateMedAbsDev <- function(dataset){

  madValues <- lapply(t(dataset), 
                      function(row) mad(as.numeric(row), na.rm=TRUE))
}



mad_cal= function(dataset){
  
  mad_values=vector("numeric", nrow(dataset))
  a=Sys.time()
  for(i in 1: length(rownames(dataset))){
    mad_values[i] = mad(as.numeric(dataset[i,]), na.rm=TRUE)
  }
  b=Sys.time()
  b - a
  cc=as.data.frame(cbind(genes=rownames(dataset), mad_values))
  #cc=cc[order(cc$mad_values, decreasing = TRUE),]
  cc$rank=dense_rank(-mad_values)
  #rownames(cc)=cc$genes
  #cc=cc[sort(rownames(cc)),]
  return(cc)
}
