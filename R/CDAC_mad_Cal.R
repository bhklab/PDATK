
mad_cal= function(dataset){
  mad_values=vector()
  
  for(i in 1: length(rownames(dataset))){
    mad_values[i] = mad(as.numeric(dataset[i,]), na.rm=TRUE)
  }
  cc=as.data.frame(cbind(genes=rownames(dataset), mad_values))
  #cc=cc[order(cc$mad_values, decreasing = TRUE),]
  cc$rank=dense_rank(-mad_values)
  #rownames(cc)=cc$genes
  #cc=cc[sort(rownames(cc)),]
  return(cc)
}


