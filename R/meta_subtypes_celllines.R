

meta_subtype_cellines <- function(dataset1){
  
  feat <- c(gg$gene1, gg$gene2)
  g1= gg$gene1
  g2= gg$gene2
  gp= gg$genepairs
  
  if(length(which(rownames(dataset1) %in% feat)) != length(unique(feat))){
  
 z1= which(gg$gene1 %in% setdiff(gg$gene1, rownames(dataset1)))
 z2= which(gg$gene2 %in% setdiff(gg$gene2, rownames(dataset1)))
  zz= c(z1,z2)
  
    g1= gg$gene1[-zz]
    g2= gg$gene2[-zz]
    gp=gg$genepairs[-zz]
     output1=output1[-zz,]
     binary_rf_model <- randomForest( t(output1), y1)
     feat <- c(g1, g2)
   }
  
  
  xx= dataset1[feat,]
  
  ll=lapply(1: dim(xx)[2], function(x)  ifelse(xx[g1,x] > xx[g2,x], 1,0))
  output <- matrix(unlist(ll), nrow = ncol(xx), byrow = TRUE)
  rownames(output) = colnames(dataset1)
  colnames(output) = gp
  
  aaa = predict(binary_rf_model, output, type="prob")
    pred_subtypes1 = as.character(predict(binary_rf_model, output))

  dd=data.frame(id= colnames(dataset1), subtypes=pred_subtypes1 , prob1= aaa[,1], prob2= aaa[,2], prob3= aaa[,3])
  
  return(dd)
}
