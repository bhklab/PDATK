
compare_interclusters <- function(name_dataset1, name_dataset2, cluster1, cluster2, dataset,comp_no){
  
  print(paste(comp_no,"Comparing intraclusters between", name_dataset1, "and", name_dataset2,sep=" "))
  
  map_rows = which(rownames(cluster1$centroid_clusters) %in% rownames(dataset))
  cluster1$centroid_clusters = cluster1$centroid_clusters[map_rows,]
  dataset=dataset[ rownames(cluster1$centroid_clusters),]
  set.seed(1987)
  
  
  
  z=clusterRepro(cluster1$centroid_clusters, dataset, Number.of.permutations = 500)
  m=vector()
  a=list()
  k=1
  for( i in 1: cluster1$optimumK){
    
    for(j in 1:cluster2$optimumK){
      ind = (cluster2$classes  == j)
      x=as.numeric( which(ind==TRUE))
      m[j]=mean(unlist(lapply(x, function(x) as.numeric(cor.test(cluster1$centroid_clusters[,i], dataset[,x], method="pearson")$estimate))))
      
    }
    if(is.na(z$p.value[i]) != TRUE){
      if(z$Actual.IGP[i] >0.5 & z$p.value[i] < 0.05 & max(m) > 0){
        c1= paste(name_dataset1, i,sep="-")
        c2= paste(name_dataset2, which.max(m), sep="-")
        
        a[[k]] = list(c1,c2,max(m), z$Actual.IGP[i], z$p.value[i])
        
        k=k+1
      }
    }
  }
  
  return(a)
}
