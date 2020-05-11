
msm_thresholds <- function(name_dataset1, name_dataset2, cluster1, cluster2, dataset,comp_no){
  
  map_rows = which(rownames(cluster1$centroid_clusters) %in% rownames(dataset))
  cluster1$centroid_clusters = cluster1$centroid_clusters[map_rows,]
  dataset=dataset[ rownames(cluster1$centroid_clusters),]
  set.seed(1987)
  m=vector()
  a=list()
  k=1
  for( i in 1: cluster1$optimumK){
    
    for(j in 1:cluster2$optimumK){
      ind = (cluster2$classes  == j)
      x=as.numeric( which(ind==TRUE))
      m[j]=mean(unlist(lapply(x, function(x) as.numeric(cor.test(cluster1$centroid_clusters[,i], dataset[,x], method="pearson")$estimate))))
      
    }}
return(m)
    
    
    }