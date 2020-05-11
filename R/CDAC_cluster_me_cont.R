

cluster_me <-function(dataset, maxK, distance, method){
  
  results = ConsensusClusterPlus(dataset,maxK=maxK,reps=1000,pItem=0.8,pFeature=1,clusterAlg=method,distance=distance,innerLinkage="complete",seed=1987,corUse="pairwise.complete.obs")
  #pairwise.complete.obs takes into account the NA by calculating the correlation within non-NA values
  
  Kvec = 2:maxK
  x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
  PAC = rep(NA,length(Kvec)) 
  names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK
  
  for(i in Kvec){
    M = results[[i]][["consensusMatrix"]]
    Fn = ecdf(M[lower.tri(M)])
    PAC[i-1] = Fn(x2) - Fn(x1)
  }#end for i
  
  optK = Kvec[which.min(PAC)]
  
  zz=table(results[[optK]]$consensusClass)
  
  
  clusters=results[[optK]][["consensusClass"]]
  dataset_centroid=sapply(unique(clusters), clust.centroid, t(dataset) , clusters)
  
  if(length(which(is.na(dataset_centroid))) >0){
       dataset_centroid =   NaRV.omit(dataset_centroid)
    }
  ret=list(optimumK=optK,cluster_table= zz, centroid_clusters= dataset_centroid, classes=clusters )
  return(ret)
}


clust.centroid = function(i, dat, clusters) {
  ind = (clusters == i)
  
  if(length( which(ind == TRUE))==1){
    dat[ind,]
  }else{
    colMeans(dat[ind,], na.rm=TRUE)
  }
}