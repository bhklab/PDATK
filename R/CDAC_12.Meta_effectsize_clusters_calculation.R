library(plyr)
library(effsize)

## META EFFECT SIZE ########################################################################################

load('../../data/common_genes_cohorts_new.RData')       
load('../data/clusters.RData')

names(clusters)
names(cohorts1)

meta_effect=vector('list', 21)
sd=vector('list', 21)



for(j in 1:length(clusters)){
  for (i in 1:3){
    z=sort(unique(clusters[[j]]$meta_classes))
    k=1:4
    if(i %in% z){
      s= setdiff(unique(clusters[[j]]$meta_classes),i)
      
      response=mapvalues(clusters[[j]]$meta_classes,s, rep(0,length(s)))
      
      meta_effect[[j]][[i]] = sapply(1:dim(cohorts1[[j]])[1], function(x) cohen.d(as.numeric(cohorts1[[j]][x,]),response,hedges.correction=TRUE)$estimate)
      sd[[j]][[i]]= sapply(1:dim(cohorts1[[j]])[1], function(x) sd(as.numeric(cohorts1[[j]][x,which(response!=0)])))
      
    }
    else{
      meta_effect[[j]][[i]] = NA
      sd[[j]][[i]]= NA
    }
  }
  print(names(cohorts1)[j])
}

n=10331  #no. of genes
mean_effect_size=vector('list', 3)

for(j in 1: 3){
  for(i in 1:n){
    
    val=sapply(1:21, function(x) meta_effect[[x]][[j]][i])
    wt=sapply(1:21, function(x) 1/sd[[x]][[j]][i])
    
    
    mean_effect_size[[j]][[i]] = weighted.mean(unlist(val),unlist(wt), na.rm = TRUE)
  } 
}


mm=matrix(unlist(mean_effect_size), ncol=3)
rownames(mm)=rownames(cohorts1$pcsi)
colnames(mm)=c("Metacluster1","Metacluster2","Metacluster3")
mm=data.frame(mm)

save(mm,file="../results/meta_effectsize.RData")  

#######

