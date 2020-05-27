#'
#'
#'
#'
#'
#'



aa=merge(gcsi_subtypes, ccle_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,6])))/dim(aa)[1]
aa[,1][which(as.character(aa[,2]) != as.character(aa[,6]))]
aa[which(as.character(aa[,2]) != as.character(aa[,6])),]

aa=merge( ccle_subtypes, ctrpv2_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,6])))/dim(aa)[1]
aa[,1][which(as.character(aa[,2]) != as.character(aa[,6]))]
aa[which(as.character(aa[,2]) != as.character(aa[,6])),]


aa=merge(gcsi_subtypes,  ctrpv2_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,6])))/dim(aa)[1]
aa[,1][which(as.character(aa[,2]) != as.character(aa[,6]))]
aa[which(as.character(aa[,2]) != as.character(aa[,6])),]

rownames( gcsi_subtypes)=gcsi_subtypes[,1]
rownames( ccle_subtypes)=ccle_subtypes[,1]
rownames( gdsc_subtypes)=gdsc_subtypes[,1]
rownames( ctrpv2_subtypes)=ctrpv2_subtypes[,1]
celline_subtypes=list(gcsi_subtypes=  gcsi_subtypes, ccle_subtypes= ccle_subtypes, gdsc_subtypes=gdsc_subtypes, ctrpv2_subtypes= ctrpv2_subtypes)

save(celline_subtypes, file="../results/pharmacogx_celllines_subtypes.RData")
