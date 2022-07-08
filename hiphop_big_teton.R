#run hiphop
library(hiphop)
args <- commandArgs(trailingOnly=TRUE)
index <- args[1]
cat("index =",index, fill=T)
geno <- as.matrix(read.csv(paste0("/project/ysctrout/wrosenth/NFS_19/parentage/hiphop/files/","ind",index,"_geno_0.9_0.15_0.01_ld_filt.csv"),header=T))
cat("geno dim = ",dim(geno),fill=T)
ind <- read.csv(paste0("/project/ysctrout/wrosenth/NFS_19/parentage/hiphop/files/","ind",index,"_ind_old_sequoia.csv"))
cat("ind dim =",dim(ind),fill=T)
ind[1,]
rownames(geno) <- ind$individual

test <- hothiphop(ind,geno)
test_results <- topmatch(test,ranking=c("hot.parents","hot.sire"),top=3)
test_results$offspring <- apply(test_results,1,function(x) ind[as.numeric(x[3]),2])
test_results$dam <- apply(test_results,1,function(x) ind[as.numeric(x[5]),2])
test_results$sire <- apply(test_results,1,function(x) ind[as.numeric(x[6]),2])

test_results2 <- topmatch(test,ranking=c("hothiphop.parents","hot.parents"),top=3)
test_results2$offspring <- apply(test_results2,1,function(x) ind[as.numeric(x[3]),2])
test_results2$dam <- apply(test_results2,1,function(x) ind[as.numeric(x[5]),2])
test_results2$sire <- apply(test_results2,1,function(x) ind[as.numeric(x[6]),2])

write.csv(test_results,file=paste0("./hiphop_results_hotparents_0.9_0.15_0.01_ld_filt_ind",index,".csv"),quote=F,row.names=F)
write.csv(test_results2,file=paste0("./hiphop_results_hothiphop_0.9_0.15_0.01_ld_filt_ind",index,".csv"),quote=F,row.names=F)