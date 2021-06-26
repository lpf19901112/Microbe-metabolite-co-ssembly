#Modified by Pengfa Li.

library(picante)
library(ape)
library(parallel)

comun<-read.csv("otutable.csv",row.names=1,header=T)
phylo<-read.tree("phylo_tree.tre")
match.phylo.comun = match.phylo.data(phylo, comun)
str(match.phylo.comun)
print(c(date(),rep))
beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.comun$data),cophenetic(match.phylo.comun$phy),abundance.weighted=T));
dim(beta.mntd.weighted);
print(c(date(),rep))
head(beta.mntd.weighted)[1:5,1:5]
identical(colnames(match.phylo.comun$data),colnames(beta.mntd.weighted)); # just a check, should be TRUE
identical(colnames(match.phylo.comun$data),rownames(beta.mntd.weighted)); # just a check, should be TRUE

library(foreach)
library(doParallel)
detectCores(logical = F)
no_cores <- detectCores(logical = F)
no_cores
cl <- makeCluster(no_cores-1)
registerDoParallel(cl)

results <- foreach(x=1:999,.packages = "picante") %dopar% 
  as.matrix(comdistnt(t(match.phylo.comun$data),taxaShuffle(cophenetic(match.phylo.comun$phy)),abundance.weighted=T,exclude.conspecifics = F))

stopCluster(cl)

beta.reps <- 999
rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.comun$data),ncol(match.phylo.comun$data),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (l in 1:length(results)) {
  rand.weighted.bMNTD.comp[,,l] <- results[[l]]
}

weighted.bNTI = matrix(c(NA),nrow=ncol(match.phylo.comun$data),ncol=ncol(match.phylo.comun$data));
dim(weighted.bNTI);
weighted.bNTI

for (columns in 1:(ncol(match.phylo.comun$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.comun$data)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(match.phylo.comun$data);
colnames(weighted.bNTI) = colnames(match.phylo.comun$data);
weighted.bNTI;
write.csv(weighted.bNTI,"bNTI.csv",quote=F)
write.csv(beta.mntd.weighted,"bMNTD.csv",quote=F);


