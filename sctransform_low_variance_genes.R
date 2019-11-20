#load datasets
library(SeuratData)
library(ifnb.SeuratData)
library(cbmc.SeuratData)
library(panc8.SeuratData)
hbrain <- readRDS("/home/butlera/Projects/ZooBrain/processed_data/human.rds")
bm <- readRDS('/home/butlera/Projects/muir/seurat_objects/citeseq.rds')
objects <- list(ifnb,cbmc,ipanc,hbrain,bm)
means <- list()
vars <- list()

#run SCTransform
for(i in 1:length(objects)) {
  objects[[i]] <- SCTransform(objects[[i]],ncells = 3000)
  print(i)
}

#find genes that appear in all datasets
allgenes=rownames(objects[[1]])
for(i in 2:5) {
  allgenes=intersect(allgenes, rownames(objects[[i]]))
}

#calculate means and residual variances in each dataset
#for variance, the ranking is relative based on other genes in the dataset
rmeans=list()
rvars=list()
for(i in 1:5) {
  tmean <- objects[[i]][["SCT"]]@meta.features$sct.gmean
  names(tmean) <- rownames(objects[[i]][["SCT"]])
  tmean <- tmean[allgenes]

  tvar <- objects[[i]][["SCT"]]@meta.features$sct.residual_variance
  names(tvar) <- rownames(objects[[i]][["SCT"]])
  tvar <- tvar[allgenes]
  
  rmeans[[i]] <- tmean
  rvars[[i]] <- rank(tvar)/length(tvar)
}

mean_matrix <- matrix(unlist(rmeans),nrow = length(rmeans[[1]]))
var_matrix <- matrix(unlist(rvars),nrow = length(rvars[[1]]))
rownames(mean_matrix) <- rownames(var_matrix) <- allgenes

#lvg = low-variance-genes, based on the mean level of dispersion across all datasets
vmin <- apply(var_matrix,1,mean)
mmin <- apply(mean_matrix,1,mean)
df <- data.frame(vmin,mmin)
df <- df[order(df$vmin),]

lvg <- rownames(df)[1:1000]





