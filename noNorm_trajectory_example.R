## No normalization showing trajectory is driven by Sequencing Depth.
setwd("~/Desktop/normjam/")


da1 = readRDS("pbmc3k.rds")
rdat <- data.matrix(GetAssayData(object = da1, slot = "counts"))

library(scater)
rdat <- SingleCellExperiment(assays = list(counts = rdat, logcounts=log(rdat+1)))
rdat <- calculateQCMetrics(rdat)
rdat <- runPCA(rdat)
set.seed(99)
rdat <- runUMAP(rdat)
plotReducedDim(rdat, use_dimred = "UMAP", colour_by="log10_total_counts")

library(scran)
rdat.norm2 <- computeSumFactors(rdat)
rdat.norm2 <- normalize(rdat.norm2)
rdat.norm2 <- runPCA(rdat.norm2)
set.seed(99)
rdat.norm2 <- runUMAP(rdat.norm2)

plotReducedDim(rdat, use_dimred = "UMAP", colour_by="log10_total_counts")
plotReducedDim(rdat.norm2, use_dimred = "UMAP", colour_by="log10_total_counts")

### Focus on a specific cluster to show the trajectory

subt <- rdat[,which(rdat@reducedDims$UMAP[,1] > -5 & rdat@reducedDims$UMAP[,2] > 10)]
dim(subt)
set.seed(99)
sdat <- runPCA(sdat)
sdat <- runUMAP(sdat)
plotReducedDim(sdat, use_dimred = "UMAP", colour_by="log10_total_counts")

library(scran)
sdat.norm2 <- rdat.norm2[,which(rdat@reducedDims$UMAP[,1] > -5 & rdat@reducedDims$UMAP[,2] > 10)]
sdat.norm2 <- runPCA(sdat.norm2)
set.seed(99)
sdat.norm2 <- runUMAP(sdat.norm2)
plotReducedDim(sdat.norm2, use_dimred = "UMAP", colour_by="log10_total_counts")

library(slingshot)
sdat.sling <- slingshot(sdat, reducedDim = 'UMAP')
XX = plotReducedDim(sdat, use_dimred = "UMAP", colour_by="log10_total_counts")
CURVE = data.frame(SlingshotDataSet(sdat.sling)@curves$curve1$s)
colnames(CURVE) <- c("x", "y")
CURVE <- CURVE[SlingshotDataSet(sdat.norm2.sling)@curves$curve1$ord,]
XX + geom_segment(CURVE, mapping = aes(x = x, xend = dplyr::lead(x), y = y, yend = dplyr::lead(y)), 
                  size = 0.5, color='black')




sdat.norm2.sling <- slingshot(sdat.norm2, reducedDim = 'UMAP')
XX = plotReducedDim(sdat.norm2, use_dimred = "UMAP", colour_by="log10_total_counts")
CURVE.norm = data.frame(SlingshotDataSet(sdat.norm2.sling)@curves$curve1$s)
colnames(CURVE.norm) <- c("x", "y")
CURVE.norm <- CURVE.norm[SlingshotDataSet(sdat.norm2.sling)@curves$curve1$ord,]
CURVE.norm <- unique(CURVE.norm)
XX + geom_segment(CURVE.norm, mapping = aes(x = x, 
                                          xend = dplyr::lead(x), 
                                          y = y, yend = dplyr::lead(y)), 
                  size = 0.5, color='black')
plot(reducedDims(sdat.norm2.sling)$UMAP,  asp = 1, pch = 16)
lines(SlingshotDataSet(sdat.norm2.sling), lwd = 3, col = 'black')





t1 <- sdat.sling$slingPseudotime_1

# for time, only look at the 100 most variable genes
Y <- (assays(sdat.sling)$logcounts)
Y <- Y[which(apply(Y, 1, function(x) sum(x!=0)) > 0),]
library(gam)
# fit a GAM with a loess term for pseudotime
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t1)
  suppressWarnings({
    tmp <- suppressWarnings(gam(z ~ lo(t1), data=d))
  })
  p <- summary(tmp)[3][[1]][2,3]
  p
})



t1 <- sdat.norm2.sling$slingPseudotime_1
# for time, only look at the 100 most variable genes
Y <- (assays(sdat.norm2.sling)$logcounts)
Y <- Y[which(apply(Y, 1, function(x) sum(x!=0)) > 0),]
library(gam)
# fit a GAM with a loess term for pseudotime
gam.pval.norm <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t1)
  suppressWarnings({
    tmp <- suppressWarnings(gam(z ~ lo(t1), data=d))
  })
  p <- summary(tmp)[3][[1]][2,3]
  p
})

gam.pval.adj <- p.adjust(gam.pval, method = 'fdr')
gam.pval.norm.adj <- p.adjust(gam.pval.norm, method = 'fdr')

sum(gam.pval.adj < .05)
sum(gam.pval.norm.adj < .05)

X = names(gam.pval.adj[gam.pval.adj < .05])
Y = names(gam.pval.adj[gam.pval.norm.adj < .05])

write.table(X, file="~/sig_nonorm.txt", quote=F, row.names = F, col.names = F)
write.table(Y, file="~/sig_withnorm.txt", quote=F, row.names = F, col.names = F)
# 
# ### Try to focus on really different depths...
# 
# subt <- rdat[,which(rdat@reducedDims$UMAP[,1] > -5 & rdat@reducedDims$UMAP[,2] > 10)]
# subt <- subt[,which(subt$log10_total_counts > 3.4 | subt$log10_total_counts < 3.1)]
# dim(subt)
# set.seed(99)
# sdat <- runPCA(sdat)
# sdat <- runUMAP(sdat)
# plotReducedDim(sdat, use_dimred = "UMAP", colour_by="log10_total_counts")
# 
# ## Doesn't break apart as two clusters.



