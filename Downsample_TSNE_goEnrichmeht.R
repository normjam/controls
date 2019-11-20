bmSeu <- readRDS("~/Downloads/bm.cite.rds")
bmTypes <- bmSeu$celltype
library(Rtsne)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(BiocParallel)

bmCnt.raw <- bmSeu@assays$RNA@counts
bmTypes <- bmTypes[colSums(bmCnt.raw) >= 2000]
bmCnt.raw <- bmCnt.raw[, colSums(bmCnt.raw) >= 2000]
bmCnt.tpm <- 1e6 * t(t(bmCnt.raw) / colSums(bmCnt.raw))
subLvl <- 2e3
binVec <- seq(0, ncol(bmCnt.raw), 200)
binVec[length(binVec)] <- ncol(bmCnt.raw)
rawList <- lapply(1:(length(binVec) - 1), function(i) {
  bmCnt.raw[, (binVec[i] + 1):binVec[i + 1]]
})
bmCnt.sub <- do.call(cbind, bplapply(rawList, function(subDat) {
  Matrix(apply(subDat, 2, function(x) {
    if(sum(x) <= subLvl) {
      return(x)
    } else {
      cntVec <- unlist(sapply(1:length(x), function(i) {rep(i, x[i])}))
      cntSamp <- sample(cntVec, subLvl)
      x[1:length(x)] <- 0
      x[sort(unique(cntSamp))] <- table(cntSamp)
      return(x)
    }
  }))
}))

subAnn <- bmTypes %in% c("CD14 Mono", "CD4 Naive")

bmTSNE.raw <- Rtsne(t(log(as.matrix(bmCnt.raw) + 1)), partial_pca = TRUE, initial_dims = 25, 
                    max_iter = 2000, num_threads = 0, perplexity = 50)
bmMuVar.raw <- data.frame(
  mu = log(rowMeans(bmCnt.raw[, subAnn])), 
  s2 = log(apply(bmCnt.raw[, subAnn], 1, var))
)
plotDat.raw <- data.frame(
  x = bmTSNE.raw$Y[, 1],
  y = bmTSNE.raw$Y[, 2],
  depth = log(colSums(bmCnt.raw)),
  cell_type = bmTypes
)
pTSNE_depth.raw <- ggplot(plotDat.raw, aes(x = x, y = y, col = depth)) + 
  theme_classic() + 
  geom_point(size = 0.25) + 
  scale_color_viridis_c(option = "magma") + 
  labs(title = "Un-normalized") + 
  theme(
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 4), 
    title = element_text(size = 6), 
    axis.title = element_text(size = 6)
  )
pTSNE_depth.raw
pTSNE_type.raw <- ggplot(plotDat.raw, aes(x = x, y = y, col = cell_type)) + 
  theme_classic() + 
  geom_point(size = 0.25) + 
  labs(title = "Un-normalized", color = "cell\ntype") + 
  guides(colour = guide_legend(override.aes = list(size=1))) + 
  theme(
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 4), 
    title = element_text(size = 6), 
    axis.title = element_text(size = 6)
  )
pTSNE_type.raw
pMuVar.raw <- ggplot(bmMuVar.raw, aes(x = mu, y = s2)) + 
  theme_classic() + 
  geom_point(size = 0.5) + 
  labs(title = "Un-normalized", x = "Mean (log)", y = "Var (log)") + 
  theme(
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 4), 
    title = element_text(size = 6), 
    axis.title = element_text(size = 6)
  ) + 
  geom_abline(slope = 1, intercept = 0)
pMuVar.raw

bmTSNE.tpm <- Rtsne(t(log(as.matrix(bmCnt.tpm) + 1)), partial_pca = TRUE, initial_dims = 25, 
                    max_iter = 2000, num_threads = 0, perplexity = 50)
bmMuVar.tpm <- data.frame(
  mu = log(rowMeans(bmCnt.tpm[, subAnn])), 
  s2 = log(apply(bmCnt.tpm[, subAnn], 1, var))
)
plotDat.tpm <- data.frame(
  x = bmTSNE.tpm$Y[, 1],
  y = bmTSNE.tpm$Y[, 2],
  depth = log(colSums(bmCnt.raw)),
  cell_type = bmTypes
)
pTSNE_depth.tpm <- ggplot(plotDat.tpm, aes(x = x, y = y, col = depth)) + 
  theme_classic() + 
  geom_point(size = 0.25) + 
  scale_color_viridis_c(option = "magma") + 
  labs(title = "TPM") + 
  theme(
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 4), 
    title = element_text(size = 6), 
    axis.title = element_text(size = 6)
  )
pTSNE_depth.tpm
pTSNE_type.tpm <- ggplot(plotDat.tpm, aes(x = x, y = y, col = cell_type)) + 
  theme_classic() + 
  geom_point(size = 0.25) + 
  labs(title = "TPM", color = "cell\ntype") + 
  guides(colour = guide_legend(override.aes = list(size=1))) + 
  theme(
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 4), 
    title = element_text(size = 6), 
    axis.title = element_text(size = 6)
  )
pTSNE_type.tpm
pMuVar.tpm <- ggplot(bmMuVar.tpm, aes(x = mu, y = s2)) + 
  theme_classic() + 
  geom_point(size = 0.5) + 
  labs(title = "TPM", x = "Mean (log)", y = "Var (log)") + 
  theme(
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 4), 
    title = element_text(size = 6), 
    axis.title = element_text(size = 6)
  ) + 
  geom_abline(slope = 1, intercept = 0)
pMuVar.tpm

bmTSNE.sub <- Rtsne(t(log(as.matrix(bmCnt.sub) + 1)), partial_pca = TRUE, initial_dims = 25, 
                    max_iter = 2000, num_threads = 0, perplexity = 50)
bmMuVar.sub <- data.frame(
  mu = log(rowMeans(bmCnt.sub[, subAnn])), 
  s2 = log(apply(bmCnt.sub[, subAnn], 1, var))
)
plotDat.sub <- data.frame(
  x = bmTSNE.sub$Y[, 1],
  y = bmTSNE.sub$Y[, 2],
  depth = log(colSums(bmCnt.raw)),
  cell_type = bmTypes
)
pTSNE_depth.sub <- ggplot(plotDat.sub, aes(x = x, y = y, col = depth)) + 
  theme_classic() + 
  geom_point(size = 0.25) + 
  scale_color_viridis_c(option = "magma") + 
  labs(title = "Down-Sampled") + 
  theme(
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 4), 
    title = element_text(size = 6), 
    axis.title = element_text(size = 6)
  )
pTSNE_depth.sub
pTSNE_type.sub <- ggplot(plotDat.sub, aes(x = x, y = y, col = cell_type)) + 
  theme_classic() + 
  geom_point(size = 0.25) + 
  labs(title = "Down-Sampled", color = "cell\ntype") + 
  guides(colour = guide_legend(override.aes = list(size=1))) + 
  theme(
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 4), 
    title = element_text(size = 6), 
    axis.title = element_text(size = 6)
  )
pTSNE_type.sub
pMuVar.sub <- ggplot(bmMuVar.sub, aes(x = mu, y = s2)) + 
  theme_classic() + 
  geom_point(size = 0.5) + 
  labs(title = "Down-Sampled", x = "Mean (log)", y = "Var (log)") + 
  theme(
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 4), 
    title = element_text(size = 6), 
    axis.title = element_text(size = 6)
  ) + 
  geom_abline(slope = 1, intercept = 0)
pMuVar.sub

grid.arrange(pTSNE_depth.raw, pTSNE_type.raw, pMuVar.raw,
             pTSNE_depth.sub, pTSNE_type.sub, pMuVar.sub, 
             nrow = 2)
grid.arrange(pTSNE_depth.raw, pTSNE_type.raw, nrow = 1)
grid.arrange(pTSNE_depth.sub, pTSNE_type.sub, nrow = 1)
grid.arrange(pMuVar.raw, pMuVar.sub, nrow = 1)


#####################
## Sig. Var. Genes ##
#####################
posRaw <- !is.infinite(bmMuVar.raw$mu) & !is.infinite(bmMuVar.raw$s2)
mod.raw <- loess(s2 ~ mu, data = bmMuVar.raw[posRaw, ])
sigGenes.raw <- residuals(mod.raw) / sd(residuals(mod.raw)) > 2
sigGenes.raw <- rownames(bmCnt.raw)[posRaw][sigGenes.raw]

posSub <- !is.infinite(bmMuVar.sub$mu) & !is.infinite(bmMuVar.sub$s2)
mod.sub <- loess(s2 ~ mu, data = bmMuVar.sub[posSub, ])
sigGenes.sub <- residuals(mod.sub) / sd(residuals(mod.sub)) > 2
sigGenes.sub <- rownames(bmCnt.sub)[posSub][sigGenes.sub]


###########################
## Sig. Genes Enrichment ##
###########################
suppressMessages(library(topGO))
suppressMessages(library(org.Hs.eg.db))
goEnrich <- function(refGenes,selGenes,mapping){
  geneVec <- factor(as.integer(refGenes%in%selGenes))
  names(geneVec) <- refGenes
  GOdata <- new("topGOdata",ontology="BP",allGenes=geneVec,
                annot=annFUN.org,mapping=mapping,ID="Symbol",nodeSize=10)
  allGO = usedGO(object = GOdata)
  if(length(allGO)==0){return(NA)}
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  topResults <- GenTable(GOdata,classicFisher=resultFisher,topNodes=length(allGO))
  topResults$FDR <- p.adjust(topResults$classicFisher,method="bonferroni")
  return(topResults)
}

TopGo.raw <- goEnrich(rownames(bmCnt.raw), sigGenes.raw, "org.Hs.eg.db")
TopGo.sub <- goEnrich(rownames(bmCnt.sub), sigGenes.sub, "org.Hs.eg.db")
TopGo.diffRS <- goEnrich(rownames(bmCnt.sub), setdiff(sigGenes.raw, sigGenes.sub), "org.Hs.eg.db")
TopGo.diffSR <- goEnrich(rownames(bmCnt.sub), setdiff(sigGenes.sub, sigGenes.raw), "org.Hs.eg.db")

cleanGOdat <- function(GOdat,nodeTrim=100) {
  GOdat$hitPct <- GOdat$Significant/GOdat$Annotated
  GOdat <- GOdat[GOdat$Annotated>=nodeTrim,]
  GOdat <- GOdat[1:10, ]
  GOdat$Term <- unlist(lapply(strsplit(GOdat$Term, " "), function(x) {
    n <- cumsum(nchar(x))
    n <- n+seq_len(length(n))-1
    if(length(x)==1){return(x)}
    ret <- x[1]
    m <- 0
    for(i in 2:length(x)) {
      if(n[i-1] >= m + 14) {
        m <- n[i-1]
        ret <- c(ret,"\n",x[i])
      } else {
        ret <- c(ret," ",x[i])
      }
    }
    return(paste(ret,collapse=""))
  }))
  GOdat$Term <- factor(as.character(GOdat$Term), 
                       levels = as.character(GOdat$Term), 
                       ordered = TRUE)
  return(GOdat)
}
plotGO <- function(GOdat,nodeTrim=10){
  GOdat <- cleanGOdat(GOdat,nodeTrim)
  p <- ggplot(GOdat,aes(x=hitPct,y=forcats::fct_rev(Term),size=Annotated,color=FDR))+
    theme_linedraw()+
    geom_point()+
    scale_size(trans="log10")+
    scale_color_viridis_c(option = "magma") + 
    theme(plot.title=element_text(size=18,face="bold.italic"),
          axis.title.x=element_text(size=14,face="bold"),
          axis.title.y=element_text(size=14,face="bold"),
          legend.title=element_text(size=8,face="bold"),
          legend.text=element_text(size=8),
          axis.text=element_text(size=10),
          axis.ticks.y=element_blank(),
          strip.text.y=element_text(size=10,face="bold",angle=0),
          legend.position="left",
          plot.background=element_rect(fill="white", colour="white"))+
    labs(x="hits (%)",y="GO term",color="FDR",size="Count",title="Top terms")
  return(p)
}
p_TopGo.raw <- plotGO(TopGo.raw)
p_TopGo.raw <- plotGO(TopGo.raw)

