library(my.utils)
library(gembox)
load("../data/Recon3D.RData")
m <- recon3d

# prepare data for genome-scale metabolic modeling (GEM) analysis on the remdesivir treatment data

exprs <- readRDS("../data/remdesivir.RDS")$tpm

tmp <- unique(colnames(exprs))
names(tmp) <- tmp
exprs.int <- lapply(tmp, function(x) {
  exprs2int(m, exprs[, colnames(exprs)==x])
})

load("../expression/de.and.gsea.res.remdesivir.RData")

dflux.int <- de.dt2dflux(m, de.res$c2vd[, .(id, log.fc=-log.fc, pval, padj)], topn=100, padj.cutoff=1.1)

save(exprs.int, dflux.int, file="data.for.gem.RData")
