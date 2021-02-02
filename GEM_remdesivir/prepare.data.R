library(my.utils)
library(gembox)
data("recon1")

# prepare data for genome-scale metabolic modeling (GEM) analysis on the remdesivir treatment data

exprs <- readRDS("../data/remdesivir.RDS")$tpm

tmp <- unique(colnames(exprs))
names(tmp) <- tmp
exprs.int <- lapply(tmp, function(x) {
  exprs2int(recon1, exprs[, colnames(exprs)==x])
})

load("../expression/de.and.gsea.res.remdesivir.RData")

tmpf <- function(x, r=-1, n=40, p=0.05) {
  de.dt2dflux(recon1, de.res[, .(id, log.fc=log.fc*r, pval, padj)], topn=n, padj.cutoff=p)
}
dflux.int <- list(
  v2c=tmpf("c2v"),
  vd2c=tmpf("c2vd"),
  d2c=tmpf("c2d"),
  vd2d=tmpf("d2vd"),
  v2vd=tmpf("v2vd",1)
)

save(exprs.int, dflux.int, file="data.for.gem.RData")

