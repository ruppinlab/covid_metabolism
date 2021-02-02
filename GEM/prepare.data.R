library(my.utils)
library(gembox)
library(Seurat)
data("recon1")

# prepare data for genome-scale metabolic modeling (GEM) analysis

load("../expression/de.and.gsea.res.RData")

exprs.int <- list()
dflux.int <- list()


# 1. Blanco-Melo et al. (the NHBE, A549-ACE2, and Calu-3 samples included in this study were used)

dat <- readRDS("../data/blanco.melo.RDS")$tpm

tmpf <- function(ctrl, trt, dat=dat) {
  ctrl <- exprs2int(recon1, dat[, colnames(dat)==ctrl])
  trt <- exprs2int(recon1, dat[, colnames(dat)==trt])
  list(ctrl=ctrl, trt=trt)
}

exprs.int$NHBE <- tmpf("Mock treated NHBE cells","SARS-CoV-2 infected NHBE cells")
exprs.int$`Calu-3` <- tmpf("Mock treated Calu-3 cells","SARS-CoV-2 infected Calu-3 cells")
exprs.int$A549 <- tmpf("Mock treated A549 cells trasnduced with a vector expressing human ACE2","SARS-CoV-2 infected A549 cells trasnduced with a vector expressing human ACE2")

dflux.int <- lapply(de.res[c("NHBE","Calu-3","A549")], function(x) {
  de.dt2dflux(recon1, x[, .(id, log.fc=-log.fc, pval, padj)], padj.cutoff=0.1)
}


# 2. Riva et al. (Vero E6 cell)

mat <- readRDS("../data/riva.RDS")$tpm
exprs.int$Vero <- tmpf("ctrl","trt",mat[,-4]) # 247 genes not in data
dflux.int$Vero <- de.dt2dflux(recon1, de.res$Vero[, .(id, log.fc=-log.fc, pval, padj)], topn=40, padj.cutoff=0.01)


# 3. Weingarten-Gabby et al. (the HEK293T cell samples included in this study were used)

dat <- readRDS("../data/weingarten.RDS")$count
tmp <- c("12h","18h","24h")
v.293 <- rowMeans(sapply(tmp, function(x) rowMeans(dat$expr[, dat$pheno[, cell=="HEK293T" & time==x]], na.rm=TRUE)))
c.293 <- rowMeans(dat$expr[, dat$pheno[, cell=="HEK293T" & treatment=="ctrl"]], na.rm=TRUE)
tmp <- list(ctrl=c.293, trt=v.293)
exprs.int$`293T` <- lapply(tmp, exprs2int, model=recon1)
dflux.int$`293T` <- de.dt2dflux(recon1, de.res$`293T`[, .(id, log.fc=-log.fc, pval, padj)], topn=100)


# 4. Bojkova et al. (Caco-2 cell), the 24h DE result provided by the authors was used

tmp <- fread("../data/bojkova.tsv")
dat1 <- data.matrix(tmp[,3:26])
rownames(dat1) <- gns
colnames(dat1) <- c(rep(c("ctrl2","ctrl6","ctrl10","ctrl24","v2","v6","v10","v24"), each=3))
exprs.int$`Caco-2` <- tmpf("ctrl24","v24",dat1)
dflux.int$`Caco-2` <- de.dt2dflux(recon1, de.res$`Caco-2`[, .(id, log.fc=-log.fc, pval, padj)], padj.cutoff=0.1)


# 5. Butler et al. (human patient swab samples), the DE result "Voom:Positive_vs_Negative:10M_samples:sva_correction_2sv" provided by the authors was used

x <- readRDS("../butler.RDS")
a <- x$baseMean - 173/580*x$log2FoldChange
b <- a + x$log2FoldChange
names(a) <- x$Gene
names(b) <- x$Gene
tmp <- list(ctrl=a, trt=b)
exprs.ave$Swab.Butler <- lapply(tmp, exprs2int, model=recon1)
dflux.int$Swab.Butler <- de.dt2dflux(recon1, x[, .(id, log.fc=-log.fc, pval=pval, padj)], topn=100, padj.cutoff=0.05)


# 6. Lieberman et al. (human patient swab samples)

dat <- readRDS("../data/lieberman.RDS")$log.tpm
tmpf <- function(x) {
  tmp <- c("C","R","T")
  ctrl <- rowMeans(sapply(tmp, function(i) rowMeans(x$expr[, x$pheno[, batch==i & virus=="neg"], drop=FALSE], na.rm=TRUE)), na.rm=TRUE)
  pat <- rowMeans(sapply(tmp, function(i) rowMeans(x$expr[, x$pheno[, batch==i & virus=="pos"], drop=FALSE], na.rm=TRUE)), na.rm=TRUE)
  res <- list(ctrl=ctrl, trt=pat)
}
exprs.int$Swab.Lieberman <- lapply(tmpf(dat), exprs2int, model=recon1)
dflux.int$Swab.Lieberman <- de.dt2dflux(recon1, de.res$Swab.Lieberman[, .(id, log.fc=-log.fc, pval, padj)], padj.cutoff=0.4, topn=60)


# 7. Xiong et al. (human patient BALF samples), the DE result provided by the authors was used

mat <- fread("../data/xiong.rpkm.tsv")
mat.ctrl <- data.matrix(mat[,11:13])
mat.trt <- data.matrix(mat[,7:10])
rownames(mat.ctrl) <- mat$Name
rownames(mat.trt) <- mat$Name
mat <- list(ctrl=mat.ctrl, trt=mat.trt)
exprs.int$BALF <- lapply(mat, function(x) exprs2int(recon1, x))
dflux.int$BALF <- de.dt2dflux(recon1, de.res$BALF[, .(id, log.fc=-log.fc, pval, padj)], topn=40)


# 8. Liao et al. (human patient single-cell data), the epithelial cells were used

dat <- readRDS("../data/liao.RDS") # the data has been pre-filtered and normalized
dat <- NormalizeData(dat)

cell.info <- fread("../data/liao.cell.annotation.txt")
dat$celltype <- cell.info[match(dat$ID, ID), celltype]
Idents(dat) <- paste(dat$celltype, dat$disease, sep="_")
tmp <- c(ctrl="Epithelial_N", patient="Epithelial_Y")
expr.ave <- lapply(tmp, function(x) {
  dat1 <- subset(dat, subset=cg1==x)
  res <- AverageExpression(dat1)$RNA[[1]]
  names(res) <- rownames(dat1)
  res
})
exprs.int$SC.Liao <- lapply(expr.ave, exprs2int, model=recon1)
dflux.int$SC.Liao <- de.dt2dflux(recon1, de.res$SC.Liao[, .(id, log.fc=-log.fc, pval, padj)], topn=80)


# 9. Chua et al. (human patient single-cell data), the basal and ciliated cells were used

dat <- readRDS("../data/chua.RDS")
dat <- NormalizeData(dat)

dat$cg1 <- paste(dat$celltype, dat$infection, sep="_")
dat$cg2 <- paste(dat$celltype, dat$severity, sep="_")

Idents(dat) <- dat$cg1
tmp <- c(ctrl="Basal_healthy", trt="Basal_SARS-CoV-2")
expr.ave <- lapply(tmp, function(x) {
  dat1 <- subset(dat, subset=cg1==x)
  res <- AverageExpression(dat1)$RNA[[1]]
  names(res) <- rownames(dat1)
  res
})
exprs.int$SC.Chua.Basal <- lapply(expr.ave, exprs2int, model=recon1)
dflux.int$SC.Chua.Basal <- de.dt2dflux(recon1, de.res$SC.Chua.Basal[, .(id, log.fc=-log.fc, pval, padj)], topn=40)

# for Ciliated cells, infected vs control failed to run, so used critical vs control instead
Idents(dat) <- dat$cg2
tmp <- c(ctrl="Ciliated_control", patient="Ciliated_critical")
expr.ave <- lapply(tmp, function(x) {
  dat1 <- subset(dat, subset=cg2==x)
  res <- AverageExpression(dat1)$RNA[[1]]
  names(res) <- rownames(dat1)
  res
})
exprs.int$SC.Chua.Ciliated<- lapply(expr.ave, exprs2int, model=recon1)
dflux.int$SC.Chua.Ciliated <- de.dt2dflux(recon1, de.res$SC.Chua.Ciliated[, .(id, log.fc=-log.fc, pval, padj)], topn=100)



res <- list(exprs.int=exprs.int, dflux.int=dflux.int)
saveRDS(res, file="data.for.gem.RDS")

