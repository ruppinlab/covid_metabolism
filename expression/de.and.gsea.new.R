library(my.utils)
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(Seurat)

# differential expression and GSEA analysis of each dataset

### DE analysis

# all results will be saved in `de.res`
de.res <- list()


# 1. Blanco-Melo et al. (the NHBE, A549-ACE2, and Calu-3 samples included in this study were used)

gse <- readRDS("../data/blanco.melo.RDS")$gse

nrow(gse) # [1] 60233
keep <- rowSums(assays(gse)$counts>=10)>=17
sum(keep) # [1] 15066
gse1 <- gse[keep,]

de1 <- function(ctrl, trt) {
  dat <- gse1[, gse1$source_name %in% c(ctrl, trt)]
  dat$group <- factor(dat$source_name, levels=c(ctrl,trt))
  dds <- DESeqDataSet(dat, design = ~ group)
  dds <- DESeq(dds)
  de.res <- results(dds, contrast=c("group",trt,ctrl))
  tids <- rownames(de.res)
  gids <- substr(tids, 1, 15)
  de.res <- as.data.table(de.res)
  gss <- mapIds(org.Hs.eg.db, keys=gids, column="SYMBOL", keytype="ENSEMBL", multiVals="first")
  de.res <- de.res[, .(id=gss, ensembl=tids, ave=baseMean, log.fc=log2FoldChange, lfc.se=lfcSE, stat, pval=pvalue, padj=padj)]
  de.res <- de.res[order(padj, pval)]
}

de.res$NHBE <- de1("Mock treated NHBE cells","SARS-CoV-2 infected NHBE cells")
de.res$`Calu-3` <- de1("Mock treated Calu-3 cells","SARS-CoV-2 infected Calu-3 cells")
de.res$A549 <- de1("Mock treated A549 cells trasnduced with a vector expressing human ACE2","SARS-CoV-2 infected A549 cells trasnduced with a vector expressing human ACE2")


# 2. Riva et al. (Vero E6 cell)

mat <- readRDS("../data/riva.RDS")$count

idx <- rowSums(mat>=10)>=3
sum(idx) # [1] 13221
mat1 <- round(mat[idx,])

phe <- data.table(group=factor(rep(c("ctrl","trt"), each=3)))

dds1 <- DESeqDataSetFromMatrix(countData=mat1[,-4], colData=phe[-4], design=~group) # discared sample #4 as it seems to be an outlier
dds1 <- DESeq(dds1)
de.res1 <- results(dds1, contrast=c("group","trt","ctrl"))
gid <- rownames(de.res1)
de.res1 <- as.data.table(de.res1)[, .(gid=gid, ave=baseMean, log.fc=log2FoldChange, lfc.se=lfcSE, pval=pvalue, padj=padj)]

tmp <- str_match(gid, "[^_]*_(.*)")[,2]
idx <- is.na(tmp) | grepl("gp|rRNA|^MIR|Y_RNA|^U[0-9]*$|^SNOR|sno|^RNU|Metazoa|7SK|^SCARNA|RNase|Vault|Telomerase", tmp)
de.res1.sub <- cbind(data.table(id=tmp[!idx]), de.res1[!idx])
de.res1.sub[, padj:=p.adjust(pval, "BH")]
de.res$Vero <- de.res1.sub[order(padj, pval)]


# 3. Weingarten-Gabby et al. (the HEK293T cell samples included in this study were used)

dat <- readRDS("../data/weingarten.RDS")$count
dat1 <- rm.low.genes(dat, rm.low.frac.gt=15/16)
dat.293t <- filter.eset(dat1, cell=="HEK293T")
dds1 <- DESeqDataSetFromMatrix(countData=round(dat.293t$expr), colData=dat.293t$pheno, design=~treatment)
dds1 <- DESeq(dds1)
de.res1 <- results(dds1, contrast=c("treatment","virus","ctrl"))
gid <- rownames(de.res1)
de.res1 <- as.data.table(de.res1)[, .(id=gid, ave=baseMean, log.fc=log2FoldChange, lfc.se=lfcSE, pval=pvalue, padj=padj)][order(padj,pval)]
de.res$`293T` <- de.res1


# 4. Bojkova et al. (Caco-2 cell), the 24h DE result provided by the authors was used

tmp <- fread("../data/bojkova.tsv")
de.res$`Caco-2` <- tmp[, .(id=ifelse(`Gene Symbol`=="",`UniProt Accession`,`Gene Symbol`), log.fc=`Ratio 24h`, pval=`P value 24h`)][, padj:=p.adjust(pval,"BH")]


# 5. Butler et al. (human patient swab samples), the DE result "Voom:Positive_vs_Negative:10M_samples:sva_correction_2sv" provided by the authors was used

de.res$Swab.Butler <- readRDS("../butler.RDS")


# 6. Lieberman et al. (human patient swab samples)

dat <- readRDS("../data/lieberman.RDS")$count
dat$pheno[, table(virus, batch)]
dat1 <- filter.eset(dat, batch %in% c("C","R","T")) # I decided to use only batch C, R and T where there are both negative and positive samples
mat <- rm.low.genes(dat1$expr, rm.low.frac.gt=1-1/nrow(dat1$pheno))
phe <- copy(dat1$pheno)
phe[, batch:=C(as.factor(batch), "contr.sum")]
dds1 <- DESeqDataSetFromMatrix(countData=mat, colData=phe, design=~virus+batch)
dds1 <- DESeq(dds1)
de.res1 <- results(dds1, contrast=c("virus","pos","neg"))
gid <- rownames(de.res1)
de.res1 <- as.data.table(de.res1)[, .(id=gid, ave=baseMean, log.fc=log2FoldChange, lfc.se=lfcSE, pval=pvalue, padj=padj)][order(padj,pval)]
de.res$Swab.Lieberman <- de.res1


# 7. Xiong et al. (human patient BALF samples), the DE result provided by the authors (it was also done with DESeq2) was used

tmp <- fread("../data/xiong.de.tsv")
de.res$BALF <- setnames(tmp, c("Name","log2FoldChange","pvalue"), c("id","log.fc","pval"))


# 8. Liao et al. (human patient single-cell data), the epithelial cells were used

dat <- readRDS("../data/liao.RDS") # the data has been pre-filtered and normalized
dat <- NormalizeData(dat)

cell.info <- fread("../data/liao.cell.annotation.txt")
dat$celltype <- cell.info[match(dat$ID, ID), celltype]
Idents(dat) <- paste(dat$celltype, dat$disease, sep="_")
tmp <- tryCatch(FindMarkers(dat, ident.1="Epithelial_Y", ident.2="Epithelial_N", test.use="MAST", min.pct=0.01, logfc.threshold=0), error=function(e) NULL)
tmp <- cbind(id=row.names(tmp), as.data.table(tmp))
setnames(tmp, c("id","pval","log.fc","pct1","pct2","padj"))
de.res$SC.Liao <- tmp


# 9. Chua et al. (human patient single-cell data), the basal and ciliated cells were used

dat <- readRDS("../data/chua.RDS")
dat <- NormalizeData(dat)

dat$cg1 <- paste(dat$celltype, dat$infection, sep="_")
dat$cg2 <- paste(dat$celltype, dat$severity, sep="_")

Idents(dat) <- dat$cg1
tmp <- tryCatch(FindMarkers(dat, ident.1="Basal_SARS-CoV-2", ident.2="Basal_healthy", test.use="MAST", min.pct=0.01, logfc.threshold=0), error=function(e) NULL)
tmp <- cbind(id=row.names(tmp), as.data.table(tmp))
setnames(tmp, c("id","pval","log.fc","pct1","pct2","padj"))
de.res$SC.Chua.Basal <- tmp

# for Ciliated cells, infected vs control failed to run, so used critical vs control instead
Idents(dat) <- dat$cg2
tmp <- tryCatch(FindMarkers(dat, ident.1="Ciliated_critical", ident.2="Ciliated_control", test.use="MAST", min.pct=0.01, logfc.threshold=0), error=function(e) NULL)
tmp <- cbind(id=row.names(tmp), as.data.table(tmp))
setnames(tmp, c("id","pval","log.fc","pct1","pct2","padj"))
de.res$SC.Chua.Ciliated <- tmp



### GSEA analysis

gs1 <- readRDS("../data/c2.cp.kegg.v7.0.symbols.RDS")
gs2 <- readRDS("../data/c2.cp.reactome.v7.0.symbols.RDS")
gsea.res <- lapply(de.res, function(x) rbind(gsea(x, gs1), gsea(x, gs2))[order(padj,pval)])


save(de.res, gsea.res, file="de.and.gsea.res.new.RData")

