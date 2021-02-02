library(my.utils)
library(DESeq2)

# differential expression and GSEA analysis of Vero E6 data with remdesivir treatment

### DE analysis

mat <- readRDS("../data/remdesivir.RDS")$count

idx <- rowSums(mat>=10)>=6
sum(idx) # [1] 13580
mat1 <- round(mat[idx,])

mat1 <- mat1[, c(1,5:12,2:4)]
colnames(mat1) <- rep(c("ctrl","drug","virus","virus_drug"), each=3)
phe <- data.table(group=factor(rep(c("ctrl","drug","virus","virus_drug"), each=3)))

dds <- DESeqDataSetFromMatrix(countData=mat1, colData=phe, design=~group)
dds <- DESeq(dds)

de.res0 <- list(
  c2v=results(dds, contrast=c("group","virus","ctrl")), # c2v: virus-infected (w/o remdesivir) vs (i.e. compared to) control (non-infected, and w/o remdesivir)
  c2vd=results(dds, contrast=c("group","virus_drug","ctrl")), # c2vd: virus+remdesivir vs control
  c2d=results(dds, contrast=c("group","drug","ctrl")), # c2d: remdesivir-only vs control
  v2vd=results(dds, contrast=c("group","virus_drug","virus")), # v2vd: virus+remdesivir vs virus-infected (w/o remdesivir)
  d2vd=results(dds, contrast=c("group","virus_drug","drug")) # d2vd: virus+remdesivir vs remdesivir-only
)

de.res0 <- lapply(de.res0, function(x) {
  gid <- rownames(x)
  as.data.table(x)[, .(gid=gid, ave=baseMean, log.fc=log2FoldChange, lfc.se=lfcSE, pval=pvalue, padj=padj)][order(padj, pval)]
})

de.res <- lapply(de.res0, function(x) {
  tmp <- str_match(x$gid, "[^_]*_(.*)")[,2]
  idx <- is.na(tmp) | grepl("gp|rRNA|^MIR|Y_RNA|^U[0-9]*$|^SNOR|sno|^RNU|Metazoa|7SK|^SCARNA|RNase|Vault|Telomerase", tmp)
  res <- cbind(data.table(id=tmp[!idx]), x[!idx])
  res[, padj:=p.adjust(pval, "BH")]
  res[order(padj, pval)]
})


### GSEA analysis

gs1 <- readRDS("../data/c2.cp.kegg.v7.0.symbols.RDS")
gs2 <- readRDS("../data/c2.cp.reactome.v7.0.symbols.RDS")
gsea.res <- lapply(de.res, function(x) rbind(gsea(x, gs1), gsea(x, gs2))[order(padj,pval)])


save(de.res, gsea.res, file="de.and.gsea.res.remdesivir.RData")

