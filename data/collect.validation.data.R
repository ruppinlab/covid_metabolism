library(my.utils)
library(GSA)
library(gembox)


##### gene sets

# validated gene sets collected by Kuleshov et al. (Maayan lab)
tmp <- GSA.read.gmt("./kuleshov.genesets.gmt")
genes.maayan <- tmp$genesets
names(genes.maayan) <- tmp$geneset.names
# these represent many different types of gene sets with different biological meanings, so I do not try to take the union or intersection

# info on virus-host PPI
genes.ppi <- list()
# Gordon et al.
genes.ppi[["ppi.gordon"]] <- fread("gordon.ppi.tsv")$genes
# Stukalov et al.
tmp <- fread("stukalov.ppi.tsv", skip=1)
tmp <- tmp[organism=="SARS-CoV-2" & is_contaminant==FALSE & p.adjust(p_value,"BH")<0.05 & !grepl("SARS_CoV",gene_name), unique(gene_name)]
genes.ppi[["ppi.stukalov"]] <- tmp
genes.ppi[["ppi.intersect"]] <- Reduce(intersect, genes.ppi)
genes.ppi[["ppi.union"]] <- unique(unlist(genes.ppi))

# CRISPR-Cas9 or RNAi screen of host factors affecting viral infection
# Wei et al.
tmp <- fread("wei.crispr.tsv")
genes.screen <- list(crispr.wei.pos=tmp[fdr<0.1 & mean_z>0, Gene], crispr.wei.neg=tmp[fdr<0.1 & mean_z<0, Gene])
# Daniloski et al.
tmp <- fread("daniloski.crispr.tsv")
cut1 <- tmp[fdr_MOI1<0.1, min(lfc_MOI1)] # [1] 0.61333
cut3 <- tmp[fdr_MOI3<0.1, min(lfc_MOI3)] # [1] 0.85435
genes.screen <- c(genes.screen,
  list(crispr.daniloski.pos=tmp[fdr_MOI1<0.1 | fdr_MOI3<0.1, gene_id],
       crispr.daniloski.neg=tmp[lfc_MOI1<(-cut1) & lfc_MOI3<(-cut3), gene_id]) # for daniloski.neg, raw p value not give so I filter like this; `|` gave too many neg hits, so used `&`
)
#genes.screen[["pos.crispr.intersect"]] <- Reduce(intersect, genes.screen[c("crispr.wei.pos","crispr.daniloski.pos")]) # only ACE2 and CTSL
genes.screen[["pos.crispr.union"]] <- Reduce(union, genes.screen[c("crispr.wei.pos","crispr.daniloski.pos")])
#genes.screen[["neg.crispr.intersect"]] <- Reduce(intersect, genes.screen[c("crispr.wei.neg","crispr.daniloski.neg")])
genes.screen[["neg.crispr.union"]] <- Reduce(union, genes.screen[c("crispr.wei.neg","crispr.daniloski.neg")])

# also CRISPR-Cas9 full log.fc data (to be used for GSEA)
# Wei et al.
tmp <- fread("genetic_screens/wei.crispr.tsv")
screen.full <- list(wei=tmp[, .(id=Gene, log.fc=mean_lfc)])
# Daniloski et al.
tmp <- fread("genetic_screens/daniloski.crispr.tsv")
screen.full$daniloski.moi1 <- tmp[, .(id=gene_id, log.fc=lfc_MOI1)]
screen.full$daniloski.moi3 <- tmp[, .(id=gene_id, log.fc=lfc_MOI3)]


##### drug sets

tmp <- readRDS(file.path(env("mypr"),"Data/drugbank/5.1.7/drugs2targets.RDS"))
d2t <- tmp[action.simp.strict=="inhibition", .(drug=drug.name, target=protein.gene.symbol)]
tmpf <- function(dat) {
  lapply(dat, function(x) merge(data.table(drug=x), d2t[drug %in% x], by="drug", all=TRUE))
}

# experimental drug sets collected by Maayan lab
tmp <- GSA.read.gmt("kuleshov.experimental.drugsets.gmt")
drugs.exp <- tmp$genesets
names(drugs.exp) <- tmp$geneset.names
drugs.exp <- tmpf(drugs.exp)

# my previously collected, there may be duplications but I simply add it to drugs.exp
tmp <- fread("drug.hits.from.literature.tsv")
tmp <- tmp[note!="MERS-CoV", -"note"]
tmp[study=="", study:="?"]
tmp <- merge(tmp, d2t[drug %in% tmp$drug], by="drug", all=TRUE)
drugs.exp[["my.old.list"]] <- tmp[, .(drug, target, study)]

# drugs from Riva et al.
tmp <- fread("riva.metab.drugs.tsv")
tmp1 <- fread("riva.enriched.drug.targets.tsv")
drugs.exp[["sumit"]] <- rbind(tmp, data.table(drug="unknown", target=tmp1$`Gene Symbol`))

# add union
tmp <- unique(rbindlist(drugs.exp, fill=TRUE)[, .(drug, target)])
drugs.exp[["union"]] <- tmp

# genes mapped to drugs.exp
drugs.exp.genes <- lapply(drugs.exp, function(x) {res <- x[!is.na(target), unique(target)]; if (length(res)==0) NULL else res})
drugs.exp.genes <- drugs.exp.genes[!sapply(drugs.exp.genes, is.null)]

# computational drug sets collected by Maayan lab
tmp <- GSA.read.gmt("kuleshov.computational.drugsets.gmt")
drugs.comp <- tmp$genesets
names(drugs.comp) <- tmp$geneset.names
drugs.comp <- tmpf(drugs.comp)

# add union
tmp <- unique(rbindlist(drugs.comp, fill=TRUE)[, .(drug, target)])
drugs.comp[["union"]] <- tmp

# genes mapped to drugs.comp
drugs.comp.genes <- lapply(drugs.comp, function(x) {res <- x[!is.na(target), unique(target)]; if (length(res)==0) NULL else res})
drugs.comp.genes <- drugs.comp.genes[!sapply(drugs.comp.genes, is.null)]

# map the gene sets to reactios in recon1
tmpf <- function(x, m=recon1) {
  tmp <- unique(unlist(x))
  tmp <- tmp[tmp %in% m$genes]
  mapp <- genes2rxns(m, tmp, out="rxns")
  res <- lapply(x, function(xx) {
    xx <- xx[xx %in% names(mapp)]
    res <- unique(unlist(mapp[xx]))
    if (length(res)==0) NULL else res
  })
  res[!sapply(res, is.null)]
}
rxns.maayan <- tmpf(genes.maayan)
rxns.ppi <- tmpf(genes.ppi)
rxns.screen <- tmpf(genes.screen)
drugs.exp.rxns <- tmpf(drugs.exp.genes)
drugs.comp.rxns <- tmpf(drugs.comp.genes)

# combine data
vgenes <- list(maayan=genes.maayan, ppi=genes.ppi, screen=genes.screen, drugs.exp=drugs.exp.genes, drugs.comp=drugs.comp.genes)
vrxns <- list(maayan=rxns.maayan, ppi=rxns.ppi, screen=rxns.screen, drugs.exp=drugs.exp.rxns, drugs.comp=drugs.comp.rxns)

# besides, for genetic screen datasets, remove intersection between positive set and negative set rxns to be used for generating ROC curve
# note: the negative sets in these datasets were the strongest negative hits; for ROC I decided not to redefine negative sets with a more loose cutoff, as it will result in hugely imbanlanced positive and negative sets
rxns.screen.for.roc <- vrxns$screen
for (x in c("crispr.wei","crispr.daniloski","crispr.union")) {
  if (x=="crispr.union") {
    a <- paste0("pos.", x)
    b <- paste0("neg.", x)
  } else {
    a <- paste0(x, ".pos")
    b <- paste0(x, ".neg")
  }
  tmp <- intersect(rxns.screen.for.roc[[a]], rxns.screen.for.roc[[b]])
  message(length(tmp))
  rxns.screen.for.roc[[a]] <- setdiff(rxns.screen.for.roc[[a]], tmp)
  rxns.screen.for.roc[[b]] <- setdiff(rxns.screen.for.roc[[b]], tmp)
}


save(vgenes, vrxns, screen.full, d2t, drugs.exp, drugs.comp, rxns.screen.for.roc, file="data.for.validation.RData")


